#include "align.h"
#include "ksw2.h"
// #include "bindings/cpp/WFAligner.hpp"  // 未使用：当前实现直接使用 C API（wavefront_align/cigar_get_CIGAR）
extern "C" {
#include "alignment/cigar.h"
#include "wavefront/wavefront_align.h"
}

// ==================================================================
// align.cpp - HAlign-4 序列比对算法封装（KSW2 / WFA2）
// ==================================================================
// 本文件实现 align.h 中声明的比对接口，核心目标是：
// 1) 将不同底层比对库（KSW2、WFA2）的调用细节封装起来，对外统一返回 cigar::Cigar_t
// 2) 统一使用项目内的 DNA 字符编码与替换矩阵（ScoreChar2Idx / dna5_simd_mat）
// 3) 明确资源所有权（哪些需要 free/delete），避免内存泄漏
// ==================================================================

namespace align
{
    // ------------------------------------------------------------------
    // globalAlignKSW2：KSW2 全局比对（end-to-end）
    // ------------------------------------------------------------------
    // 输入：ref/query 原始序列（字符串，允许包含 A/C/G/T/N，大小写均可）
    // 输出：cigar::Cigar_t（压缩 CIGAR），描述 query 如何对齐到 ref
    //
    // 实现分 4 步：
    // 1) 序列编码：将字符映射到 0..4（A/C/G/T/N）供 KSW2 使用
    // 2) 配置参数：替换矩阵、gap 罚分、band/zdrop 等
    // 3) 调用 ksw_extz2_sse：获得 ez（包含 CIGAR 指针与长度）
    // 4) 拷贝 CIGAR 并释放 ez.cigar（KSW2 内部 malloc）
    //
    // 正确性/一致性：
    // - 本实现启用 KSW_EZ_GENERIC_SC，确保使用 dna5_simd_mat 完整替换矩阵
    // - 启用 KSW_EZ_RIGHT，让 gap 右对齐，提高 CIGAR 的"规范化"程度，便于后续比较/测试
    // ------------------------------------------------------------------
    cigar::Cigar_t globalAlignKSW2(const std::string& ref, const std::string& query)
    {

        align::KSW2AlignConfig cfg;
        return globalAlignKSW2(ref, query, cfg);
    }


    cigar::Cigar_t globalAlignKSW2(const std::string& ref, const std::string& query, align::KSW2AlignConfig cfg)
    {
        // 边界：任意一条序列为空时，直接返回纯 I / 纯 D 的 CIGAR（避免进入 KSW2）
        if (ref.size() == 0 || query.size() == 0) {
            cigar::Cigar_t cigar; // 返回用的 CIGAR 容器
            if (ref.size() == 0 && query.size() > 0) {
                cigar.push_back(cigar::cigarToInt('I', static_cast<uint32_t>(query.size()))); // query 全插入
            } else if (query.size() == 0 && ref.size() > 0) {
                cigar.push_back(cigar::cigarToInt('D', static_cast<uint32_t>(ref.size())));   // ref 全删除
            }
            return cigar;
        }

        /* ---------- 1. 编码序列 ---------- */
        std::vector<uint8_t> ref_enc(ref.size());   // ref 的 DNA5 编码（A/C/G/T/N -> 0..4）
        std::vector<uint8_t> qry_enc(query.size()); // query 的 DNA5 编码

        for (size_t i = 0; i < ref.size(); ++i)
            ref_enc[i] = align::ScoreChar2Idx[static_cast<uint8_t>(ref[i])]; // 字符->索引
        for (size_t i = 0; i < query.size(); ++i)
            qry_enc[i] = align::ScoreChar2Idx[static_cast<uint8_t>(query[i])]; // 字符->索引

        cfg.band_width = align::auto_band(ref.size(), query.size()); // 启发式带宽（-1 表示不限制）

        /* ---------- 3. 调用 KSW2 ---------- */
        ksw_extz_t ez{}; // KSW2 输出（含 ez.cigar/ez.n_cigar 等；ez.cigar 需要 free）

        // 调用 ksw_extz2_sse：输出写入 ez（注意参数顺序：先 query 再 ref，这是 ksw2 的接口约定）
        ksw_extz2_sse(0,
            static_cast<int>(qry_enc.size()), qry_enc.data(),
            static_cast<int>(ref_enc.size()), ref_enc.data(),
            cfg.alphabet_size, cfg.mat,
            cfg.gap_open, cfg.gap_extend,
            cfg.band_width, cfg.zdrop, cfg.end_bonus,
            cfg.flag, &ez);

        /* ---------- 4. 拷贝 / 释放 CIGAR ---------- */
        cigar::Cigar_t cigar;          // 统一返回格式（压缩 uint32_t 操作序列）
        cigar.reserve(ez.n_cigar);     // 预分配：避免 push_back 扩容
        for (int i = 0; i < ez.n_cigar; ++i)
            cigar.push_back(ez.cigar[i]);  // 拷贝 KSW2 的 CIGAR 单元（uint32_t）

        free(ez.cigar); // 释放：KSW2 内部用 malloc 分配
        return cigar;   // 返回：query->ref 的 CIGAR
    }

    // ------------------------------------------------------------------
    // extendAlignKSW2：KSW2 延伸比对（extension / ends-free 风格）
    // ------------------------------------------------------------------
    // 与 globalAlignKSW2 的差异：
    // - 使用 zdrop 以及 EXTZ_ONLY/APPROX_DROP 等 flag，允许在得分下降时提前终止
    // - 用 end_bonus 奖励末端延伸，有利于从 seed 位置向外扩展
    //
    // 参数：
    // @param zdrop: Z-drop 阈值，越大越不容易提前停止（更慢但更完整）
    // ------------------------------------------------------------------
    cigar::Cigar_t extendAlignKSW2(const std::string& ref,
        const std::string& query,
        int zdrop)
    {
        /* ---------- 1. 序列编码 ---------- */
        std::vector<uint8_t> ref_enc(ref.size());
        std::vector<uint8_t> qry_enc(query.size());
        for (size_t i = 0; i < ref.size(); ++i) ref_enc[i] = align::ScoreChar2Idx[(uint8_t)ref[i]];
        for (size_t i = 0; i < query.size(); ++i) qry_enc[i] = align::ScoreChar2Idx[(uint8_t)query[i]];

        // ------------------------------------------------------------------
        // 2. 配置参数
        // ------------------------------------------------------------------
        // 说明：
        // - 这里保留了历史上尝试过的参数组合（注释块），方便后续调参回归。
        // - 当前启用的 cfg：在 extension 场景下使用 EXTZ_ONLY + APPROX_DROP + RIGHT。
        //   * EXTZ_ONLY：启用 ends-free extension 模式（更适合延伸）
        //   * APPROX_DROP：在 approximate 模式下触发 z-drop 就中断（加速）
        //   * RIGHT：gap 右对齐，CIGAR 更稳定
        // - end_bonus=50：给末端延伸一定奖励，鼓励更长的 alignment。
        // - band_width 启发式估计：加速同时尽量不损失太多准确性。
        align::KSW2AlignConfig cfg;
        cfg.mat = align::dna5_simd_mat;
        cfg.zdrop = zdrop;
        cfg.flag = KSW_EZ_EXTZ_ONLY | KSW_EZ_RIGHT | KSW_EZ_APPROX_DROP;
        cfg.end_bonus = 50;
        cfg.alphabet_size = 5;
        cfg.gap_open = 6;
        cfg.gap_extend = 2;
        cfg.band_width = align::auto_band(ref.size(), query.size());


        /* ---------- 3. 调用 KSW2 ---------- */
        ksw_extz_t ez{};
        // 注意：这里第一个参数传 nullptr 与项目历史实现保持一致。
        ksw_extz2_sse(nullptr,
            static_cast<int>(qry_enc.size()), qry_enc.data(),
            static_cast<int>(ref_enc.size()), ref_enc.data(),
            cfg.alphabet_size, cfg.mat,
            cfg.gap_open, cfg.gap_extend,
            cfg.band_width, cfg.zdrop, cfg.end_bonus,
            cfg.flag, &ez);

        /* ---------- 4. 拷贝 & 释放 ---------- */
        cigar::Cigar_t cigar;
        cigar.reserve(ez.n_cigar);
        for (int i = 0; i < ez.n_cigar; ++i)
            cigar.push_back(ez.cigar[i]);

        free(ez.cigar);                    // ksw2 使用 malloc
        return cigar;                      // 返回的 CIGAR 即延伸片段
    }


    // ------------------------------------------------------------------
    // globalAlignWFA2：WFA2 全局比对
    // ------------------------------------------------------------------
    // WFA2（Wavefront Alignment）特点：
    // - 对高相似度序列（编辑距离小）通常显著快于 DP
    // - 可配置不同 memory_mode（速度/内存权衡）与 heuristic（剪枝策略）
    //
    // 本实现使用 gap_affine（仿射 gap 罚分）：
    // - mismatch=2
    // - gap_opening=3
    // - gap_extension=1
    // 这些参数的量级与 KSW2 的 5/-4 + gap(6,2) 并不完全等价，
    // 但对下游而言我们只需要一个可用且一致的 CIGAR 表示。
    // ------------------------------------------------------------------
    cigar::Cigar_t globalAlignWFA2(const std::string& ref,
        const std::string& query)
    {
        // 1) 构建 WFA2 属性（attributes）
        wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default; // 从默认配置开始（减少遗漏）
        attributes.distance_metric = gap_affine;                              // 选择仿射 gap 代价模型
        attributes.affine_penalties.mismatch = 3;                             // mismatch 罚分（>0）
        attributes.affine_penalties.gap_opening = 4;                          // gap open 罚分（>=0）
        attributes.affine_penalties.gap_extension = 1;                        // gap extend 罚分（>0）

        // memory_mode：当前实现实际选择 high（更快但更占内存）；注释需与赋值一致
        attributes.memory_mode = wavefront_memory_high;

        // heuristic：这里保留参数模板但默认不启用（启用会改变剪枝/速度/最优性权衡）
        // attributes.heuristic.strategy = wf_heuristic_banded_adaptive;
        // attributes.heuristic.min_k = -200;
        // attributes.heuristic.max_k = +200;
        // attributes.heuristic.steps_between_cutoffs = 1;

        // 2) 创建 aligner（堆对象，需要 delete）
        wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);

        // 3) 执行比对
        // wavefront_align 会把结果写入 wf_aligner->cigar
        wavefront_align(wf_aligner, ref.c_str(), ref.length(), query.c_str(), query.length());

        // 4) 从 WFA2 输出中取 CIGAR
        // cigar_get_CIGAR：将内部 cigar 转换为 BAM 风格的 uint32_t buffer
        // - cigar_buffer：由 WFA2 管理（通过 wf_aligner->cigar 释放），这里不手动 free
        // - cigar_length：操作数
        uint32_t* cigar_buffer = nullptr; // 指向内部 buffer（不要 free）
        int cigar_length = 0;             // cigar 操作数量
        cigar_get_CIGAR(wf_aligner->cigar, false, &cigar_buffer, &cigar_length);

        cigar::Cigar_t cigar; // 统一返回格式
        cigar.reserve(static_cast<std::size_t>(cigar_length)); // 预分配
        for (int i = 0; i < cigar_length; ++i)
            cigar.push_back(cigar_buffer[i]); // 拷贝每个 CIGAR 单元

        wavefront_aligner_delete(wf_aligner); // 释放 aligner（同时释放内部 cigar 等资源）

        return cigar; // 返回：query->ref 的 CIGAR
    }

    // ------------------------------------------------------------------
    // extendAlignWFA2：WFA2 延伸比对
    // ------------------------------------------------------------------
    // WFA2 的延伸比对（ends-free extension）：
    // - 通常用于从种子位置向外延伸，快速获得局部比对结果
    // - 可配置 zdrop 阈值，控制延伸过程中的提前终止
    //
    // 本实现使用 wavefront_aligner_attr_default 的高内存模式：
    // - 适合长序列比对，但可能导致较高的内存占用
    // ------------------------------------------------------------------
    // cigar::Cigar_t extendAlignWFA2(const std::string& ref,
    //     const std::string& query, int zdrop)
    // {
    //     wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    //     attributes.distance_metric = gap_affine;
    //     attributes.affine_penalties.mismatch = 2;      // X > 0
    //     attributes.affine_penalties.gap_opening = 3;   // O >= 0
    //     attributes.affine_penalties.gap_extension = 1; // E > 0
    //     attributes.memory_mode = wavefront_memory_high;
    //     attributes.heuristic.strategy = wf_heuristic_zdrop;
    //     attributes.heuristic.zdrop = zdrop;
    //     attributes.heuristic.steps_between_cutoffs = 1;
    //     //// Create a WFAligner
    //     //
    //     wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
    //
    //     wavefront_align(wf_aligner, ref.c_str(), ref.length(), query.c_str(), query.length());
    //     /*wfa::WFAlignerGapAffine aligner(2, 3, 1, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryUltralow);
    //
    //     aligner.alignEnd2End(ref, query);*/
    //
    //     uint32_t* cigar_buffer; // Buffer to hold the resulting CIGAR operations.
    //     int cigar_length = 0; // Length of the CIGAR string.
    //     // Retrieve the CIGAR string from the wavefront aligner.
    //     cigar_get_CIGAR(wf_aligner->cigar, true, &cigar_buffer, &cigar_length);
    //
    //     /* ---------- 4. 拷贝 / 释放 CIGAR ---------- */
    //     cigar::Cigar_t cigar;
    //
    //     for (int i = 0; i < cigar_length; ++i)
    //         cigar.push_back(cigar_buffer[i]);
    //
    //     wavefront_aligner_delete(wf_aligner);
    //
    //     return cigar;
    // }

    // ------------------------------------------------------------------
    // globalAlignMM2：基于锚点（anchors）的分段全局比对（minimap2 风格）
    // ------------------------------------------------------------------
    // 功能概述（精炼）：
    // 使用预先计算的锚点将一对序列（ref / query）拆分为若干可比对的片段，
    // 对每个片段调用指定的全局比对函数（align_func），最后合并各段的 CIGAR
    // 为完整的 query->ref CIGAR。该函数优先保证正确性（CIGAR 消耗长度匹配），
    // 在极端边界上提供稳健的兜底逻辑以避免崩溃或最终不一致性报错。
    //
    // 输入参数：
    // - ref: 参考序列字符串（A/C/G/T/N），长度 ref_len
    // - query: 查询序列字符串，长度 qry_len
    // - anchors: 预先生成的锚点列表（anchor::Anchors），可以是任意顺序，函数内部会排序并链化
    // ------------------------------------------------------------------
    cigar::Cigar_t globalAlignMM2(const std::string& ref,
                                  const std::string& query,
                                  const anchor::Anchors& anchors)
    {
        // 说明：本实现把“锚点链”当作分段边界参考，但推进位置严格由 CIGAR 消耗长度驱动（更稳健）。

        align::KSW2AlignConfig cfg;       // 普通段配置（默认 KSW2 配置）
        align::KSW2AlignConfig first_cfg; // 左端段配置（可用于更保守/不同 flag）
        first_cfg.flag = KSW_EZ_GENERIC_SC; // 显式启用完整替换矩阵（保持与其它路径一致）

        const std::size_t ref_len = ref.size(); // ref 总长度
        const std::size_t qry_len = query.size(); // query 总长度

        // 1) 链化 anchors：得到按得分选择的“最佳链”（可能为空）
        anchor::Anchors sorted_anchors = anchors;                        // 拷贝：避免修改输入
        anchor::ChainParams chain_params = anchor::default_chain_params(); // 默认链化参数
        anchor::Anchors chain_anchors = anchor::chainAnchors(sorted_anchors, chain_params); // 链化
        if (chain_anchors.empty()) {
            return globalAlignKSW2(ref, query); // 无有效锚点：直接全局比对兜底
        }

        // 2) 统一按 query 坐标排序（保证从左到右分段处理）
        std::sort(chain_anchors.begin(), chain_anchors.end(),
                  [](const anchor::Anchor& a, const anchor::Anchor& b) {
                      if (a.pos_qry != b.pos_qry) return a.pos_qry < b.pos_qry;
                      return a.pos_ref < b.pos_ref;
                  });

        cigar::Cigar_t result;                     // 最终 CIGAR
        result.reserve(chain_anchors.size() * 2 + 2); // 估算 reserve：span + gap 段数量

        std::size_t ref_pos = 0; // 当前已消费到的 ref 位置（由 CIGAR 推进）
        std::size_t qry_pos = 0; // 当前已消费到的 query 位置（由 CIGAR 推进）

        auto append_segment = [&](std::size_t ref_start, std::size_t ref_end,
                                  std::size_t qry_start, std::size_t qry_end, align::KSW2AlignConfig seg_cfg) {
            // 边界裁剪：保证 substr 不越界（即使锚点坐标异常也不崩）
            ref_start = std::min(ref_start, ref_len);
            ref_end = std::min(ref_end, ref_len);
            qry_start = std::min(qry_start, qry_len);
            qry_end = std::min(qry_end, qry_len);

            if (ref_end < ref_start) ref_end = ref_start; // 纠正反向区间
            if (qry_end < qry_start) qry_end = qry_start;

            const std::string seg_ref = ref.substr(ref_start, ref_end - ref_start);   // ref 子串
            const std::string seg_qry = query.substr(qry_start, qry_end - qry_start); // query 子串

            cigar::Cigar_t seg_cigar = globalAlignKSW2(seg_ref, seg_qry, seg_cfg); // 本段 end-to-end 全局比对

            // 用 CIGAR 反推消耗长度（关键）：不依赖锚点 span 的“理想长度”
            const std::size_t seg_ref_len = seg_ref.size();
            const std::size_t seg_qry_len = seg_qry.size();
            const std::size_t c_ref = cigar::getRefLength(seg_cigar);   // CIGAR 实际消耗的 ref 长度
            const std::size_t c_qry = cigar::getQueryLength(seg_cigar); // CIGAR 实际消耗的 query 长度

            if (c_ref != seg_ref_len || c_qry != seg_qry_len) {
#ifdef _DEBUG
                spdlog::warn("globalAlignMM2(seg): segment cigar mismatch (expected ref:{}/qry:{}, got ref:{}/qry:{}); forcing robust fallback for this segment",
                             seg_ref_len, seg_qry_len, c_ref, c_qry);
#endif
                // 兜底策略：强制让 Query 先全部 I、Ref 再全部 D，保证一定回到对角线
                cigar::Cigar_t forced_cigar;
                if (seg_qry_len > 0) {
                    forced_cigar.push_back(cigar::cigarToInt('I', static_cast<uint32_t>(seg_qry_len)));
                }
                if (seg_ref_len > 0) {
                    forced_cigar.push_back(cigar::cigarToInt('D', static_cast<uint32_t>(seg_ref_len)));
                }
                cigar::appendCigar(result, forced_cigar); // 合并到总 CIGAR

                ref_pos = ref_end; // 指针按“预期区间”推进（因为 forced 已覆盖全段）
                qry_pos = qry_end;
                return;
            }

            cigar::appendCigar(result, seg_cigar); // 合并本段 CIGAR（会自动合并相邻同 op）

            ref_pos = ref_start + c_ref; // 指针按 CIGAR 消耗推进
            qry_pos = qry_start + c_qry;
        };

        // 3) 左端：序列起点到第一个锚点开始（可能为空段）
        {
            const auto& first = chain_anchors.front();
            append_segment(ref_pos, first.pos_ref, qry_pos, first.pos_qry, first_cfg);
        }

        // 4) 逐锚点：先处理锚点 span，再处理到下一个锚点的 gap
        for (std::size_t i = 0; i < chain_anchors.size(); ++i) {
            const auto& a = chain_anchors[i];

            const std::size_t a_ref_start = static_cast<std::size_t>(a.pos_ref);           // 锚点 ref 起点
            const std::size_t a_qry_start = static_cast<std::size_t>(a.pos_qry);           // 锚点 qry 起点
            const std::size_t a_ref_end = a_ref_start + static_cast<std::size_t>(a.span); // 锚点 ref 终点
            const std::size_t a_qry_end = a_qry_start + static_cast<std::size_t>(a.span); // 锚点 qry 终点

            append_segment(ref_pos, a_ref_end, qry_pos, a_qry_end, cfg); // span 段

            if (i + 1 < chain_anchors.size()) {
                const auto& b = chain_anchors[i + 1];
                append_segment(ref_pos, b.pos_ref, qry_pos, b.pos_qry, cfg); // gap 段
            }
        }

        // 5) 右端：最后一个锚点之后到序列末尾
        append_segment(ref_pos, ref_len, qry_pos, qry_len, cfg);

        // 6) 最终一致性检查：如不匹配则退化为全局比对（保证输出合法）
        const std::size_t total_ref = cigar::getRefLength(result);
        const std::size_t total_qry = cigar::getQueryLength(result);
        if (total_ref != ref_len || total_qry != qry_len) {
            spdlog::error("globalAlignMM2: final cigar mismatch (ref:{}/{}, qry:{}/{}), fallback to global",
                         total_ref, ref_len, total_qry, qry_len);
            return globalAlignKSW2(ref, query);
        }

        return result; // 返回最终拼接后的 CIGAR
    }

} // namespace align
