#include "align.h"
#include "ksw2.h"
#include "bindings/cpp/WFAligner.hpp"
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
//
// 重要约定（正确性）：
// - 所有接口都返回"query 相对于 ref"的 CIGAR（即描述 query 如何对齐到 ref）。
// - CIGAR 的压缩编码与 BAM/SAM 规范一致：cigar_unit=(len<<4)|op。
// - 本文件只做"调用封装"与"参数配置"，不改变算法本身的比对逻辑。
//
// 性能说明：
// - KSW2：经典 DP（SIMD 加速），对中短序列稳定；对长序列可用 band 降低复杂度。
// - WFA2：Wavefront，编辑距离小（高相似度）时非常快；相似度低时可能退化。
// - 编码阶段会产生两个临时 vector（ref_enc/qry_enc），属于必要开销；
//   若未来成为热点，可考虑复用缓冲（但会引入状态/线程安全复杂度，本次不改行为）。
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
        /* ---------- 1. 编码序列 ---------- */
        // 说明：KSW2 期望输入为整数编码序列。
        // ScoreChar2Idx 的映射规则见 align.h：A/C/G/T -> 0..3，其它 -> 4(N)
        std::vector<uint8_t> ref_enc(ref.size());
        std::vector<uint8_t> qry_enc(query.size());

        for (size_t i = 0; i < ref.size(); ++i)
            ref_enc[i] = align::ScoreChar2Idx[static_cast<uint8_t>(ref[i])];
        for (size_t i = 0; i < query.size(); ++i)
            qry_enc[i] = align::ScoreChar2Idx[static_cast<uint8_t>(query[i])];

        /* ---------- 2. 配置 KSW2 参数 ---------- */
        // dna5_simd_mat 是 constexpr 编译期常量，无需运行时初始化。
        // 这里 cfg 的所有字段都显式填写，避免依赖未初始化内存。
        align::KSW2AlignConfig cfg;
        cfg.mat = align::dna5_simd_mat;
        cfg.alphabet_size = 5;

        // Gap 惩罚参数（KSW2 内部会转为负数，这里传正值）
        // - gap_open=6：开启一个 gap 的基础代价
        // - gap_extend=2：每延长 1bp 的额外代价
        // - 总代价公式：gap of length L 的代价 = -(gap_open + gap_extend * L)
        // - 例如：1bp gap 代价 = -(6+2) = -8，与单个错配代价 -9 接近
        cfg.gap_open = 6;
        cfg.gap_extend = 2;

        // end_bonus：全局比对通常不需要 ends-free 奖励
        cfg.end_bonus = 0;

        // zdrop：全局比对应完整跑完，一般不依赖 z-drop 提前终止
        // 这里传 -1 表示禁用（与项目原有逻辑保持一致）
        cfg.zdrop = -1;

        // band_width：限制 DP 计算在对角线附近的带宽，可显著加速长序列。
        // 注意：带宽过小会导致比对失败或产生次优路径。
        // 这里用 auto_band 做启发式估计，属于"性能/准确性折中"的经验值。
        cfg.band_width = align::auto_band(ref.size(), query.size());

        // Flag 组合说明：
        // - KSW_EZ_GENERIC_SC：启用完整替换矩阵（必须！否则会忽略 dna5_simd_mat）
        // - KSW_EZ_RIGHT：gap 右对齐（CIGAR 标准化，便于比较）
        cfg.flag = KSW_EZ_GENERIC_SC | KSW_EZ_RIGHT;

        /* ---------- 3. 调用 KSW2 ---------- */
        // ez：KSW2 输出结构体。
        // 关键字段：
        // - ez.cigar：uint32_t*，malloc 分配，需要调用者 free
        // - ez.n_cigar：CIGAR 操作数量
        ksw_extz_t ez{};

        // 注意：ksw_extz2_sse 的第一个参数通常是 qlen/tlen 的额外参数（如 score-only 模式），
        // 这里传 0 与现有逻辑保持一致。
        ksw_extz2_sse(0,
            static_cast<int>(qry_enc.size()), qry_enc.data(),
            static_cast<int>(ref_enc.size()), ref_enc.data(),
            cfg.alphabet_size, cfg.mat,
            cfg.gap_open, cfg.gap_extend,
            cfg.band_width, cfg.zdrop, cfg.end_bonus,
            cfg.flag, &ez);


        /* ---------- 4. 拷贝 / 释放 CIGAR ---------- */
        // 输出转换：把 KSW2 的 uint32_t CIGAR 拷贝到项目统一的 cigar::Cigar_t
        // 性能：reserve(ez.n_cigar) 避免 push_back 扩容
        cigar::Cigar_t cigar;
        cigar.reserve(ez.n_cigar);
        for (int i = 0; i < ez.n_cigar; ++i)
            cigar.push_back(ez.cigar[i]);

        // 资源释放：KSW2 用 malloc() 分配 ez.cigar -> 需要 free()
        free(ez.cigar);
        return cigar;
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
        // 注意：wavefront_aligner_attr_default 会填充大量默认值。
        // 本函数只覆盖我们关心的部分，保持行为稳定。
        wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
        attributes.distance_metric = gap_affine;
        attributes.affine_penalties.mismatch = 2;      // X > 0
        attributes.affine_penalties.gap_opening = 3;   // O >= 0
        attributes.affine_penalties.gap_extension = 1; // E > 0

        // memory_mode：ultralow 表示尽量降低内存占用，适合长序列但可能更慢
        attributes.memory_mode = wavefront_memory_ultralow;

        // heuristic：自适应 band 策略，限制波前带宽以加速
        // 注意：heuristic 的启发式会影响速度与最优性，本项目当前选择偏速度。
        attributes.heuristic.strategy = wf_heuristic_banded_adaptive;
        attributes.heuristic.min_k = -50;
        attributes.heuristic.max_k = +50;
        attributes.heuristic.steps_between_cutoffs = 1;

        // 2) 创建 aligner（堆对象，需要 delete）
        wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);

        // 3) 执行比对
        // wavefront_align 会把结果写入 wf_aligner->cigar
        wavefront_align(wf_aligner, ref.c_str(), ref.length(), query.c_str(), query.length());

        // 4) 从 WFA2 输出中取 CIGAR
        // cigar_get_CIGAR：将内部 cigar 转换为 BAM 风格的 uint32_t buffer
        // - cigar_buffer：由 WFA2 管理（通过 wf_aligner->cigar 释放），这里不手动 free
        // - cigar_length：操作数
        uint32_t* cigar_buffer;
        int cigar_length = 0;
        cigar_get_CIGAR(wf_aligner->cigar, false, &cigar_buffer, &cigar_length);

        // 5) 拷贝到项目统一格式
        cigar::Cigar_t cigar;
        cigar.reserve(static_cast<std::size_t>(cigar_length));
        for (int i = 0; i < cigar_length; ++i)
            cigar.push_back(cigar_buffer[i]);

        // 6) 释放 aligner（会同时释放内部 cigar 结构）
        wavefront_aligner_delete(wf_aligner);

        return cigar;
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
}