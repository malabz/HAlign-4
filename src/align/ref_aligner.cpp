#include "align.h"
#include "config.hpp"
#include <algorithm>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>

namespace align {

    // ------------------------------------------------------------------
    // 构造函数1：直接传入参数初始化
    // 说明：
    // 1. 如果提供了 consensus_string，则直接使用它作为 reference
    // 2. 否则从 ref_fasta_path 读取参考序列
    // 3. consensus_string 会被保存到成员变量 consensus_seq 中
    // 4. keep_first_length 和 keep_all_length 会被保存到成员变量
    // ------------------------------------------------------------------
    RefAligner::RefAligner(const FilePath& work_dir, const FilePath& ref_fasta_path, int kmer_size, int window_size,
                           int sketch_size, bool noncanonical, std::string consensus_string,
                           bool keep_first_length, bool keep_all_length)
        : work_dir(work_dir),
          kmer_size(kmer_size),
          window_size(window_size),
          sketch_size(sketch_size),
          noncanonical(noncanonical),
          consensus_seq(std::move(consensus_string)),  // 保存 consensus_string 到成员变量
          keep_first_length(keep_first_length),        // 保存 keep_first_length 到成员变量
          keep_all_length(keep_all_length)             // 保存 keep_all_length 到成员变量
    {
        // 从文件读取参考序列并构建索引
        seq_io::KseqReader reader(ref_fasta_path);
        seq_io::SeqRecord rec;
        while (reader.next(rec))
        {
            // 关键修复：必须在 move(rec) 之前计算 sketch 和 minimizer
            auto sketch = mash::sketchFromSequence(rec.seq, kmer_size, sketch_size, noncanonical, random_seed);
            auto minimizer = minimizer::extractMinimizerHash(rec.seq, kmer_size, window_size, noncanonical);

            ref_sequences.push_back(std::move(rec));
            ref_sketch.push_back(std::move(sketch));
            ref_minimizers.push_back(std::move(minimizer));
        }
    }

    // ------------------------------------------------------------------
    // 构造函数2：基于 Options 结构体初始化
    // 说明：
    // 1. 从 Options 中提取相关参数，委托给第一个构造函数
    // 2. 将 consensus_string 传递给第一个构造函数并保存到成员变量
    // 3. 从 opt 中提取 keep_first_length 和 keep_all_length
    // ------------------------------------------------------------------
    RefAligner::RefAligner(const Options& opt, const FilePath& ref_fasta_path, std::string consensus_string)
        : RefAligner(
            opt.workdir,                // work_dir：工作目录
            ref_fasta_path,             // 参考序列文件路径
            opt.kmer_size,              // kmer_size：k-mer 大小
            opt.kmer_window,            // window_size：minimizer 窗口大小
            opt.sketch_size,            // sketch_size：sketch 大小
            true,                       // noncanonical：是否使用非标准模式（固定为 true）
            std::move(consensus_string),// consensus_string：共识序列（移动语义）
            opt.keep_first_length,      // keep_first_length：从 opt 提取
            opt.keep_all_length         // keep_all_length：从 opt 提取
        )
    {
        // 委托构造函数已完成所有初始化工作
    }

    void RefAligner::alignOneQueryToRef(const seq_io::SeqRecord& q, int tid) const
    {

        auto& out = *outs[tid];
        auto& out_insertion = *outs_with_insertion[static_cast<std::size_t>(tid)];
        // 1) 计算 query sketch
        const mash::Sketch qsk = mash::sketchFromSequence(
            q.seq,
            static_cast<std::size_t>(kmer_size),
            static_cast<std::size_t>(sketch_size),
            noncanonical,
            random_seed);

        // 2) 选择最相似 reference（线性扫描）
        double best_j = -1.0;
        std::size_t best_r = 0;

        for (std::size_t r = 0; r < ref_sketch.size(); ++r) {
            const double j = mash::jaccard(qsk, ref_sketch[r]);



            if (j > best_j) {
                best_j = j;
                best_r = r;
            }
        }



        const auto& best_ref = ref_sequences[best_r];

        // 3) 执行全局比对（默认使用 WFA2，可切换为 KSW2）
        cigar::Cigar_t cigar;
        if (true)
        {
            cigar = globalAlignWFA2(best_ref.seq, q.seq);
        }
        else
        {
            cigar= globalAlignKSW2(best_ref.seq, q.seq);
        }

        // 4) 将 CIGAR 从压缩格式转为 SAM 字符串格式（例如 "100M5I95M"）
        const std::string cigar_str = cigar::cigarToString(cigar);

        // makeSamRecord获得sam_rec
        const auto sam_rec = seq_io::makeSamRecord(
            q,
            best_ref.id,
            cigar_str,
            1,      // pos，简化处理为 1
            60,     // mapq，简化处理为 60
            0       // flag，简化处理为 0（正向链、未比对）
        );


        // 6) 根据是否存在插入，决定写入哪个输出文件
        // 说明：
        // - 如果 CIGAR 中存在插入操作（'I'），需要进行二次判断
        // - 二次判断依据：根据 keep_first_length 标志，选择不同的参考序列进行比对
        //   * keep_first_length = true:  与 ref_seq[0]（第一条参考序列）比对
        //   * keep_first_length = false: 与 consensus_seq（共识序列）比对
        // - 比对后，如果新的 CIGAR 中仍然有插入，则写入 out_insertion；否则写入 out
        if (cigar::hasInsertion(cigar)) {
            cigar::Cigar_t recheck_cigar;

            if (keep_first_length) {
                // 策略1：与第一条参考序列（ref_seq[0]）进行比对
                // 说明：当需要保持第一条序列长度不变时，以它为基准判断插入
                if (!ref_sequences.empty()) {
                    recheck_cigar = globalAlignWFA2(ref_sequences[0].seq, q.seq);
                } else {
                    // 如果没有参考序列，回退到原始 CIGAR
                    recheck_cigar = cigar;
                }
            } else {
                // 策略2：与共识序列（consensus_seq）进行比对
                // 说明：当使用共识序列作为参考时，以它为基准判断插入
                if (!consensus_seq.empty()) {
                    recheck_cigar = globalAlignWFA2(consensus_seq, q.seq);
                } else {
                    // 如果没有共识序列，回退到原始 CIGAR
                    recheck_cigar = cigar;
                }
            }

            // 重新检查二次比对的 CIGAR 是否有插入
            // 如果仍然有插入，说明这是真正的插入变异，写入 insertion 文件
            // 如果没有插入，说明可能是比对策略导致的差异，写入普通文件
            if (cigar::hasInsertion(recheck_cigar)) {
                // 重新生成 SAM 记录（使用二次比对的 CIGAR）
                const std::string recheck_cigar_str = cigar::cigarToString(recheck_cigar);
                const auto recheck_sam_rec = seq_io::makeSamRecord(
                    q,
                    keep_first_length ? ref_sequences[0].id : "consensus",
                    recheck_cigar_str,
                    1,
                    60,
                    0
                );
                out_insertion.writeSam(recheck_sam_rec);
            } else {
                // 二次比对后无插入，写入普通文件（使用二次比对的 CIGAR）
                const std::string recheck_cigar_str = cigar::cigarToString(recheck_cigar);
                const auto recheck_sam_rec = seq_io::makeSamRecord(
                    q,
                    keep_first_length ? ref_sequences[0].id : "consensus",
                    recheck_cigar_str,
                    1,
                    60,
                    0
                );
                out.writeSam(recheck_sam_rec);
            }
        } else {
            // 没有插入，直接写入普通输出文件
            out.writeSam(sam_rec);
        }

    }

    void RefAligner::alignQueryToRef(const FilePath& qry_fasta_path, int threads, std::size_t batch_size)
    {
        // =====================================================================
        // OpenMP 并行骨架（流式读取 + 每线程独立 writer）
        //
        // 核心思路：
        // - 单线程按 chunk 流式读取 query；
        // - 对一个 chunk 用 OpenMP 并行 for；
        // - 每线程有稳定 thread_id (= omp_get_thread_num)；
        // - 每线程一个独占输出文件，避免加锁。
        //
        // 我们把“相似度计算 +（占位）比对 + 写出”抽象成一个函数 process_one()。
        // 这样后续你只需要替换 process_one() 的内部逻辑即可，不影响并行框架。
        // =====================================================================

        if (ref_sequences.empty() || ref_sketch.empty()) {
            throw std::runtime_error("RefAligner::alignQueryToRef: reference is empty");
        }

        // batch_size==0 没有意义，这里兜底为 1
        if (batch_size == 0) batch_size = 1;

        // 线程数：
        // - threads<=0：尊重 OpenMP 默认（OMP_NUM_THREADS 或 runtime 配置）
        // - threads>0 ：显式设置
        if (threads > 0) {
            omp_set_num_threads(threads);
        }
        const int nthreads = std::max(1, omp_get_max_threads());

        FilePath result_dir = work_dir / RESULTS_DIR;
        file_io::ensureDirectoryExists(result_dir, "result directory");

        // 每线程一个输出
        // writer（SAM 模式，占位输出也能保持格式合法）
        outs.clear();
        outs_with_insertion.clear();
        outs.resize(static_cast<std::size_t>(nthreads));
        outs_with_insertion.resize(static_cast<std::size_t>(nthreads));  // 关键修复：必须 resize 以避免越界访问

        for (int tid = 0; tid < nthreads; ++tid) {
            FilePath out_path = result_dir / ("thread" + std::to_string(tid) + ".sam");
            FilePath out_path_insertion = result_dir / ("thread" + std::to_string(tid) + "_insertion.sam");

            auto tmp = seq_io::SeqWriter::Sam(out_path);
            auto tmp_insertion = seq_io::SeqWriter::Sam(out_path_insertion);
            outs[static_cast<std::size_t>(tid)] =
                std::make_unique<seq_io::SeqWriter>(std::move(tmp));

            outs[static_cast<std::size_t>(tid)]->writeSamHeader("@HD\tVN:1.6\tSO:unknown");
            outs_with_insertion[static_cast<std::size_t>(tid)] =
                std::make_unique<seq_io::SeqWriter>(std::move(tmp_insertion));
            outs_with_insertion[static_cast<std::size_t>(tid)]->writeSamHeader("@HD\tVN:1.6\tSO:unknown");
        }

        // ---- 流式读取（单线程）+ chunk/batch 并行处理 ----
        seq_io::KseqReader reader(qry_fasta_path);
        std::vector<seq_io::SeqRecord> chunk;
        chunk.reserve(batch_size);

        while (true)
        {
            chunk.clear();
            seq_io::SeqRecord rec;
            for (std::size_t i = 0; i < batch_size; ++i) {
                if (!reader.next(rec)) break;
                chunk.push_back(std::move(rec));
            }
            if (chunk.empty()) break;

            // OpenMP 规约：
            // - default(none) 强制显式标注共享/私有变量，避免无意的数据竞争；
            // - outs 与 chunk 是跨线程共享只读/线程索引访问的；tid/out/i 为线程私有。
            #pragma omp parallel default(none) shared(outs, chunk)
            {
                const int tid = omp_get_thread_num();

                // 这里每个 query 的处理耗时可能不均匀（不同长度/不同参考命中），用 dynamic(1) 先保证负载均衡。
                // 若后续确认为均匀负载，可调整为 guided/static 并基准测试。
                #pragma omp for schedule(dynamic, 1)
                for (std::int64_t i = 0; i < static_cast<std::int64_t>(chunk.size()); ++i)
                {
                    alignOneQueryToRef(chunk[static_cast<std::size_t>(i)], tid);
                }
            }

            for (auto& w : outs) {
                w->flush();
            }
            for (auto& w : outs_with_insertion) {
                w->flush();
            }
        }
    }

} // namespace align

