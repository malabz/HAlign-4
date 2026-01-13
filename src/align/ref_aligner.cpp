#include "align.h"
#include "config.hpp"
#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <omp.h>

namespace align {

    RefAligner::RefAligner(const FilePath& work_dir, const FilePath& ref_fasta_path, int kmer_size, int window_size, int sketch_size, bool noncanonical)
        : work_dir(work_dir), kmer_size(kmer_size), window_size(window_size), sketch_size(sketch_size), noncanonical(noncanonical)
    {
        seq_io::KseqReader reader(ref_fasta_path);
        seq_io::SeqRecord rec;
        while (reader.next(rec))
        {
            ref_sequences.push_back(std::move(rec));
            ref_sketch.push_back(mash::sketchFromSequence(rec.seq, kmer_size, sketch_size, noncanonical, random_seed));
            ref_minimizers.push_back(minimizer::extractMinimizerHash(rec.seq, kmer_size, window_size, noncanonical));
        }
    }

    void RefAligner::alignOneQueryToRef(const seq_io::SeqRecord& q, seq_io::SeqWriter& out) const
    {

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

        // 3) TODO: 这里做真正的“比对”
        // - minimizer chaining / DP / 外部 aligner
        // - 产出 POS / CIGAR / MAPQ 等字段

        // 4) 写出（占位 SAM）：
        // - QNAME=query id
        // - RNAME=best ref id
        // - OPT=写入一些调试 tag（jaccard、ref index、长度）
        seq_io::SeqWriter::SamRecord sr;
        sr.qname = q.id;
        sr.flag = 0;
        sr.rname = best_ref.id;
        sr.pos = 0;
        sr.mapq = 0;
        sr.cigar = "*";
        sr.rnext = "*";
        sr.pnext = 0;
        sr.tlen = 0;
        sr.seq = "*";
        sr.qual = "*";

        std::ostringstream opt;
        opt << "JI:f:" << std::fixed << std::setprecision(6) << best_j
            << "\tRI:i:" << best_r
            << "\tQL:i:" << q.seq.size()
            << "\tRL:i:" << best_ref.seq.size();
        const std::string opt_s = opt.str();
        sr.opt = opt_s;

        out.writeSam(sr);
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
        std::vector<std::unique_ptr<seq_io::SeqWriter>> outs;
        outs.resize(static_cast<std::size_t>(nthreads));
        for (int tid = 0; tid < nthreads; ++tid) {
            FilePath out_path = result_dir / ("thread" + std::to_string(tid) + ".sam");

            auto tmp = seq_io::SeqWriter::Sam(out_path);
            outs[static_cast<std::size_t>(tid)] =
                std::make_unique<seq_io::SeqWriter>(std::move(tmp));

            outs[static_cast<std::size_t>(tid)]->writeSamHeader("@HD\tVN:1.6\tSO:unknown");
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
                auto& out = *outs[static_cast<std::size_t>(tid)];

                // 这里每个 query 的处理耗时可能不均匀（不同长度/不同参考命中），用 dynamic(1) 先保证负载均衡。
                // 若后续确认为均匀负载，可调整为 guided/static 并基准测试。
                #pragma omp for schedule(dynamic, 1)
                for (std::int64_t i = 0; i < static_cast<std::int64_t>(chunk.size()); ++i)
                {
                    alignOneQueryToRef(chunk[static_cast<std::size_t>(i)], out);
                }
            }

            for (auto& w : outs) {
                w->flush();
            }
        }
    }

} // namespace align

