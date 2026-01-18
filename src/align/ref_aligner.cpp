#include "align.h"
#include "config.hpp"
#include "preprocess.h"    // alignConsensusSequence
#include "consensus.h"     // generateConsensusSequence
#include <algorithm>
#include <cstddef>
#include <fstream>         // std::ifstream, std::ofstream
#include <filesystem>      // std::filesystem::remove
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>

namespace align {

    // ------------------------------------------------------------------
    // 构造函数1：直接传入参数初始化
    // 说明：
    // 1. 如果 keep_first_length == true，使用 ref_sequences[0] 作为参考
    // 2. 否则，调用 MSA 生成共识序列作为参考
    // 3. threads 和 msa_cmd 参数用于共识序列生成
    // 4. keep_first_length 和 keep_all_length 会被保存到成员变量
    // ------------------------------------------------------------------
    RefAligner::RefAligner(const FilePath& work_dir, const FilePath& ref_fasta_path, int kmer_size, int window_size,
                           int sketch_size, bool noncanonical, int threads, std::string msa_cmd,
                           bool keep_first_length, bool keep_all_length)
        : work_dir(work_dir),
          kmer_size(kmer_size),
          window_size(window_size),
          sketch_size(sketch_size),
          noncanonical(noncanonical),
          threads(threads),
          msa_cmd(std::move(msa_cmd)),          // 使用 move 语义避免拷贝
          keep_first_length(keep_first_length),
          keep_all_length(keep_all_length)
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

        if (keep_first_length)
        {
            consensus_seq = ref_sequences.front();
        }
        else
        {
            FilePath consensus_unaligned_file = ref_fasta_path;
            FilePath consensus_aligned_file = FilePath(work_dir) / WORKDIR_DATA / DATA_CLEAN/ CLEAN_CONS_ALIGNED;
            FilePath consensus_file = FilePath(work_dir) / WORKDIR_DATA / DATA_CLEAN/ CLEAN_CONS_FASTA;
            FilePath consensus_json_file = FilePath(work_dir) / WORKDIR_DATA / DATA_CLEAN/ CLEAN_CONS_JSON;

            const std::size_t batch_size = 4096;
            alignConsensusSequence(consensus_unaligned_file, consensus_aligned_file, msa_cmd, threads);
            std::string consensus_string = consensus::generateConsensusSequence(
                consensus_aligned_file,
                consensus_file,
                consensus_json_file,
                0, // 不限制数量
                threads,
                batch_size
            );
            consensus_seq.id = "consensus";
            consensus_seq.seq = std::move(consensus_string);
        }


    }

    // ------------------------------------------------------------------
    // 构造函数2：基于 Options 结构体初始化
    // 说明：
    // 1. 从 Options 中提取相关参数，委托给第一个构造函数
    // 2. 构造函数1会根据 keep_first_length 标志自动生成或选择参考序列
    // 3. 从 opt 中提取 keep_first_length 和 keep_all_length
    // ------------------------------------------------------------------
    RefAligner::RefAligner(const Options& opt, const FilePath& ref_fasta_path)
        : RefAligner(
            opt.workdir,                // work_dir：工作目录
            ref_fasta_path,             // 参考序列文件路径
            opt.kmer_size,              // kmer_size：k-mer 大小
            opt.kmer_window,            // window_size：minimizer 窗口大小
            opt.sketch_size,            // sketch_size：sketch 大小
            true,                       // noncanonical：是否使用非标准模式（固定为 true）
            opt.threads,                // threads：线程数（用于共识序列生成）
            opt.msa_cmd,                // msa_cmd：MSA 命令模板
            opt.keep_first_length,      // keep_first_length：从 opt 提取
            opt.keep_all_length         // keep_all_length：从 opt 提取
        )
    {
        // 委托构造函数已完成所有初始化工作（包括共识序列生成）
    }


    // ------------------------------------------------------------------
    // 辅助函数：写入SAM记录
    // ------------------------------------------------------------------
    void RefAligner::writeSamRecord(const seq_io::SeqRecord& q, const cigar::Cigar_t& cigar,
                                    std::string_view ref_name, seq_io::SeqWriter& out) const
    {
        const std::string cigar_str = cigar::cigarToString(cigar);
        const auto sam_rec = seq_io::makeSamRecord(q, ref_name, cigar_str, 1, 60, 0);
        out.writeSam(sam_rec);
    }

    // ------------------------------------------------------------------
    // 辅助函数：mergeConsensusAndSamToFasta
    // 功能：将共识序列和多个 SAM 文件合并为一个 FASTA 文件
    //
    // 实现说明：
    // 1. 创建一个 SeqWriter 打开输出 FASTA 文件
    // 2. 先写入共识序列（consensus_seq）
    // 3. 逐个读取 SAM 文件，提取序列并追加写入到同一个 writer
    // 4. 使用同一个 writer 实例，避免文件追加模式的复杂性和性能损失
    // 5. 所有写入操作在同一个 writer 的生命周期内完成
    //
    // 性能优化：
    // - 使用大缓冲区（默认 8MiB）批量写入，减少系统调用
    // - 流式处理 SAM 文件，内存占用与文件大小无关
    // - 避免文件的重复打开和关闭
    // ------------------------------------------------------------------
    std::size_t RefAligner::mergeConsensusAndSamToFasta(
        const std::vector<FilePath>& sam_paths,
        const FilePath& fasta_path,
        std::size_t line_width) const
    {
        // 创建 FASTA writer（非追加模式，会覆盖已存在的文件）
        seq_io::SeqWriter writer(fasta_path, line_width);

        // 1. 先写入共识序列
        writer.writeFasta(consensus_seq);
        writer.flush();

        // 2. 统计信息（用于返回值和调试）
        std::size_t total_count = 1;  // 已写入共识序列，计数从 1 开始
        std::size_t file_idx = 0;

        // 3. 逐个处理每个 SAM 文件
        for (const auto& sam_path : sam_paths) {
            // 打开当前 SAM 文件
            // 说明：
            // - 每个 SamReader 独立打开和关闭
            // - 使用 RAII 确保文件在处理完毕后自动关闭
            // - 如果文件打开失败，会抛出异常并中止合并
            seq_io::SamReader reader(sam_path);

            // 逐条读取并追加写入
            seq_io::SeqRecord rec;
            std::size_t file_count = 0;

            while (reader.next(rec)) {
                // 写入 FASTA 格式（只保留 id 和 seq，丢弃 qual）
                writer.writeFasta(rec);
                ++file_count;
                ++total_count;
            }

            // 调试信息：记录每个文件的处理进度
            #ifdef _DEBUG
            spdlog::debug("mergeConsensusAndSamToFasta: processed file {} ({}/{}): {} records from {}",
                         file_idx + 1, file_idx + 1, sam_paths.size(), file_count, sam_path.string());
            #endif

            ++file_idx;
        }

        // 4. 确保所有数据已刷新到磁盘
        writer.flush();

        // 调试信息：记录合并统计
        #ifdef _DEBUG
        spdlog::debug("mergeConsensusAndSamToFasta: merged {} SAM files ({} total records) to {}",
                     sam_paths.size(), total_count, fasta_path.string());
        #endif

        return total_count;
    }

    // ------------------------------------------------------------------
    // 主函数：alignOneQueryToRef
    // 功能：对单个query执行比对并写入SAM文件
    //
    // 流程：
    // 1. 计算 query sketch
    // 2. 找到最相似的 reference
    // 3. 执行全局比对
    // 4. 根据是否有插入进行二次判断并写入相应文件
    //
    // 参数说明：
    //   - q: 待比对的查询序列
    //   - out: 普通输出文件的 writer（无插入或二次比对后无插入的序列）
    //   - out_insertion: 插入序列输出文件的 writer（二次比对后仍有插入的序列）
    // ------------------------------------------------------------------
    void RefAligner::alignOneQueryToRef(const seq_io::SeqRecord& q,
                                       seq_io::SeqWriter& out,
                                       seq_io::SeqWriter& out_insertion) const
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

        // 3) 执行全局比对（使用 WFA2）
        const cigar::Cigar_t initial_cigar = globalAlignWFA2(best_ref.seq, q.seq);

        // 4) 根据是否存在插入，决定写入哪个输出文件
        if (!cigar::hasInsertion(initial_cigar)) {
            // 无插入：直接写入普通文件
            writeSamRecord(q, initial_cigar, best_ref.id, out);
            return;
        }

        // 有插入：进行二次比对判断
        const cigar::Cigar_t recheck_cigar = globalAlignWFA2(consensus_seq.seq, q.seq);;

        // 使用二次比对结果（如果失败则使用初始结果）
        const cigar::Cigar_t& final_cigar = recheck_cigar.empty() ? initial_cigar : recheck_cigar;

        // 根据二次比对结果决定输出文件
        if (cigar::hasInsertion(final_cigar)) {
            writeSamRecord(q, final_cigar, consensus_seq.id, out_insertion);
        } else {
            writeSamRecord(q, final_cigar, consensus_seq.id, out);
        }

    }

    void RefAligner::alignQueryToRef(const FilePath& qry_fasta_path, std::size_t batch_size)
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
        outs_path.clear();
        outs_path.resize(static_cast<std::size_t>(nthreads));
        outs_with_insertion_path.clear();
        outs_with_insertion_path.resize(static_cast<std::size_t>(nthreads));


        std::vector<std::unique_ptr<seq_io::SeqWriter>> outs;
        std::vector<std::unique_ptr<seq_io::SeqWriter>> outs_with_insertion;

        outs.clear();
        outs_with_insertion.clear();
        outs.resize(static_cast<std::size_t>(nthreads));
        outs_with_insertion.resize(static_cast<std::size_t>(nthreads));  // 关键修复：必须 resize 以避免越界访问

        for (int tid = 0; tid < nthreads; ++tid) {
            FilePath out_path = result_dir / ("thread" + std::to_string(tid) + ".sam");
            FilePath out_path_insertion = result_dir / ("thread" + std::to_string(tid) + "_insertion.sam");

            outs_path[static_cast<std::size_t>(tid)] = out_path;
            outs_with_insertion_path[static_cast<std::size_t>(tid)] = out_path_insertion;


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
            // - outs/outs_with_insertion 与 chunk 是跨线程共享只读/线程索引访问的；tid/i 为线程私有。
            #pragma omp parallel default(none) shared(outs, outs_with_insertion, chunk)
            {
                const int tid = omp_get_thread_num();
                auto& out = *outs[static_cast<std::size_t>(tid)];
                auto& out_insertion = *outs_with_insertion[static_cast<std::size_t>(tid)];

                // 这里每个 query 的处理耗时可能不均匀（不同长度/不同参考命中），用 dynamic(1) 先保证负载均衡。
                // 若后续确认为均匀负载，可调整为 guided/static 并基准测试。
                #pragma omp for schedule(dynamic, 1)
                for (std::int64_t i = 0; i < static_cast<std::int64_t>(chunk.size()); ++i)
                {
                    alignOneQueryToRef(chunk[static_cast<std::size_t>(i)], out, out_insertion);
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

    void RefAligner::mergeAlignedResults(const FilePath& aligned_consensus_path, const std::string& msa_cmd)
    {
        // 将多个线程的 SAM 文件合并为一个文件
        // 所有的序列都比对到多个参考序列上了，输入的参数为这些参考序列比对好的文件
        // 因此要解析这些参考的互相对应关系，然后把它们合并到一个文件中

        // 1. 第一步 星比对合并或者调用外部软件比对有插入的序列
        // 首先把insertion的sam文件合并为fasta
        FilePath result_dir = work_dir / RESULTS_DIR;

        bool using_other_align_insertion = true;

        FilePath aligned_insertion_fasta = result_dir / "aligned_insertion.fasta";

        if (using_other_align_insertion)
        {
            // 将所有 insertion SAM 文件路径收集到一个 vector
            std::vector<FilePath> insertion_sam_paths;
            for (const auto& path : outs_with_insertion_path)
            {
                insertion_sam_paths.push_back(path);
            }

            // 定义输出 FASTA 文件路径
            FilePath insertion_fasta_path = result_dir / "all_insertion.fasta";

            // 调用辅助函数：将共识序列和所有 SAM 文件合并为一个 FASTA 文件
            // 说明：
            // 1. 该函数会先写入共识序列（consensus_seq）
            // 2. 然后逐个读取 SAM 文件并提取序列追加写入
            // 3. 使用同一个 SeqWriter 实例，避免追加模式的复杂性
            // 4. 返回值为合并的总序列数（包括共识序列）
            const std::size_t total_sequences = mergeConsensusAndSamToFasta(
                insertion_sam_paths,
                insertion_fasta_path,
                80  // FASTA 行宽
            );

            #ifdef _DEBUG
            spdlog::info("mergeAlignedResults: merged {} sequences (1 consensus + {} from SAM files) to {}",
                        total_sequences, total_sequences - 1, insertion_fasta_path.string());
            #endif

            alignConsensusSequence(insertion_fasta_path, aligned_insertion_fasta, msa_cmd, threads);
        }
        else
        {

        }
    }

} // namespace align

