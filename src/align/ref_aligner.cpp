#include "align.h"

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

    void RefAligner::alignQueryToRef(const FilePath& qry_fasta_path)
    {
        // 读取每一条序列
        // 首先计算和所有参考序列的相似度，选择最相似的参考序列进行比对
        // 使用并行优化，有多少个线程就创建多少个文件，将比对结果写出到线程对应的文件
    }



} // namespace align


