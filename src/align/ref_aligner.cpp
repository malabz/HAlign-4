#include "align.h"

namespace align {

    RefAligner::RefAligner(const FilePath& work_dir, const FilePath& ref_fasta_path, int kmer_size, int window_size, int sketch_size, bool noncanonical, bool keep_first_length, bool keep_all_length)
        : work_dir(work_dir), kmer_size(kmer_size), window_size(window_size), sketch_size(sketch_size), noncanonical(noncanonical), keep_first_length(keep_first_length), keep_all_length(keep_all_length)
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


} // namespace align


