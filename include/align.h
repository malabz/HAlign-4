#ifndef HALIGN4_ALIGN_H
#define HALIGN4_ALIGN_H
#include "utils.h"
#include "mash.h"
#include "seed.h"
#include <filesystem>
#include <string>
#include <vector>
#include <functional>

namespace align {
    using SeedHit = minimizer::MinimizerHit;
    using SeedHits = std::vector<SeedHit>;
    static constexpr seed::SeedKind kSeedKind = seed::SeedKind::minimizer;

    class RefAligner
    {
        public:
        // 初始化函数
        RefAligner(const FilePath& work_dir, const FilePath& ref_fasta_path, int kmer_size = 21, int window_size = 10,
                    int sketch_size = 2000, bool noncanonical = true);

        void alignQueryToRef(const FilePath& qry_fasta_path);


        private:
        FilePath work_dir;
        seq_io::SeqRecords ref_sequences;
        mash::Sketches ref_sketch;
        std::vector<SeedHits> ref_minimizers;

        int kmer_size = 21;
        int window_size = 10;
        int sketch_size = 2000;
        int random_seed = 42;

        bool keep_first_length = false;
        bool keep_all_length = false;
        bool noncanonical = true;


    };

} // namespace align

#endif //HALIGN4_ALIGN_H

