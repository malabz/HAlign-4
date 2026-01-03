#ifndef HALIGN4_PREPROCESS_H
#define HALIGN4_PREPROCESS_H

#include <cstddef>
#include "config.hpp"
#include "utils.h"
#include "consensus.h"

// 预处理输入 FASTA，并返回处理的序列数量（total records processed）
uint_t preprocessInputFasta(const std::string input_path, const std::string workdir, const int cons_n = 1000);

// 对外：对输入 fasta 文件进行 MSA（consensus 对齐），将结果写入 output 文件
void alignConsensusSequence(const FilePath& input_file, const FilePath& output_file,
                            const std::string& msa_cmd, const std::string& workdir, int threads);


#endif //HALIGN4_PREPROCESS_H
