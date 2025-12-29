#include "Arguments.hpp"
//std::string arguments::in_file_name = "/home/zqzhoutong/stmsa-BWT/data/people13.fasta";
//std::string arguments::out_file_name = "/home/zqzhoutong/stmsa-BWT/data/test/people13ed.fasta";
std::string arguments::in_file_name;
// add for viral msa --- star---
std::string arguments::refer_file_name;
// add for viral msa --- end---
std::string arguments::out_file_name;
std::string arguments::tmp_file_name;
std::string arguments::score_file;
std::string arguments::snp_file_name;
size_t arguments::ALL_LEN = 0;
bool arguments::output_matrix;
