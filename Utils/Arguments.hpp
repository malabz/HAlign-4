#pragma once
#include <string>

namespace arguments
{
    extern std::string in_file_name;
    extern std::string out_file_name;
    // add for viral msa --- star---
    extern std::string refer_file_name;
    // add for viral msa --- end---
    extern std::string tmp_file_name;
    extern std::string score_file;
    extern std::string snp_file_name;
    extern bool output_matrix;
    extern size_t ALL_LEN;
}