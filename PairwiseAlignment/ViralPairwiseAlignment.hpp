#pragma once

#include <iostream>
#include <string>
#include <vector>

std::string call_wfa_to_get_cigar(std::vector<unsigned char> viral_sequence, std::vector<unsigned char> reference_sequence);
std::vector<std::string> cigar_to_vcf(const std::string& ref_id, const std::string& query_id, const std::vector<unsigned char>& ref_seq,
                                      const std::vector<unsigned char>& query_seq, const std::string& cigar);
std::vector<unsigned char> cigar_to_fasta(std::vector<unsigned char> viral_sequence, std::vector<unsigned char> reference_sequence,
                                          const std::string& cigar);
std::vector<unsigned char> apply_simple_cigar(std::vector<unsigned char> viral_sequence, const std::string& cigar);
// std::vector<unsigned char> fix_alignment_length(const std::vector<unsigned char>& viral_sequence,
//                                                const std::vector<unsigned char>& reference_sequence,
//                                                const std::string& cigar);

