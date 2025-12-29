//
// Created by zhaiyixiao on 2025/6/24.
//

#include "ViralPairwiseAlignment.hpp"
#include "WFA2-lib/bindings/cpp/WFAligner.hpp"
#include <vector>
#include <iostream>
#include <memory>
#include <stdexcept>

std::string call_wfa_to_get_cigar(std::vector<unsigned char> viral_sequence, std::vector<unsigned char> reference_sequence)
{

    std::string cigar = "";

    // Convert sequences to strings
    std::string refer(reference_sequence.begin(), reference_sequence.end());
    std::string curr_viral(viral_sequence.begin(), viral_sequence.end());

    // std::cout << "Sequence length: Reference=" << refer_length << ", Viral=" << curr_length << std::endl;

    // Create WFA aligner with error handling
    try {
        // Try MemoryLow first 
        thread_local wfa::WFAlignerGapAffine aligner_low(4,6,2,wfa::WFAligner::Alignment,wfa::WFAligner::MemoryLow);
        aligner_low.alignEnd2End(refer, curr_viral);
        cigar = aligner_low.getAlignment();

        // if (!cigar.empty()) {
        //     std::cout << "Successfully using MemoryLow mode, CIGAR length: " << cigar.size() << std::endl;
        // }
        
    } catch (const std::exception& e) {
        std::cout << "MemoryHigh mode failed: " << e.what() << std::endl;
        try {
            // Fallback to MemoryHigh
            thread_local wfa::WFAlignerGapAffine aligner_high(4,6,2,wfa::WFAligner::Alignment,wfa::WFAligner::MemoryHigh);
            aligner_high.alignEnd2End(refer, curr_viral);
            cigar = aligner_high.getAlignment();

            if (!cigar.empty()) {
                std::cout << "Successfully using MemoryHigh mode, CIGAR length: " << cigar.size() << std::endl;
            }
        } catch (const std::exception& e) {
            std::cout << "MemoryLow mode failed: " << e.what() << std::endl;
            throw;
        }
    }

    return cigar;
}

std::vector<std::string> cigar_to_vcf(const std::string& ref_id,
                                      const std::string& query_id,
                                      const std::vector<unsigned char>& ref_seq,
                                      const std::vector<unsigned char>& query_seq,
                                      const std::string& cigar) {
    std::vector<std::string> vcf_lines;
    int ref_pos = 0;   // 参考序列当前位置（0-based）
    int query_pos = 0; // 查询序列当前位置（0-based）

    for (size_t i = 0; i < cigar.size(); ++i) {
        char op = cigar[i];

        if (op == 'M' || op == 'X' || op == '=') {
            // 匹配或错配
            if (ref_pos < (int)ref_seq.size() && query_pos < (int)query_seq.size()) {
                char rbase = static_cast<char>(ref_seq[ref_pos]);
                char qbase = static_cast<char>(query_seq[query_pos]);
                if (rbase != qbase) {
                    // 这是一个 SNP
                    std::string line = ref_id + "\t" +
                                       std::to_string(ref_pos + 1) + "\t.\t" + // VCF 1-based 坐标
                                       std::string(1, rbase) + "\t" +
                                       std::string(1, qbase) + "\t.\tPASS\tSEQID=" + query_id + ", TYPE=SNP";
                    vcf_lines.push_back(line);
                }
            }
            ref_pos++;
            query_pos++;
        } else if (op == 'I') {
            // 插入
            size_t ins_start = i;
            while (i < cigar.size() && cigar[i] == 'I') i++;
            size_t ins_len = i - ins_start;
            i--;

            if (ref_pos > 0 && query_pos + ins_len <= query_seq.size()) {
                std::string ref(1, static_cast<char>(ref_seq[ref_pos - 1]));
                std::string ins(query_seq.begin() + query_pos,
                                query_seq.begin() + query_pos + ins_len);
                std::string alt = ref + ins;

                std::string line = ref_id + "\t" +
                                   std::to_string(ref_pos) + "\t.\t" +
                                   ref + "\t" + alt + "\t.\tPASS\tSEQID=" + query_id + ", TYPE=INS";
                vcf_lines.push_back(line);
            }
            query_pos += ins_len;
        } else if (op == 'D') {
            // 缺失
            size_t del_start = i;
            while (i < cigar.size() && cigar[i] == 'D') i++;
            size_t del_len = i - del_start;
            i--;

            if (ref_pos > 0 && ref_pos + del_len <= (int)ref_seq.size()) {
                std::string ref;
                for (size_t k = 0; k < del_len + 1; k++) {
                    ref.push_back(static_cast<char>(ref_seq[ref_pos - 1 + k]));
                }
                std::string alt(1, static_cast<char>(ref_seq[ref_pos - 1]));

                std::string line = ref_id + "\t" +
                                   std::to_string(ref_pos) + "\t.\t" +
                                   ref + "\t" + alt + "\t.\tPASS\tSEQID=" + query_id + ", TYPE=DEL";
                vcf_lines.push_back(line);
            }
            ref_pos += del_len;
        }
    }

    return vcf_lines;
}


std::vector<unsigned char> cigar_to_fasta(std::vector<unsigned char> viral_sequence, std::vector<unsigned char> reference_sequence,
                                          const std::string& cigar)
{
    try {
        std::vector<unsigned char> curr_aligned_sequence;

        // Convert sequences to strings
        std::string refer(reference_sequence.begin(), reference_sequence.end());
        std::string curr_viral(viral_sequence.begin(), viral_sequence.end());

        const int refer_length = refer.size();
        const int curr_length = curr_viral.size();

        // std::cout << "Sequence length: Reference=" << refer_length << ", Viral=" << curr_length << std::endl;

        if (cigar.empty()) {
            std::cerr << "Warning: Failed to get alignment, returning padded sequence" << std::endl;
            // Return a sequence padded to reference length
            curr_aligned_sequence = viral_sequence;
            while (curr_aligned_sequence.size() < reference_sequence.size()) {
                curr_aligned_sequence.push_back('-');
            }
            return curr_aligned_sequence;
        }

        // Apply CIGAR and get aligned sequence
        curr_aligned_sequence = apply_simple_cigar(viral_sequence, cigar);

        // Ensure output length matches reference length
        if (curr_aligned_sequence.size() != reference_sequence.size()) {
            std::cout << "Warning: Length mismatch. CIGAR length=" << cigar.size()
                      << ", Output length=" << curr_aligned_sequence.size()
                      << ", Expected length=" << reference_sequence.size() << std::endl;

            // Fix length
            if (curr_aligned_sequence.size() < reference_sequence.size()) {
                while (curr_aligned_sequence.size() < reference_sequence.size()) {
                    curr_aligned_sequence.push_back('-');
                }
            } else {
                curr_aligned_sequence.resize(reference_sequence.size());
            }
        }

        return curr_aligned_sequence;

    } catch (const std::exception& e) {
        std::cerr << "Error in alignment: " << e.what() << std::endl;
        // Return a sequence padded to reference length
        std::vector<unsigned char> error_sequence = viral_sequence;
        while (error_sequence.size() < reference_sequence.size()) {
            error_sequence.push_back('-');
        }
        return error_sequence;
    }
}


std::vector<unsigned char> apply_simple_cigar(std::vector<unsigned char> viral_sequence, const std::string& cigar)
{
    try {
        std::vector<unsigned char> aligned_sequence;
        aligned_sequence.reserve(viral_sequence.size() * 2); // Pre-allocate memory

        size_t seq_pos = 0;

        for (char op : cigar) {
            switch (op) {
                case 'M': // Match or mismatch
                case 'X': // Mismatch
                case '=': // Match
                    if (seq_pos < viral_sequence.size()) {
                        aligned_sequence.push_back(viral_sequence[seq_pos]);
                        seq_pos++;
                    } else {
                        aligned_sequence.push_back('-');
                    }
                    break;

                case 'I': // Insertion (skip in query)
                    if (seq_pos < viral_sequence.size()) {
                        seq_pos++;
                    }
                    break;

                case 'D': // Deletion (gap in query)
                    aligned_sequence.push_back('-');
                    break;

                default:
                    std::cerr << "Warning: Unknown CIGAR operation: " << op << std::endl;
                    break;
            }
        }

        return aligned_sequence;

    } catch (const std::exception& e) {
        std::cerr << "Error in CIGAR processing: " << e.what() << std::endl;
        return viral_sequence; // Return original sequence on error
    }
}

