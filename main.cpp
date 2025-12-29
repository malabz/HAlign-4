//
// Created by zhaiyixiao on 2025/7/2.
//

#include "PairwiseAlignment/ViralPairwiseAlignment.hpp"
#include "Utils/Arguments.hpp"
#include "Utils/CommandLine.hpp"
#include "Utils/Utils.hpp"
#include "multi-thread/ThreadPoolPerThread.hpp"

#include <filesystem>
#include <fstream>
#include <mutex>
#include <iostream>
#include <chrono>
#include <random>

#ifdef __linux__
#include <string>
void print_peak_memory_usage() {
    std::ifstream status_file("/proc/self/status");
    std::string line;
    while (std::getline(status_file, line)) {
        if (line.find("VmPeak:") == 0 || line.find("VmHWM:") == 0) {
            std::cout << "[MEM]  " << line << "\n";
        }
    }
}
#else
void print_peak_memory_usage() {
    std::cout << "[MEM]  Memory usage check is only supported on Linux.\n";
}
#endif

// 生成随机临时目录
std::string create_temp_dir() {
    std::string base = "tmp_";
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(10000, 99999);
    std::string dir = base + std::to_string(dist(gen));

    std::filesystem::create_directory(dir);
    return dir;
}

// 删除临时目录及内容
void remove_temp_dir(const std::string& dir) {
    try {
        std::filesystem::remove_all(dir);
    } catch (std::exception& e) {
        std::cerr << "[WARN] Failed to remove temp dir " << dir << ": " << e.what() << "\n";
    }
}

int main(int argc, char* argv[]) {
    std::cout << "[INFO] Program started.\n";
    auto start_point = std::chrono::high_resolution_clock::now();

    // ==== Step 1. 解析参数 ====
    SmpCommandLine userCommands(argc, argv);

    // ==== Step 2. 检查帮助或参数不足 ====
    if (userCommands.helpMessageWanted(4)) {   // 3个必需参数 + 程序名
        userCommands.showHelpMessage();
        return 0;
    }

    // ==== Step 3. 解析位置参数 ====
    arguments::in_file_name  = userCommands.getString(1, "", "Input file/folder path [.fasta or folder]");
    arguments::refer_file_name = userCommands.getString(2, "", "Reference file path [.fasta]");
    arguments::out_file_name = userCommands.getString(3, "", "Output file prefix");

    // ==== Step 4. 解析可选参数 ====
    int numThreads = userCommands.getInteger("t", "threads", 1, "Number of threads");
    bool save_vcf  = userCommands.getBoolean("s", "save-vcf", "Enable VCF output");

    // ==== Step 5. 参数检查 ====
    if (arguments::in_file_name.empty() || arguments::refer_file_name.empty() || arguments::out_file_name.empty()) {
        std::cerr << "[ERROR] Missing required input/output files.\n";
        userCommands.showHelpMessage();
        return 1;
    }

    // ==== Step 6. 输入输出路径检查 ====
    std::filesystem::path inputPath = std::filesystem::absolute(arguments::in_file_name);
    if (!std::filesystem::exists(inputPath)) {
        std::cerr << "[ERROR] Input file/folder does not exist.\n";
        return 1;
    }
    arguments::in_file_name = inputPath.generic_string();
    std::replace(arguments::in_file_name.begin(), arguments::in_file_name.end(), '\\', '/');
    if (std::filesystem::is_directory(inputPath) && arguments::in_file_name.back() != '/')
        arguments::in_file_name += '/';

    arguments::out_file_name = std::filesystem::absolute(arguments::out_file_name).generic_string();
    std::replace(arguments::out_file_name.begin(), arguments::out_file_name.end(), '\\', '/');

    std::cout << "[INFO] Input_file   : " << arguments::in_file_name << "\n";
    std::cout << "[INFO] Reference    : " << arguments::refer_file_name << "\n";
    std::cout << "[INFO] Output_file  : " << arguments::out_file_name << "\n";
    std::cout << "[INFO] Threads      : " << numThreads << "\n";
    std::cout << "[INFO] Output VCF   : " << (save_vcf ? "YES (-s enabled)" : "NO (default)") << "\n";

    // ==== Step 7. 临时文件夹 ====
    std::string temp_dir = create_temp_dir();
    std::cout << "[INFO] Temporary directory: " << temp_dir << "\n";

    // ==== Step 8. 读取参考序列 ====
    std::ifstream ref_ifs(arguments::refer_file_name);
    if (!ref_ifs) {
        std::cerr << "[ERROR] Cannot access reference file: " << arguments::refer_file_name << "\n";
        return 1;
    }
    std::string ref_id;
    std::vector<unsigned char> ref_seq;
    utils::read_single_fasta_sequence(ref_ifs, ref_id, ref_seq);

    // ==== Step5. 初始化线程池 & 临时输出文件 ====
    ThreadPoolPerThread pool(numThreads);
    std::vector<std::ofstream> thread_outputs(numThreads);
    std::vector<std::ofstream> thread_vcf;

    if (save_vcf)
        thread_vcf.resize(numThreads);

    for (int i = 0; i < numThreads; ++i) {
        std::string fasta_tmp = temp_dir + "/temp_out_" + std::to_string(i) + ".fasta.tmp";
        thread_outputs[i].open(fasta_tmp, std::ios::binary);

        if (save_vcf) {
            std::string vcf_tmp = temp_dir + "/thread_" + std::to_string(i) + ".vcf";
            thread_vcf[i].open(vcf_tmp);
        }
    }

    // ==== Step 9. 读取并分批处理输入序列 ====
    std::ifstream data_ifs(arguments::in_file_name);
    if (!data_ifs) {
        std::cerr << "[ERROR] Cannot open input file: " << arguments::in_file_name << "\n";
        return 1;
    }

    size_t batch_size = numThreads * 100;
    size_t total_sequences = 0;
    size_t batch_id = 0;

    while (true) {
        std::vector<std::pair<std::string, std::vector<unsigned char>>> batch;
        if (!utils::read_batch_of_n_sequences_with_ids(data_ifs, batch, batch_size))
            break;
        total_sequences += batch.size();

        for (size_t i = 0; i < batch.size(); ++i) {
            // std::cerr << "[DEBUG] before enqueue: "
            //   << "id=" << batch[i].first
            //   << " len=" << batch[i].second.size() << "\n";
            int tid = i % numThreads;
            pool.enqueue(tid, [seq_id = batch[i].first,
                               seq_data = batch[i].second,
                               tid, &ref_seq, &ref_id,
                               &thread_outputs, &thread_vcf, save_vcf]() {

                if (seq_id.empty() || seq_data.empty()) {
                    std::cerr << "[DEBUG] Empty seq_id or seq_data in thread " << tid << "\n";
                }
                
                std::string cigar = call_wfa_to_get_cigar(seq_data, ref_seq);
                if (cigar.empty()) {
                    std::cerr << "[DEBUG] Empty CIGAR for seq_id: " << seq_id << "\n";
                }
                    

                // 仅在需要时计算VCF
                if (save_vcf) {
                    auto vcf_results = cigar_to_vcf(ref_id, seq_id, ref_seq, seq_data, cigar);
                    auto& vcf_out = thread_vcf[tid];
                    for (auto& line : vcf_results)
                        vcf_out << line << "\n";
                }

                auto result = cigar_to_fasta(seq_data, ref_seq, cigar);
                if (result.empty()) {
                    std::cerr << "[DEBUG] cigar_to_fasta returned empty for " << seq_id << "\n";
                }
                auto& out = thread_outputs[tid];
                out << ">" << seq_id << "\n";
                for (auto c : result) out << c;
                out << "\n";
            });
        }
        ++batch_id;
        std::cout << "[INFO] Batch " << batch_id << " processed, total: " << total_sequences << "\n";
    }

    pool.wait_for_all();
    for (auto& ofs : thread_outputs) ofs.close();
    if (save_vcf) for (auto& of_vcf : thread_vcf) of_vcf.close();

    // ==== Step 10. 合并fasta输出 ====
    std::string fasta_output_name = arguments::out_file_name + ".fasta";
    std::ofstream fasta_out(fasta_output_name, std::ios::binary | std::ios::out);
    fasta_out << ">" << ref_id << "\n";
    for (auto c : ref_seq) fasta_out << c;
    fasta_out << "\n";

    for (int i = 0; i < numThreads; ++i) {
        std::ifstream tmp_in(temp_dir + "/temp_out_" + std::to_string(i) + ".fasta.tmp", std::ios::binary);
        fasta_out << tmp_in.rdbuf();
        fasta_out << "\n";
    }
    fasta_out.close();

    // ==== Step 11. （可选）合并VCF输出 ====
    if (save_vcf) {
        std::string vcf_output_name = arguments::out_file_name + ".vcf";
        std::ofstream vcf_out(vcf_output_name);
        vcf_out << "##fileformat=VCFv4.1\n";
        vcf_out << "##source=" << arguments::in_file_name << "\n";
        vcf_out << "##reference=" << arguments::refer_file_name << "\n";
        vcf_out << "##INFO=<ID=SEQID,Number=1,TYPE=SNP/INS/DEL,Type=String,Description=\"Query sequence ID\">\n";
        vcf_out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

        for (int tid = 0; tid < numThreads; ++tid) {
            std::string fname = temp_dir + "/thread_" + std::to_string(tid) + ".vcf";
            std::ifstream ifs(fname);
            if (!ifs.is_open()) continue;

            std::string line;
            while (std::getline(ifs, line)) {
                if (line.empty()) continue;
                vcf_out << line << "\n";
            }
            ifs.close();
        }
        vcf_out.close();
        std::cout << "[INFO] VCF output merged successfully.\n";
    }

    // ==== Step 12. 清理临时目录 ====
    remove_temp_dir(temp_dir);

    std::cout << "[INFO] Completed. Total sequences processed: " << total_sequences << "\n";
    auto end_point = std::chrono::high_resolution_clock::now();
    std::cout << "[TIME] Total runtime: "
              << std::chrono::duration_cast<std::chrono::seconds>(end_point - start_point).count()
              << " seconds.\n";
    return 0;
}
