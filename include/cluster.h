#ifndef HALIGN4_CLUSTER_H
#define HALIGN4_CLUSTER_H

#include <filesystem>
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <mutex>
#include <functional>

namespace cluster {

// FilePath 与项目中其它地方保持一致（使用 std::filesystem::path）
using FilePath = std::filesystem::path;

// Cluster 类用于：
// - 将一个“共识序列”（consensus）与一个输入的 FASTA 文件中每条序列进行比较，
//   比较基于 minimizer（或其他 k-mer 索引）来快速估计相似度；
// - 根据相似度阈值将输入序列聚类到不同组；
// - 将每个聚类写出到一个 FASTA 文件（或一个合并的输出文件，每个序列在 header 中带上聚类 id）。
//
// 新增：支持按百分比阈值分箱（例如 99%, 98%, 每1%）——实现思路是计算序列与共识的相似度 s ∈ [0,1]，
// 将其映射为整数百分比 p = floor(s*100)，然后把序列分配到用户指定的阈值列表中：
// 例如阈值列表 {99,98,95} 表示：
//  - 若 p >= 99 -> 归入 99 桶
//  - 否则若 p >= 98 -> 归入 98 桶
//  - 否则若 p >= 95 -> 归入 95 桶
//  - 否则 -> 可归入一个默认的 "<95" 桶（实现可选择）
// 如果不提供阈值列表（默认 empty），类将自动按每 1% 建立阈值（100,99,...,0），相当于每1%为一类。

class Cluster {
public:
    // 可注入的 minimizer 提取函数类型：输入序列字符串，返回一组 minimizer（用 64 位哈希或其他表示）
    using MinimizerExtractor = std::function<std::vector<std::uint64_t>(const std::string&)>;

    Cluster(const FilePath& workdir,
            std::size_t k = 15,
            std::size_t w = 10);

    ~Cluster();

    std::size_t runCluster(const std::string& consensus,
                    const FilePath& fasta_path,
                    int thread = 0,
                    const std::vector<int>& percent_bins = {});

    std::vector<std::uint64_t> extractMinimizers(const std::string& seq) const;

    // 估计两个 minimizer 集合的相似度（0..1）。
    double estimateSimilarity(const std::vector<std::uint64_t>& a,
                              const std::vector<std::uint64_t>& b) const;

    // 为单元测试或外部调用暴露：手动设置 minimizer 提取器
    void setMinimizerExtractor(MinimizerExtractor extractor);

private:
    // 成员变量
    FilePath workdir;               // 工作目录：构造时传入
    std::size_t k = 15;             // k-mer 大小
    std::size_t w = 10;             // 窗口大小
    int thread = 0;                // 构造时的线程数设置（0 表示自动选择）

    std::string consensus;          // 当前 run 使用的共识序列（在 run() 中设置）
    MinimizerExtractor extractor;   // minimizer 提取器（用户可注入）
};

} // namespace cluster

#endif //HALIGN4_CLUSTER_H

