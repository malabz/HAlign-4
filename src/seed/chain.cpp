// ==================================================================
// chain.cpp - 基于 minimap2 的锚点收集与链化实现
// ==================================================================
//
// 功能说明：
// 本文件实现了序列比对中的锚点（anchor）收集功能，参考 minimap2 的设计。
// 锚点是 ref 和 query 序列之间共享的种子匹配位置，用于后续的链化（chaining）
// 和精确比对。
//
// 参考来源：
// - minimap2/seed.c: mm_collect_matches() - 收集种子匹配
// - minimap2/hit.c: 锚点处理和过滤
// - minimap2/lchain.c: 链化算法（DP 动态规划）
//
// 核心数据流程（minimap2 风格）：
// 1. 将 ref_hits 按 hash 排序，构建 hash -> position 的索引
// 2. 遍历 qry_hits，通过 hash 查找匹配的 ref_hits
// 3. 对于每个匹配，生成一个 Anchor
// 4. 返回所有锚点（未排序，调用者可根据需要排序）
//
// 设计说明：
// - 使用 C++ 标准库（std::sort, std::unordered_map）替代 minimap2 的 ksort/khash
// - 保持与 minimap2 相同的语义：hash 匹配表示潜在的同源区域
// - 模板化设计：支持任意继承 SeedHitBase 的种子类型
//
// 性能优化：
// - 预先对 ref_hits 按 hash 排序，使用二分查找加速匹配
// - 避免重复的 hash 计算（使用缓存的 hash 值）
// - 预分配输出容器，减少内存重新分配
// ==================================================================

#include "seed.h"
#include <algorithm>
#include <unordered_map>
#include <cstdint>

namespace seed {

// ==================================================================
// 内部辅助结构：用于 ref_hits 的 hash 索引
// ==================================================================
// 说明：
// 为了高效查找 ref_hits 中与 qry_hit 具有相同 hash 的元素，
// 我们构建一个 hash -> (start_index, count) 的索引表。
// 这类似于 minimap2 的 mm_idx_get() 函数的功能。
//
// 实现：
// 1. 先对 ref_hits 按 hash 排序（稳定排序保持原有位置顺序）
// 2. 遍历一次，记录每个 hash 的起始位置和数量
// 3. 查询时通过 unordered_map O(1) 查找
// ==================================================================
struct HashIndex {
    std::size_t start;  // 在排序后数组中的起始索引
    std::size_t count;  // 具有该 hash 的元素数量
};

// ==================================================================
// 函数：collect_anchors（模板实现）
// ==================================================================
// 功能：
// 从 ref_hits 和 qry_hits 中收集所有共享相同 hash 的锚点对
//
// 算法流程（参考 minimap2 的 mm_collect_matches）：
// 1. 对 ref_hits 按 hash 排序
// 2. 构建 hash -> (start, count) 的索引
// 3. 遍历 qry_hits，对每个 qry_hit：
//    a. 在索引中查找具有相同 hash 的 ref_hits
//    b. 对于每个匹配的 ref_hit，生成一个 Anchor
// 4. 返回所有锚点
//
// 参数：
//   @param ref_hits - 参考序列的种子命中列表
//   @param qry_hits - 查询序列的种子命中列表
//
// 返回：
//   std::vector<Anchor> - 锚点列表
//
// 复杂度：
// - 时间：O(R log R + Q * avg_matches)
//   - R = ref_hits.size(), Q = qry_hits.size()
//   - avg_matches = 每个 query hit 平均匹配的 ref hits 数量
// - 空间：O(R + A)，R 为 ref_hits 大小，A 为锚点数量
//
// 注意事项：
// 1. 输入的 ref_hits 和 qry_hits 不会被修改（使用拷贝排序）
// 2. 输出的锚点按照遍历 qry_hits 的顺序生成
// 3. 如果需要按 (pos_ref, pos_qry) 或对角线排序，调用者需自行排序
// 4. 对于高频 hash（高重复区域），可能生成大量锚点；建议后续过滤
// ==================================================================

// ------------------------------------------------------------------
// 显式实例化声明（模板定义在头文件中）
// 这里提供 MinimizerHit 的具体实现
// ------------------------------------------------------------------

template <>
std::vector<Anchor> collect_anchors<minimizer::MinimizerHit>(
    const std::vector<minimizer::MinimizerHit>& ref_hits,
    const std::vector<minimizer::MinimizerHit>& qry_hits)
{
    std::vector<Anchor> anchors;

    // ------------------------------------------------------------------
    // 边界条件：如果任一输入为空，直接返回空结果
    // ------------------------------------------------------------------
    if (ref_hits.empty() || qry_hits.empty()) {
        return anchors;
    }

    // ------------------------------------------------------------------
    // 步骤1：对 ref_hits 按 hash 排序（参考 minimap2 的 radix_sort_128x）
    // ------------------------------------------------------------------
    // 说明：
    // minimap2 使用基数排序对 mm128_t 按 x（包含 hash）排序
    // 我们使用 std::sort + 自定义比较器实现相同效果
    //
    // 性能说明：
    // - std::sort 通常是 O(N log N) 的 introsort
    // - 对于大规模数据，基数排序可能更快，但这里优先选择标准库
    // ------------------------------------------------------------------
    std::vector<minimizer::MinimizerHit> sorted_ref = ref_hits;  // 拷贝以避免修改输入
    std::sort(sorted_ref.begin(), sorted_ref.end(),
        [](const minimizer::MinimizerHit& a, const minimizer::MinimizerHit& b) {
            return a.hash() < b.hash();
        });

    // ------------------------------------------------------------------
    // 步骤2：构建 hash 索引表（参考 minimap2 的 mm_idx_get）
    // ------------------------------------------------------------------
    // 说明：
    // minimap2 使用 khash 构建 hash 到 (offset, count) 的映射
    // 我们使用 std::unordered_map 实现相同功能
    //
    // 数据结构：
    // - key: hash 值（64-bit）
    // - value: HashIndex{start, count}
    //   - start: 在 sorted_ref 中的起始索引
    //   - count: 具有该 hash 的连续元素数量
    // ------------------------------------------------------------------
    std::unordered_map<hash_t, HashIndex> hash_index;
    hash_index.reserve(sorted_ref.size());  // 预分配，减少 rehash

    if (!sorted_ref.empty()) {
        std::size_t start = 0;
        hash_t current_hash = sorted_ref[0].hash();

        for (std::size_t i = 1; i <= sorted_ref.size(); ++i) {
            // 当遇到新的 hash 或到达末尾时，记录上一个 hash 的索引
            if (i == sorted_ref.size() || sorted_ref[i].hash() != current_hash) {
                hash_index[current_hash] = HashIndex{start, i - start};
                if (i < sorted_ref.size()) {
                    start = i;
                    current_hash = sorted_ref[i].hash();
                }
            }
        }
    }

    // ------------------------------------------------------------------
    // 步骤3：遍历 qry_hits，收集锚点（参考 minimap2 的 mm_collect_matches）
    // ------------------------------------------------------------------
    // 说明：
    // 对于每个 qry_hit，查找具有相同 hash 的所有 ref_hits
    // 对于每个匹配对 (ref_hit, qry_hit)，生成一个 Anchor
    //
    // Anchor 字段说明（参考 minimap2 的 mm128_t 布局）：
    // - hash: 共享的 hash 值（用于验证和调试）
    // - rid_ref: ref 序列 ID（来自 ref_hit.rid()）
    // - pos_ref: ref 上的位置（来自 ref_hit.pos()）
    // - rid_qry: query 序列 ID（来自 qry_hit.rid()）
    // - pos_qry: query 上的位置（来自 qry_hit.pos()）
    // - span: 覆盖长度，取 min(ref_hit.span(), qry_hit.span())
    // - is_rev: 是否反向互补，通过 ref_hit.strand() XOR qry_hit.strand() 判断
    //
    // 性能优化：
    // - 预估锚点数量，预分配 anchors 容器
    // - 使用 emplace_back 避免拷贝
    // ------------------------------------------------------------------

    // 预估锚点数量：假设每个 qry_hit 平均匹配 1 个 ref_hit
    // 实际数量可能更多（重复区域）或更少（无匹配）
    anchors.reserve(qry_hits.size());

    for (const auto& qry_hit : qry_hits) {
        const hash_t qry_hash = qry_hit.hash();

        // 在索引中查找具有相同 hash 的 ref_hits
        auto it = hash_index.find(qry_hash);
        if (it == hash_index.end()) {
            // 没有匹配的 ref_hit，跳过
            continue;
        }

        const HashIndex& idx = it->second;

        // 遍历所有匹配的 ref_hits，生成锚点
        for (std::size_t i = idx.start; i < idx.start + idx.count; ++i) {
            const auto& ref_hit = sorted_ref[i];

            // ------------------------------------------------------------------
            // 生成 Anchor（参考 minimap2 的锚点布局）
            // ------------------------------------------------------------------
            // 说明：
            // minimap2 的 mm128_t 用于存储锚点：
            //   a[].x: rev<<63 | tid<<32 | tpos
            //   a[].y: flags<<40 | q_span<<32 | q_pos
            //
            // 我们使用更清晰的结构体表示：
            //   Anchor.hash:     共享的 hash 值
            //   Anchor.rid_ref:  ref 序列 ID
            //   Anchor.pos_ref:  ref 位置
            //   Anchor.rid_qry:  query 序列 ID
            //   Anchor.pos_qry:  query 位置
            //   Anchor.span:     覆盖长度
            //   Anchor.is_rev:   是否反向
            // ------------------------------------------------------------------
            Anchor anchor;
            anchor.hash = qry_hash;
            anchor.rid_ref = ref_hit.rid();
            anchor.pos_ref = ref_hit.pos();
            anchor.rid_qry = qry_hit.rid();
            anchor.pos_qry = qry_hit.pos();

            // span 取两者较小值（保守估计覆盖范围）
            anchor.span = std::min(ref_hit.span(), qry_hit.span());

            // 方向判断：如果 ref 和 query 的 strand 不同，则为反向互补匹配
            // 这对应 minimap2 中的 rev 标志
            anchor.is_rev = (ref_hit.strand() != qry_hit.strand());

            anchors.emplace_back(anchor);
        }
    }

    return anchors;
}

// ==================================================================
// 辅助函数：按对角线排序锚点（参考 minimap2 的 radix_sort_128x）
// ==================================================================
// 功能：
// 将锚点按对角线 (pos_ref - pos_qry) 和 pos_ref 排序
// 这是链化算法的常见预处理步骤
//
// 排序规则（minimap2 风格）：
// 1. 首先按 rid_ref 排序（不同参考序列分开处理）
// 2. 其次按对角线 (pos_ref - pos_qry) 排序
// 3. 最后按 pos_ref 排序
//
// 参数：
//   @param anchors - 输入/输出锚点列表（原地排序）
// ==================================================================
void sortAnchorsByDiagonal(std::vector<Anchor>& anchors)
{
    std::sort(anchors.begin(), anchors.end(),
        [](const Anchor& a, const Anchor& b) {
            // 1. 首先按 rid_ref 排序
            if (a.rid_ref != b.rid_ref) {
                return a.rid_ref < b.rid_ref;
            }
            // 2. 其次按对角线排序
            // 对角线 = pos_ref - pos_qry（注意：可能为负数，需要有符号比较）
            const int64_t diag_a = static_cast<int64_t>(a.pos_ref) - static_cast<int64_t>(a.pos_qry);
            const int64_t diag_b = static_cast<int64_t>(b.pos_ref) - static_cast<int64_t>(b.pos_qry);
            if (diag_a != diag_b) {
                return diag_a < diag_b;
            }
            // 3. 最后按 pos_ref 排序
            return a.pos_ref < b.pos_ref;
        });
}

// ==================================================================
// 辅助函数：按位置排序锚点（参考 minimap2 的链化前排序）
// ==================================================================
// 功能：
// 将锚点按 (rid_ref, pos_ref, pos_qry) 排序
// 这是另一种常见的排序方式，用于某些链化算法
//
// 参数：
//   @param anchors - 输入/输出锚点列表（原地排序）
// ==================================================================
void sortAnchorsByPosition(std::vector<Anchor>& anchors)
{
    std::sort(anchors.begin(), anchors.end(),
        [](const Anchor& a, const Anchor& b) {
            // 1. 首先按 rid_ref 排序
            if (a.rid_ref != b.rid_ref) {
                return a.rid_ref < b.rid_ref;
            }
            // 2. 其次按 pos_ref 排序
            if (a.pos_ref != b.pos_ref) {
                return a.pos_ref < b.pos_ref;
            }
            // 3. 最后按 pos_qry 排序
            return a.pos_qry < b.pos_qry;
        });
}

// ==================================================================
// 辅助函数：过滤高频锚点（参考 minimap2 的 mm_seed_select）
// ==================================================================
// 功能：
// 移除出现次数超过阈值的锚点（通常是重复区域导致的噪声）
// minimap2 使用类似策略来处理高频 minimizer
//
// 参数：
//   @param anchors - 输入/输出锚点列表（原地修改）
//   @param max_occ - 单个 hash 允许的最大出现次数（默认 500）
//
// 说明：
// 如果某个 hash 在 anchors 中出现超过 max_occ 次，
// 则移除所有具有该 hash 的锚点（认为是重复区域，不可靠）
// ==================================================================
void filterHighFrequencyAnchors(std::vector<Anchor>& anchors, std::size_t max_occ)
{
    if (anchors.empty() || max_occ == 0) {
        return;
    }

    // 1. 统计每个 hash 的出现次数
    std::unordered_map<hash_t, std::size_t> hash_count;
    for (const auto& anchor : anchors) {
        hash_count[anchor.hash]++;
    }

    // 2. 移除高频锚点
    auto new_end = std::remove_if(anchors.begin(), anchors.end(),
        [&hash_count, max_occ](const Anchor& anchor) {
            return hash_count[anchor.hash] > max_occ;
        });

    anchors.erase(new_end, anchors.end());
}

} // namespace seed

