// ==================================================================
// chain.cpp - minimap2 风格的锚点收集（seeding）与预处理
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
// 核心数据流程（minimap2 风格，参考 README Algorithm overview）：
// 1) 对参考序列 minimizers 建索引（hash -> occurrences）
// 2) 对每个 query minimizer 查索引；若该参考 minimizer 不在 top -f 最频繁集合内，
//    则收集其在参考序列的所有出现位置作为 seeds/anchors
// 3) seeds/anchors 随后按参考坐标排序并做 DP chaining
//
// 重要：
// - minimap2 的频次过滤（-f/-U/--q-occ-frac/-e）发生在“展开 occurrences 之前”，
//   否则重复区域会导致 anchors 数量爆炸。
//
// 设计说明：
// - 使用 C++ 标准库（std::sort, std::unordered_map）替代 minimap2 的 ksort/khash
// - 保持与 minimap2 相同的语义：hash 匹配表示潜在的同源区域
// - 模板化设计：支持任意继承 SeedHitBase 的种子类型
//
// 性能优化：
// - 预先对 ref_hits 按 hash 排序
// - 避免重复的 hash 计算（使用缓存的 hash 值）
// - 预分配输出容器，减少内存重新分配
// ==================================================================

#include "seed.h"
#include <algorithm>
#include <unordered_map>
#include <cstdint>
#include <vector>
#include <cmath>
#include <limits>

namespace anchor {

// ==================================================================

// ------------------------------------------------------------------
// 函数：compute_occ_cutoff_top_frac
// 功能：根据 distinct minimizer 的 occurrence 列表，计算“最高频 top f_top_frac”对应的频次阈值。
//
// 输入：
// - occs: 每个 distinct minimizer 在参考中出现的次数（例如：{1,1,2,10,3,...}）
// - f_top_frac: top fraction（例如 2e-4 表示忽略最频繁的 0.02% minimizers）
//
// 输出：
// - 返回 occ_cutoff：当 a minimizer 的 occurrence >= occ_cutoff 时，可认为其属于 top 高频区域
//   （上层可据此决定是否过滤/稀疏采样）。
//
// 关键语义（对齐 minimap2 的直觉）：
// - f_top_frac==0：不做基于 top fraction 的过滤 => 返回 +inf
// - occs 很少或算出来 n_skip==0：也相当于“不过滤” => 返回 +inf
// - f_top_frac>=1：极端情况，几乎全部都要过滤 => 返回 1
//
// 实现细节：
// - 我们需要“第 n_skip 大的 occurrence 值”，用 nth_element 做 O(N) 期望时间选择。
// - 使用 std::greater<std::size_t>() 让 nth_element 得到降序意义下的第 n_skip 大。
//
// 复杂度：
// - 时间：O(N) 期望（nth_element），最坏 O(N log N)
// - 空间：O(N)（拷贝一份 tmp，避免修改输入）
// ------------------------------------------------------------------
std::size_t compute_occ_cutoff_top_frac(const std::vector<std::size_t>& occs,
                                        double f_top_frac)
{
    if (occs.empty()) return std::numeric_limits<std::size_t>::max();
    if (f_top_frac <= 0.0) return std::numeric_limits<std::size_t>::max();
    if (f_top_frac >= 1.0) return 1; // 极端情况：几乎全丢

    // top f fraction of DISTINCT minimizers
    const std::size_t n = occs.size();

    // n_skip 表示需要“忽略掉的 distinct minimizer 数量”
    // 例：n=10000, f=2e-4 => n_skip=f*n=2 => 忽略出现次数排名前 2 的 minimizer
    const std::size_t n_skip = static_cast<std::size_t>(std::floor(f_top_frac * static_cast<double>(n)));
    if (n_skip == 0) return std::numeric_limits<std::size_t>::max();

    // nth_element 会重排 tmp，使得 tmp[n_skip-1] 是“第 n_skip 大”的元素
    std::vector<std::size_t> tmp = occs;
    std::nth_element(tmp.begin(), tmp.begin() + static_cast<std::ptrdiff_t>(n_skip - 1), tmp.end(),
                     std::greater<std::size_t>());
    return tmp[n_skip - 1];
}

// ------------------------------------------------------------------
// 函数：compute_ref_occ_threshold
// 功能：计算最终“参考端 occurrence 阈值”。
//
// minimap2 里常见逻辑可理解为：
// - 先用 -f 估计一个阈值 f_cutoff（忽略 top fraction 的高频 minimizer）
// - 再用 -U 设置上下界，把阈值限制到 [u_floor, u_ceil] 范围内
// - 最终阈值 = max(u_floor, min(u_ceil, f_cutoff))
//
// 这样做的好处：
// - 当参考序列很大且重复很多，-f 能自适应地忽略最频繁的 minimizer
// - 同时 -U 能确保阈值不会过小（误杀正常 seeds）也不会过大（让 repeats 爆炸）
//
// 复杂度：
// - 主要来自 compute_occ_cutoff_top_frac 的 nth_element
// ------------------------------------------------------------------
std::size_t compute_ref_occ_threshold(const std::vector<std::size_t>& occs,
                                      const SeedFilterParams& p)
{
    const std::size_t f_cutoff = compute_occ_cutoff_top_frac(occs, p.f_top_frac);
    const std::size_t capped = std::min(p.u_ceil, f_cutoff);
    return std::max(p.u_floor, capped);
}

// ==================================================================
// 辅助函数：按对角线排序锚点（minimap2 chaining 预处理常用）
// ==================================================================
//
// 为什么要按对角线排序？
// - chaining 本质是把“在相同对角线附近”的 anchors 串起来形成一条链。
// - 对角线（diagonal）可理解为 ref_pos - qry_pos（正向）
//   对于同源区域，锚点通常会聚集在相近的 diagonal 上。
//
// 注意反向链（is_rev==true）的对角线定义不同：
// - minimap2 内部会把反向链 query 坐标映射到反向互补坐标系。
// - 在不知道 query 总长度 qlen 的情况下，我们无法直接算 (ref - q_rc)。
// - 但排序只需要“单调等价”的 key：
//   q_rc = qlen - (q + span)  => ref - q_rc = ref + q + span - qlen
//   qlen 是常数，可以忽略，因此用 (ref + q + span) 作为反向链的 diagonal key。
//
// 排序键：
// 1) rid_ref（不同参考序列分开）
// 2) is_rev（正向/反向分开，避免混链）
// 3) diagonal key（见上）
// 4) pos_ref / pos_qry（稳定化）
//
// 复杂度：O(A log A)，A=anchors.size()
// ------------------------------------------------------------------
void sortAnchorsByDiagonal(Anchors& anchors)
{
    std::sort(anchors.begin(), anchors.end(),
        [](const Anchor& a, const Anchor& b) {
            if (a.rid_ref != b.rid_ref) return a.rid_ref < b.rid_ref;
            if (a.is_rev != b.is_rev)   return a.is_rev < b.is_rev;

            // 正向：diag = r - q
            // 反向：等价于 r - q_rc，其中 q_rc = qlen - (q + span)
            //       qlen 为常数，排序可用 (r + q + span) 代替
            const int64_t diag_a = a.is_rev
                ? (static_cast<int64_t>(a.pos_ref) + static_cast<int64_t>(a.pos_qry) + static_cast<int64_t>(a.span))
                : (static_cast<int64_t>(a.pos_ref) - static_cast<int64_t>(a.pos_qry));
            const int64_t diag_b = b.is_rev
                ? (static_cast<int64_t>(b.pos_ref) + static_cast<int64_t>(b.pos_qry) + static_cast<int64_t>(b.span))
                : (static_cast<int64_t>(b.pos_ref) - static_cast<int64_t>(b.pos_qry));
            if (diag_a != diag_b) return diag_a < diag_b;

            if (a.pos_ref != b.pos_ref) return a.pos_ref < b.pos_ref;
            return a.pos_qry < b.pos_qry;
        });
}

// ==================================================================
// 辅助函数：按位置排序锚点
// ==================================================================
void sortAnchorsByPosition(Anchors& anchors)
{
    std::sort(anchors.begin(), anchors.end(),
        [](const Anchor& a, const Anchor& b) {
            if (a.rid_ref != b.rid_ref) return a.rid_ref < b.rid_ref;
            if (a.is_rev != b.is_rev)   return a.is_rev < b.is_rev;
            if (a.pos_ref != b.pos_ref) return a.pos_ref < b.pos_ref;
            return a.pos_qry < b.pos_qry;
        });
}

// ==================================================================
// 辅助函数：过滤高频锚点（后过滤；不等价于 minimap2 的 -f/-U）
// ==================================================================
void filterHighFrequencyAnchors(Anchors& anchors, std::size_t max_occ)
{
    if (anchors.empty() || max_occ == 0) return;

    std::unordered_map<hash_t, std::size_t> hash_count;
    for (const auto& anchor : anchors) {
        hash_count[anchor.hash]++;
    }

    // 注意：这属于“后过滤”，语义上不同于 minimap2 的 -f/-U（后者在展开 occurrences 前过滤）。
    auto new_end = std::remove_if(anchors.begin(), anchors.end(),
        [&hash_count, max_occ](const Anchor& anchor) {
            return hash_count[anchor.hash] > max_occ;
        });

    anchors.erase(new_end, anchors.end());
}

} // namespace anchor
