#ifndef HALIGN4_ANCHOR_H
#define HALIGN4_ANCHOR_H

#include <cstdint>
#include <type_traits>
#include <string_view>
#include <stdexcept>
#include <utility>
#include <vector>
#include <string>
#include <cstddef>
#include <algorithm>
#include <limits>
#include <cmath>
#include "hash.h"

// ================================================================
// anchor 命名空间：锚点（Anchor）相关的数据结构与工具函数
// ================================================================
// 设计动机：
// - Anchor 及其配套的过滤/排序函数，本质属于“锚点/链化阶段”的公共数据结构，
//   不应该和 seed 抽象接口（SeedHitBase/traits）强耦合在同一个命名空间里。
// - 将其拆到 anchor::，可以让 seed:: 只负责“seed/hit 的抽象与提取”，
//   anchor:: 专注“如何用 hits 形成 anchors、如何排序/过滤”。
//
// 重要：
// - 这里仍然使用全局 hash_t（来自 include/hash.h），保持与现有代码一致。
// - 目前 minimizer::collect_anchors 返回 anchor::Anchors；后续如果支持 syncmer/strobemer，
//   也可以复用同一套 anchor 工具函数。
// ================================================================
namespace anchor
{
    // ------------------------------------------------------------------
    // 结构体：Anchor
    // ------------------------------------------------------------------
    // 语义：描述 ref 与 query 在某个 seed/hash 上的一次“锚点匹配”。
    // 下游 chaining 会把一组 anchors 串成一条链（候选比对区域）。
    // ------------------------------------------------------------------
    struct Anchor
    {
        hash_t hash{};              // seed hash（期望 ref/query 相同）
        std::uint32_t rid_ref{};    // ref 序列 id
        std::uint32_t pos_ref{};    // ref 上位置（0-based）
        std::uint32_t rid_qry{};    // query 序列 id（很多场景固定 0）
        std::uint32_t pos_qry{};    // query 上位置（0-based, forward 坐标系）
        std::uint32_t span{};       // 覆盖长度（可用 min(ref.span, qry.span)）
        bool is_rev{};              // ref/query 是否为"相反链"
        // 可选：预先缓存对角线，chaining 常用
        // int32_t diag{}; // (int32_t)pos_ref - (int32_t)pos_qry
    };

    using Anchors = std::vector<Anchor>;

    // ==================================================================
    // 内部辅助结构：用于 ref_hits 的 hash 索引
    // ==================================================================
    struct HashIndex {
        std::size_t start;  // 在排序后数组中的起始索引
        std::size_t count;  // 具有该 hash 的元素数量
    };

    // ==================================================================
    // minimap2 风格的 seeding 过滤参数（默认值与 minimap2 CLI 相同）
    // ==================================================================
    struct SeedFilterParams {
        double f_top_frac = 2e-4;                 // -f
        std::size_t u_floor = 10;                 // -U lower
        std::size_t u_ceil  = 1000000;            // -U upper
        double q_occ_frac  = 0.01;                // --q-occ-frac
        std::size_t sample_every_bp = 500;        // -e
    };

    static inline SeedFilterParams default_mm2_params() {
        return SeedFilterParams{};
    }

    // 计算 -f (fraction) 对应的 occurrence 阈值：忽略 top f_top_frac 最频繁 minimizers
    // 返回值：occ_cutoff（>=1）。当 distinct minimizers 很少或 f_top_frac==0 时，返回 +inf。
    std::size_t compute_occ_cutoff_top_frac(const std::vector<std::size_t>& occs,
                                            double f_top_frac);

    // 计算最终 reference occurrence 阈值：max{u_floor, min{u_ceil, -f}}
    std::size_t compute_ref_occ_threshold(const std::vector<std::size_t>& occs,
                                          const SeedFilterParams& p);

    // =====================================================================
    // sortAnchorsByDiagonal - 按对角线排序锚点（用于链化算法）
    // =====================================================================
    void sortAnchorsByDiagonal(Anchors& anchors);

    // =====================================================================
    // sortAnchorsByPosition - 按位置排序锚点
    // =====================================================================
    void sortAnchorsByPosition(Anchors& anchors);

    // =====================================================================
    // filterHighFrequencyAnchors - 过滤高频锚点（参考 minimap2）
    // =====================================================================
    void filterHighFrequencyAnchors(Anchors& anchors, std::size_t max_occ = 500);

} // namespace anchor


#endif //HALIGN4_ANCHOR_H