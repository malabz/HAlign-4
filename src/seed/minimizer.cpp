#include "seed.h"

#include <cstdint>
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <cmath>
#include <limits>
#include "anchor.h"

namespace minimizer
{
    // splitmix64：非常常用的 64-bit mixer，速度快且分布不错。
    // 在 minimap2 里也会对 k-mer 编码做 hash64 混洗，理念相同。
    static inline constexpr std::uint64_t splitmix64(std::uint64_t x) noexcept
    {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }

    // Cand：窗口中的候选 (hash, pos)
    // pos 是 k-mer 的起始位置（0-based）
    struct Cand
    {
        std::uint64_t h;
        std::uint32_t pos;
    };

    // =========================================================
    // RingMinQueue：固定容量的“环形单调队列”（维护窗口最小值）
    // =========================================================
    class RingMinQueue
    {
    public:
        explicit RingMinQueue(std::uint32_t capacity)
            : buf_(capacity), cap_(capacity)
        {
        }

        void clear() noexcept
        {
            head_ = 0;
            size_ = 0;
        }

        bool empty() const noexcept { return size_ == 0; }

        void push(std::uint64_t h, std::uint32_t pos) noexcept
        {
            const Cand c{h, pos};
            // 弹出所有 >= 当前 hash 的队尾元素，保持单调递增。
            while (size_ && back().h >= c.h) pop_back();
            buf_[idx(size_)] = c;
            ++size_;
        }

        void popExpired(std::uint32_t win_start) noexcept
        {
            while (size_ && front().pos < win_start) pop_front();
        }

        std::uint64_t minHash() const noexcept { return front().h; }
        std::uint32_t minPos() const noexcept { return front().pos; }

    private:
        std::vector<Cand> buf_;
        std::uint32_t cap_{0};
        std::uint32_t head_{0};
        std::uint32_t size_{0};

        std::uint32_t idx(std::uint32_t off) const noexcept
        {
            std::uint32_t i = head_ + off;
            if (i >= cap_) i -= cap_;
            return i;
        }

        Cand& front() noexcept { return buf_[head_]; }
        const Cand& front() const noexcept { return buf_[head_]; }

        Cand& back() noexcept { return buf_[idx(size_ - 1)]; }
        const Cand& back() const noexcept { return buf_[idx(size_ - 1)]; }

        void pop_front() noexcept
        {
            head_ = idx(1);
            --size_;
        }

        void pop_back() noexcept
        {
            --size_;
        }
    };


    // =============================================================
    // 目标：从单条序列中提取 minimizer（带位置 hit）。
    // =============================================================
    MinimizerHits extractMinimizer(const std::string& seq,
                                   std::size_t k,
                                   std::size_t w,
                                   bool non_canonical)
    {
        MinimizerHits out;

        const std::uint32_t n = static_cast<std::uint32_t>(seq.size());
        if (k == 0 || w == 0 || n < k) return out;
        if (k > 31) return out;
        if (w >= 256) return out;

        const auto& nt4_table_ref = minimizer::nt4_table;

        const std::uint32_t total_kmer = n - static_cast<std::uint32_t>(k) + 1;
        const std::uint32_t win = std::min<std::uint32_t>(static_cast<std::uint32_t>(w), total_kmer);
        if (win == 0) return out;

        out.reserve(std::max<std::uint32_t>(1, n / win));

        const std::uint64_t mask = (1ULL << (2 * k)) - 1ULL;
        const std::uint64_t shift = 2ULL * (k - 1);

        std::uint64_t fwd = 0;
        std::uint64_t rev = 0;
        std::uint32_t valid = 0;

        RingMinQueue q(win);

        Cand last_out{0, 0};
        bool has_last = false;

        for (std::uint32_t i = 0; i < n; ++i) {
            const std::uint8_t c = nt4_table_ref[static_cast<std::uint8_t>(seq[i])];
            if (c >= 4) {
                fwd = rev = 0;
                valid = 0;
                q.clear();
                has_last = false;
                continue;
            }

            fwd = ((fwd << 2) | c) & mask;
            rev = (rev >> 2) | (std::uint64_t(3U ^ c) << shift);

            if (valid < k) ++valid;
            if (valid < k) continue;

            const std::uint32_t pos = i + 1 - static_cast<std::uint32_t>(k);

            const std::uint64_t code = non_canonical ? fwd : std::min(fwd, rev);
            const std::uint64_t h64 = splitmix64(code);
            const std::uint64_t h56 = (h64 >> 8); // 取高 56bit，留出低 8bit 给 span

            q.push(h56, pos);

            if (pos + 1 < win) continue;
            const std::uint32_t win_start = pos + 1 - win;
            q.popExpired(win_start);

            if (!q.empty()) {
                const Cand cur{q.minHash(), q.minPos()};
                if (!has_last || cur.h != last_out.h || cur.pos != last_out.pos) {
                    out.emplace_back(cur.h, cur.pos,
                                     /*rid*/ 0,
                                     non_canonical ? true : (fwd <= rev),
                                     static_cast<std::uint8_t>(k));
                    last_out = cur;
                    has_last = true;
                }
            }
        }

        return out;
    }


// ------------------------------------------------------------------
// 函数：collect_anchors
// ------------------------------------------------------------------
// 功能：
// 从 ref_hits 和 qry_hits 中收集锚点（Anchor）列表。
//
// 这里的“锚点”含义：
// - 一个锚点表示：ref 与 query 在某个 minimizer hash 上发生了匹配；
// - 因为同一个 hash 在参考/查询中可能出现多次，所以一个 qry_hit 可能会展开为多个 anchors。
// - anchors 是后续 chaining（链化）的输入。
//
// 设计与 minimap2 对齐的关键点（非常重要）：
// 1) 先统计/过滤，再展开（occurrence expansion）
//    - minimap2 的 -f/-U/--q-occ-frac/-e 等策略，本质是在“展开 occurrences 前”抑制重复区域。
//    - 如果不这么做，重复区域会让一个高频 minimizer 展开成 O(occ_ref * occ_qry) 个 anchors，
//      直接导致内存/时间爆炸。
//
// 2) 本实现的过滤参数来自 anchor::SeedFilterParams（默认值模仿 minimap2 CLI 默认）：
//    - f_top_frac（-f）：忽略参考端最频繁的 top fraction minimizers（按 distinct 计数）
//    - u_floor/u_ceil（-U）：对 occurrence 阈值做上下限夹逼
//    - q_occ_frac（--q-occ-frac）：query 端过于高频（并且高于 reference 阈值）则丢弃
//    - sample_every_bp（-e）：对高频 minimizer 做位置稀疏采样（而不是全部展开）
//
// 3) 输出 anchors 未排序：
//    - minimap2 后续会按 (rid, strand, diagonal, ref_pos, qry_pos) 排序再进行 DP chaining。
//    - 这里仅负责收集，排序由调用方（anchor::sortAnchorsByDiagonal 等）完成。
//
// 输入：
// - ref_hits：reference 的 minimizer hits（可能来自全参考、或某段参考）
// - qry_hits：query 的 minimizer hits
//
// 输出：
// - anchor::Anchors：Anchor 列表（每个元素都记录 ref/qry 的 rid/pos/span/is_rev）
//
// 复杂度：
// - 排序 ref_hits：O(R log R)
// - 统计 qry occurrence：O(Q)
// - 生成 anchors：O(Q * avg_occ_ref_for_hash)
//   （重复区域会被过滤/稀疏，避免最坏情况爆炸）
// ------------------------------------------------------------------
anchor::Anchors collect_anchors(const MinimizerHits& ref_hits, const MinimizerHits& qry_hits, anchor::SeedFilterParams params)
{
    anchor::Anchors anchors;

    // 边界条件：任意一侧为空则不可能产生锚点
    if (ref_hits.empty() || qry_hits.empty()) {
        return anchors;
    }


    // ------------------------------------------------------------------
    // Step 1：对 ref_hits 排序 + 构建 hash -> (start, count) 索引
    // ------------------------------------------------------------------
    // 说明：
    // - ref_hits 里同一个 hash 可能出现多次（不同位置/不同 rid）。
    // - 将其按 (hash, rid, pos, strand) 排序后，同一 hash 的 hits 会变成连续区间。
    // - 我们用一个 unordered_map 保存每个 hash 对应的连续区间 (start, count)，
    //   这样 qry 端查到一个 hash 后，就可以 O(occ_ref) 展开。
    std::vector<minimizer::MinimizerHit> sorted_ref = ref_hits;
    std::sort(sorted_ref.begin(), sorted_ref.end()); // 使用 SeedHitBase 的 operator<

    std::unordered_map<hash_t, anchor::HashIndex> hash_index;
    hash_index.reserve(sorted_ref.size());

    // ref_occs：记录每个 distinct hash 在 reference 的出现次数（occurrence），
    // 用于后续根据 -f（top fraction）估算过滤阈值。
    std::vector<std::size_t> ref_occs;
    ref_occs.reserve(sorted_ref.size() / 2 + 1);

    if (!sorted_ref.empty()) {
        std::size_t start = 0;
        hash_t current_hash = sorted_ref[0].hash();

        for (std::size_t i = 1; i <= sorted_ref.size(); ++i) {
            // 遇到 hash 变化或到达末尾 => 结算上一段 hash 的区间
            if (i == sorted_ref.size() || sorted_ref[i].hash() != current_hash) {
                const std::size_t occ = i - start;
                hash_index[current_hash] = anchor::HashIndex{start, occ};
                ref_occs.push_back(occ);
                if (i < sorted_ref.size()) {
                    start = i;
                    current_hash = sorted_ref[i].hash();
                }
            }
        }
    }

    // ------------------------------------------------------------------
    // Step 2：根据 ref_occs + (-f/-U) 计算“参考端高频阈值”
    // ------------------------------------------------------------------
    // 该阈值用于判断某个 hash 是否属于重复区域（高频）。
    // 注意：这是“过滤发生在展开 occurrences 前”的关键点。
    const std::size_t ref_occ_thr = anchor::compute_ref_occ_threshold(ref_occs, params);

    // ------------------------------------------------------------------
    // Step 3：统计 query 端每个 hash 的 occurrence（用于 --q-occ-frac）
    // ------------------------------------------------------------------
    // minimap2 的直觉：
    // - 如果某个 hash 在 query 里也极其高频，那么它很可能来自低复杂度/重复区域
    // - 这些 seeds 对定位帮助不大，却会产生大量 anchors
    std::unordered_map<hash_t, std::size_t> qry_occ;
    qry_occ.reserve(qry_hits.size());
    for (const auto& qh : qry_hits) {
        ++qry_occ[qh.hash()];
    }

    // q_occ_limit：把 --q-occ-frac 从“比例”变成“数量阈值”
    const double q_occ_limit = params.q_occ_frac > 0.0
        ? (params.q_occ_frac * static_cast<double>(qry_hits.size()))
        : std::numeric_limits<double>::infinity();

    // 预估：通常 anchors 数量与 qry_hits 同量级（重复会被压制）
    anchors.reserve(qry_hits.size());

    // ------------------------------------------------------------------
    // Step 4：遍历 qry_hits，查 ref 索引并生成 anchors
    // ------------------------------------------------------------------
    for (const auto& qry_hit : qry_hits) {
        const hash_t qry_hash = qry_hit.hash();

        // ref 端不存在该 hash => 无法形成锚点
        auto it = hash_index.find(qry_hash);
        if (it == hash_index.end()) continue;

        const anchor::HashIndex& idx = it->second;
        const std::size_t ref_occ = idx.count;

        // ---- 4.1 --q-occ-frac：query 端过于高频时直接丢弃（降低爆炸风险）
        // minimap2 语义：当 query 端某个 hash 出现次数超过阈值时，直接丢弃该 hash
        // 这是独立于 ref 端过滤的机制，用于应对 query 端的低复杂度区域
        if (params.q_occ_frac > 0.0) {
            const std::size_t qocc = qry_occ[qry_hash];
            if (static_cast<double>(qocc) > q_occ_limit) {
                continue;
            }
        }

        // ---- 4.2 -f/-U + -e：reference 端高频 minimizer 采用稀疏采样
        // 语义：
        // - 如果 ref_occ > ref_occ_thr，说明该 hash 在参考中非常“重复”。
        // - minimap2 会对高频 minimizer 进行稀疏/筛选，避免把所有 occurrences 展开。
        // - 我们这里用最简单的“按 query 位置取模采样”的方式来近似 -e 行为：
        //   仅当 qry_hit.pos % sample_every_bp == 0 时才展开。
        //
        // 注意：这不是唯一实现方式，但可保持“热点路径不爆炸”的关键特性。
        if (ref_occ > ref_occ_thr) {
            if (params.sample_every_bp == 0) continue;
            if ((static_cast<std::size_t>(qry_hit.pos()) % params.sample_every_bp) != 0) {
                continue;
            }
        }

        // ---- 4.3 展开 reference 的 occurrences，生成 anchors
        // 生成时只做字段拷贝/轻量计算，避免额外分配。
        for (std::size_t i = idx.start; i < idx.start + idx.count; ++i) {
            const auto& ref_hit = sorted_ref[i];

            anchor::Anchor anchor;
            anchor.hash = qry_hash;
            anchor.rid_ref = ref_hit.rid();
            anchor.pos_ref = ref_hit.pos();
            anchor.rid_qry = qry_hit.rid();
            anchor.pos_qry = qry_hit.pos();

            // span 取二者较小值：保守估计“可用匹配长度”
            // （下游 chaining 的 gap/penalty 模型通常只需要一个 span 量级）
            anchor.span = std::min(ref_hit.span(), qry_hit.span());

            // 方向：ref XOR qry（与 minimap2 的 rev 语义一致）
            anchor.is_rev = (ref_hit.strand() != qry_hit.strand());

            anchors.emplace_back(anchor);
        }
    }

    return anchors;
}

} // namespace minimizer

