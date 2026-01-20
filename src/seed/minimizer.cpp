#include "seed.h"

#include <cstdint>
#include <string>
#include <vector>

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

} // namespace minimizer

