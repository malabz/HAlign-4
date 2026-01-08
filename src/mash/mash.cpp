#include "mash.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include "MurmurHash3.h"

namespace mash
{
    static inline std::uint8_t nt4(unsigned char c) noexcept
    {
        switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        case 'U': case 'u': return 3;
        default: return 4;
        }
    }

    static inline double clamp01(double x) noexcept
    {
        if (x < 0.0) return 0.0;
        if (x > 1.0) return 1.0;
        return x;
    }

    std::size_t intersectionSizeSortedUnique(const std::vector<hash_t>& a,
                                             const std::vector<hash_t>& b) noexcept
    {
        std::size_t i = 0, j = 0, inter = 0;
        while (i < a.size() && j < b.size()) {
            const auto av = a[i];
            const auto bv = b[j];
            if (av == bv) {
                ++inter;
                ++i;
                ++j;
            } else if (av < bv) {
                ++i;
            } else {
                ++j;
            }
        }
        return inter;
    }

    static inline hash_t murmur64(const void* data, int len, std::uint32_t seed) noexcept
    {
        std::uint64_t out[2] = {0, 0};
        MurmurHash3_x64_128(data, len, seed, out);
        return static_cast<hash_t>(out[0]);
    }

    // canonical=false: 使用 is_forward 决定方向
    // canonical=true : 无视 is_forward，对每个 k-mer 取 min(fwd, revcomp)
    // 这里为了接口兼容，is_forward=true 表示 forward-only；is_forward=false 表示 revcomp-only。
    Sketch sketchFromSequence(const std::string& seq,
                             std::size_t k,
                             std::size_t /*w*/,
                             std::size_t sketch_size,
                             bool is_forward)
    {
        Sketch sk;
        sk.k = k;

        if (k == 0 || sketch_size == 0) return sk;
        if (seq.size() < k) return sk;

        // 简化：仅支持 k<=31 的 2-bit rolling（与 Mash 默认 DNA k 限制类似）
        if (k > 31) return sk;

        const std::uint64_t mask = (1ULL << (2 * k)) - 1ULL;
        const std::uint64_t shift = 2ULL * (k - 1);
        std::uint64_t fwd = 0;
        std::uint64_t rev = 0;
        std::size_t valid = 0;

        // 使用 Mash 的思路：对 k-mer 字符串做 MurmurHash3
        // 这里不分配 substr；用 rolling code 的 8-byte 表示作为 key（更快/更少内存）。
        const std::uint32_t seed = 42; // 轻量实现固定 seed；如需可后续放到参数里

        sk.hashes.reserve(std::min<std::size_t>(sketch_size, seq.size() - k + 1));

        for (std::size_t i = 0; i < seq.size(); ++i) {
            const std::uint8_t c = nt4(static_cast<unsigned char>(seq[i]));
            if (c >= 4) {
                fwd = rev = 0;
                valid = 0;
                continue;
            }

            fwd = ((fwd << 2) | c) & mask;
            rev = (rev >> 2) | (std::uint64_t(3U ^ c) << shift);

            if (valid < k) ++valid;
            if (valid < k) continue;

            const std::uint64_t code = is_forward ? fwd : rev;
            const hash_t h = murmur64(&code, static_cast<int>(sizeof(code)), seed);
            sk.hashes.push_back(h);
        }

        std::sort(sk.hashes.begin(), sk.hashes.end());
        sk.hashes.erase(std::unique(sk.hashes.begin(), sk.hashes.end()), sk.hashes.end());
        if (sk.hashes.size() > sketch_size) sk.hashes.resize(sketch_size);

        return sk;
    }

    double jaccard(const Sketch& a, const Sketch& b)
    {
        if (a.k != b.k) throw std::invalid_argument("mash::jaccard: mismatched k");

        if (a.hashes.empty() && b.hashes.empty()) return 1.0;
        if (a.hashes.empty() || b.hashes.empty()) return 0.0;

        const std::size_t inter = intersectionSizeSortedUnique(a.hashes, b.hashes);
        const std::size_t uni = a.hashes.size() + b.hashes.size() - inter;
        if (uni == 0) return 1.0;
        return static_cast<double>(inter) / static_cast<double>(uni);
    }

    double mashDistanceFromJaccard(double j, std::size_t k)
    {
        if (k == 0) throw std::invalid_argument("mash::mashDistanceFromJaccard: k must be > 0");
        if (!(j > 0.0)) return std::numeric_limits<double>::infinity();
        if (j >= 1.0) return 0.0;

        const double x = (2.0 * j) / (1.0 + j);
        if (!(x > 0.0)) return std::numeric_limits<double>::infinity();
        return -std::log(x) / static_cast<double>(k);
    }

    double aniFromJaccard(double j, std::size_t k)
    {
        if (k == 0) throw std::invalid_argument("mash::aniFromJaccard: k must be > 0");
        if (!(j > 0.0)) return 0.0;
        if (j >= 1.0) return 1.0;

        const double x = (2.0 * j) / (1.0 + j);
        if (!(x > 0.0)) return 0.0;

        return clamp01(std::pow(x, 1.0 / static_cast<double>(k)));
    }

    double aniFromMashDistance(double d)
    {
        if (!std::isfinite(d)) return 0.0;
        if (d <= 0.0) return 1.0;
        return clamp01(std::exp(-d));
    }

} // namespace mash
