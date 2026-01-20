#include <doctest/doctest.h>

#include <seed.h>

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

// ================================================================
// 这个文件把 extractMinimizer 的性能测试接入 doctest：
// - 默认跳过（避免 CI/普通单测跑很久）
// - 通过环境变量 HALIGN4_RUN_PERF=1 显式启用
// - 测试里不做严格性能断言，只打印耗时/吞吐，用于人工对比实现改动前后表现
// ================================================================

static bool perfEnabled() {
    const char* v = std::getenv("HALIGN4_RUN_PERF");
    return (v != nullptr) && (*v != '\0') && (std::string(v) != "0");
}
static bool shouldSkipPerf() { return !perfEnabled(); }

static std::string makeRandomDna(std::size_t len, std::uint32_t seed)
{
    static constexpr char bases[4] = {'A', 'C', 'G', 'T'};
    std::uint64_t x = seed;

    // 简单的 xorshift64*，比 mt19937 更轻量，足够生成测试数据
    auto next = [&]() {
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        return x * 2685821657736338717ULL;
    };

    std::string s;
    s.resize(len);
    for (std::size_t i = 0; i < len; ++i) {
        s[i] = bases[static_cast<std::size_t>(next() & 3ULL)];
    }
    return s;
}

TEST_SUITE("minimizer" * doctest::skip(shouldSkipPerf()))
{
    TEST_CASE("extractMinimizer - throughput")
    {
        // 可通过环境变量覆盖不同规模，便于你本地调参
        // HALIGN4_MINIMIZER_SEQ_LEN=10000
        // HALIGN4_MINIMIZER_NUM_SEQS=200
        // HALIGN4_MINIMIZER_ROUNDS=5
        auto getenv_u64 = [](const char* name, std::uint64_t defv) {
            if (const char* p = std::getenv(name); p && *p) return static_cast<std::uint64_t>(std::strtoull(p, nullptr, 10));
            return defv;
        };

        const std::size_t seq_len  = static_cast<std::size_t>(getenv_u64("HALIGN4_MINIMIZER_SEQ_LEN",  30000));
        const std::size_t num_seqs = static_cast<std::size_t>(getenv_u64("HALIGN4_MINIMIZER_NUM_SEQS", 100000));
        const std::size_t rounds   = static_cast<std::size_t>(getenv_u64("HALIGN4_MINIMIZER_ROUNDS",   1));

        // 参数：k/w
        const std::size_t k = static_cast<std::size_t>(getenv_u64("HALIGN4_MINIMIZER_K", 15));
        const std::size_t w = static_cast<std::size_t>(getenv_u64("HALIGN4_MINIMIZER_W", 10));

        MESSAGE("seq_len=" << seq_len << " num_seqs=" << num_seqs << " rounds=" << rounds << " k=" << k << " w=" << w);

        std::vector<std::string> seqs;
        seqs.reserve(num_seqs);
        for (std::size_t i = 0; i < num_seqs; ++i) {
            seqs.emplace_back(makeRandomDna(seq_len, static_cast<std::uint32_t>(1234 + i)));
        }

        // 预热：避免第一次运行的 cache/branch predictor 影响过大
        std::uint64_t checksum = 0;
        for (const auto& s : seqs) {
            auto mz = minimizer::extractMinimizer(s, k, w, false);
            checksum += static_cast<std::uint64_t>(mz.size());
            if (!mz.empty()) checksum ^= mz.front().hash();
        }

        const auto t0 = std::chrono::steady_clock::now();

        std::uint64_t total_minimizers = 0;
        for (std::size_t r = 0; r < rounds; ++r) {
            for (const auto& s : seqs) {
                auto mz = minimizer::extractMinimizer(s, k, w, false);
                total_minimizers += static_cast<std::uint64_t>(mz.size());
                if (!mz.empty()) checksum ^= mz.back().hash();
            }
        }

        const auto t1 = std::chrono::steady_clock::now();
        const double sec = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0).count();

        const double total_bp = static_cast<double>(seq_len) * static_cast<double>(num_seqs) * static_cast<double>(rounds);
        const double bp_per_s = total_bp / sec;
        const double seq_per_s = (static_cast<double>(num_seqs) * static_cast<double>(rounds)) / sec;

        MESSAGE("elapsed_s=" << sec);
        MESSAGE("throughput_bp_per_s=" << bp_per_s);
        MESSAGE("throughput_seq_per_s=" << seq_per_s);
        MESSAGE("total_minimizers=" << total_minimizers);
        MESSAGE("checksum=" << checksum);

        // 不对性能做硬性断言，避免机器差异/负载导致 CI 不稳定。
        CHECK(sec > 0.0);
    }
}
