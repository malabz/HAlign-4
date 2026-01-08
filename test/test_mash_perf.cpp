#include <doctest/doctest.h>
#include <mash.h>

#include <chrono>
#include <random>
#include <string>
#include <vector>
#include <iostream>

// This perf test is guarded by the environment variable RUN_MASH_PERF.
// If not set to "1", the test will be skipped quickly so normal CI/tests stay fast.
// When enabled, it will generate N random DNA sequences of length L and measure time
// to compute sketches for each sequence.

static std::string random_dna(std::mt19937_64 &rng, std::size_t len) {
    static const char bases[4] = {'A','C','G','T'};
    std::string s;
    s.reserve(len);
    for (std::size_t i = 0; i < len; ++i) s.push_back(bases[rng() & 3]);
    return s;
}

TEST_CASE("mash_perf") {
    const char* env = std::getenv("RUN_MASH_PERF");
    if (!env || std::string(env) != "1") {
        DOCTEST_INFO("mash_perf skipped; set RUN_MASH_PERF=1 to enable");
        return;
    }

    // Default parameters for heavy run. You can override via env vars for convenience.
    const std::size_t N = []{
        const char* e = std::getenv("MASH_PERF_N");
        return e ? static_cast<std::size_t>(std::stoull(e)) : 10000ULL; // 10k sequences
    }();
    const std::size_t L = []{
        const char* e = std::getenv("MASH_PERF_L");
        return e ? static_cast<std::size_t>(std::stoull(e)) : 30000ULL; // length 30k
    }();
    const std::size_t k = 31;
    const std::size_t w = 10;
    const std::size_t sketch_size = 200;

    std::mt19937_64 rng(123456);

    std::vector<std::string> seqs;
    seqs.reserve(N);
    for (std::size_t i = 0; i < N; ++i) seqs.push_back(random_dna(rng, L));

    auto t0 = std::chrono::steady_clock::now();
    for (const auto &s : seqs) {
        volatile auto sk = mash::sketchFromSequence(s, k, w, sketch_size);
        (void)sk;
    }
    auto t1 = std::chrono::steady_clock::now();

    const double seconds = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "mash_perf: N=" << N << " L=" << L << " took " << seconds << "s\n";

    // Keep doctest happy
    CHECK(true);
}
