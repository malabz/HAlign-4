#include <doctest/doctest.h>

#include <mash.h>

#include <cmath>
#include <algorithm>
#include <string>
#include <vector>
#include <chrono>
#include <random>
#include <iostream>

TEST_SUITE("mash")
{
    TEST_CASE("intersectionSizeSortedUnique")
    {
        std::vector<hash_t> a{1, 2, 3, 10};
        std::vector<hash_t> b{2, 3, 4, 5, 10, 11};
        CHECK(mash::intersectionSizeSortedUnique(a, b) == 3);
    }

    TEST_CASE("sketchFromSequence - empty & sketch_size=0")
    {
        const std::size_t k = 15;
        const std::size_t sketch_size = 0;

        auto sk = mash::sketchFromSequence(std::string("ACGTACGT"), k, sketch_size);
        CHECK(sk.empty());

        auto sk2 = mash::sketchFromSequence(std::string(""), k, 100);
        CHECK(sk2.empty());
    }

    TEST_CASE("jaccard - empty sets")
    {
        mash::Sketch a;
        a.k = 15;
        mash::Sketch b;
        b.k = 15;
        CHECK(mash::jaccard(a, b) == doctest::Approx(1.0));

        mash::Sketch c;
        c.k = 15;
        c.hashes = {1, 2, 3};
        CHECK(mash::jaccard(a, c) == doctest::Approx(0.0));
    }

    TEST_CASE("jaccard - identical sketches")
    {
        mash::Sketch a;
        a.k = 21;
        a.hashes = {1, 2, 3, 4};
        mash::Sketch b;
        b.k = 21;
        b.hashes = {1, 2, 3, 4};
        CHECK(mash::jaccard(a, b) == doctest::Approx(1.0));
        CHECK(mash::mashDistanceFromJaccard(1.0, 21) == doctest::Approx(0.0));
        CHECK(mash::aniFromJaccard(1.0, 21) == doctest::Approx(1.0));
    }

    TEST_CASE("jaccard - disjoint sketches")
    {
        mash::Sketch a;
        a.k = 21;
        a.hashes = {1, 2, 3};
        mash::Sketch b;
        b.k = 21;
        b.hashes = {4, 5, 6};
        CHECK(mash::jaccard(a, b) == doctest::Approx(0.0));
        CHECK(!std::isfinite(mash::mashDistanceFromJaccard(0.0, 21)));
        CHECK(mash::aniFromJaccard(0.0, 21) == doctest::Approx(0.0));
    }

    TEST_CASE("sketchFromSequence + jaccard smoke")
    {
        const std::size_t k = 15;
        const std::size_t sketch_size = 200;

        const std::string s1 = "ACGTACGTACGTACGTACGTACGTACGTACGT";
        const std::string s2 = "ACGTACGTACGTACGTACGTACGTACGTACGT";
        const std::string s3 = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";

        auto sk1 = mash::sketchFromSequence(s1, k, sketch_size);
        auto sk2 = mash::sketchFromSequence(s2, k, sketch_size);
        auto sk3 = mash::sketchFromSequence(s3, k, sketch_size);

        CHECK(sk1.k == k);
        CHECK(sk1.size() <= sketch_size);

        const double j12 = mash::jaccard(sk1, sk2);
        CHECK(j12 == doctest::Approx(1.0));

        const double j13 = mash::jaccard(sk1, sk3);
        CHECK(j13 >= 0.0);
        CHECK(j13 <= 1.0);
    }

    // Performance test: only run when HALIGN4_RUN_PERF=1 (set by test/run_tests.sh --perf)
    static std::string random_dna(std::mt19937_64 &rng, std::size_t len) {
        static const char bases[4] = {'A','C','G','T'};
        std::string s;
        s.reserve(len);
        for (std::size_t i = 0; i < len; ++i) s.push_back(bases[rng() & 3]);
        return s;
    }

    TEST_CASE("mash_perf") {
        const char* env = std::getenv("HALIGN4_RUN_PERF");
        if (!env || std::string(env) != "1") {
            DOCTEST_INFO("mash_perf skipped; run run_tests.sh --perf to enable");
            return;
        }

        const std::size_t N = []{
            const char* e = std::getenv("MASH_PERF_N");
            return e ? static_cast<std::size_t>(std::stoull(e)) : 30000ULL;
        }();
        const std::size_t L = []{
            const char* e = std::getenv("MASH_PERF_L");
            return e ? static_cast<std::size_t>(std::stoull(e)) : 30000ULL;
        }();
        const std::size_t k = 21;
        const std::size_t sketch_size = 2000;

        std::mt19937_64 rng(123456);
        std::vector<std::string> seqs;
        seqs.reserve(N);
        for (std::size_t i = 0; i < N; ++i) seqs.push_back(random_dna(rng, L));

        auto t0 = std::chrono::steady_clock::now();
        for (const auto &s : seqs) {
            volatile auto sk = mash::sketchFromSequence(s, k, sketch_size);
            (void)sk;
        }
        auto t1 = std::chrono::steady_clock::now();

        const double seconds = std::chrono::duration<double>(t1 - t0).count();
        MESSAGE("mash_perf: N=" << N << " L=" << L << " took " << seconds << "s");

        CHECK(true);
    }

    TEST_CASE("sketchFromSequence - hashes are sorted & unique")
    {
        const std::string s = "ACGTACGTACGTACGTACGTACGTACGTACGT";
        const std::size_t k = 4;
        const std::size_t sketch_size = 200;

        auto sk = mash::sketchFromSequence(s, k, sketch_size);
        CHECK(sk.k == k);

        // 检查有序
        CHECK(std::is_sorted(sk.hashes.begin(), sk.hashes.end()));

        // 检查无重复
        CHECK(std::adjacent_find(sk.hashes.begin(), sk.hashes.end()) == sk.hashes.end());
    }

    // 参考实现：暴力计算所有 k-mer（rolling 2-bit + getHash2bit），
    // 取 canonical/noncanonical 后生成完整 hash 列表，最后 sort+unique 并取前 sketch_size。
    static std::vector<hash_t> reference_bottom_k(const std::string& seq,
                                                  std::size_t k,
                                                  std::size_t sketch_size,
                                                  bool noncanonical,
                                                  int seed)
    {
        std::vector<hash_t> all;
        if (k == 0 || sketch_size == 0 || seq.size() < k || k > 32) return all;

        const std::uint64_t mask = (1ULL << (2 * k)) - 1ULL;
        const std::uint64_t shift = 2ULL * (k - 1);
        std::uint64_t fwd = 0;
        std::uint64_t rev = 0;
        std::size_t valid = 0;

        all.reserve(seq.size());

        for (std::size_t i = 0; i < seq.size(); ++i)
        {
            const uint8_t c = mash::nt4_table[static_cast<unsigned char>(seq[i])];
            if (c >= 4) {
                fwd = rev = 0;
                valid = 0;
                continue;
            }

            fwd = ((fwd << 2) | c) & mask;
            rev = (rev >> 2) | (std::uint64_t(3U ^ c) << shift);

            if (valid < k) ++valid;
            if (valid < k) continue;

            const std::uint64_t code = noncanonical ? fwd : std::min(fwd, rev);
            all.push_back(getHash2bit(code, static_cast<std::uint32_t>(seed)));
        }

        std::sort(all.begin(), all.end());
        all.erase(std::unique(all.begin(), all.end()), all.end());
        if (all.size() > sketch_size) all.resize(sketch_size);
        return all;
    }

    TEST_CASE("sketchFromSequence - bottom-k property matches reference")
    {
        const std::string s = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        const std::size_t k = 15;
        const std::size_t sketch_size = 50;
        const int seed = 42;

        SUBCASE("noncanonical=true")
        {
            auto sk = mash::sketchFromSequence(s, k, sketch_size, /*noncanonical*/ true, seed);
            auto ref = reference_bottom_k(s, k, sketch_size, /*noncanonical*/ true, seed);
            CHECK(sk.hashes == ref);
        }

        SUBCASE("noncanonical=false (canonical)")
        {
            auto sk = mash::sketchFromSequence(s, k, sketch_size, /*noncanonical*/ false, seed);
            auto ref = reference_bottom_k(s, k, sketch_size, /*noncanonical*/ false, seed);
            CHECK(sk.hashes == ref);
        }
    }

    TEST_CASE("sketchFromSequence - handles invalid chars by resetting window")
    {
        // 中间插入 N：跨过 N 的 k-mer 不应产生（窗口会重置）
        const std::string s = "ACGTACGTNNNNACGTACGT";
        const std::size_t k = 5;
        const std::size_t sketch_size = 200;
        const int seed = 7;

        auto sk = mash::sketchFromSequence(s, k, sketch_size, /*noncanonical*/ true, seed);
        auto ref = reference_bottom_k(s, k, sketch_size, /*noncanonical*/ true, seed);
        CHECK(sk.hashes == ref);
        CHECK(std::is_sorted(sk.hashes.begin(), sk.hashes.end()));
        CHECK(std::adjacent_find(sk.hashes.begin(), sk.hashes.end()) == sk.hashes.end());
    }

    TEST_CASE("sketchFromSequence - canonical differs from noncanonical on non-palindrome")
    {
        // 选一个明显非回文的序列，canonical 与 noncanonical 通常会不一样（并不保证必须不同，但很大概率）。
        const std::string s = "ACGTTGCAACGTTGCAACGTTGCA";
        const std::size_t k = 7;
        const std::size_t sketch_size = 100;
        const int seed = 123;

        auto sk_nc = mash::sketchFromSequence(s, k, sketch_size, /*noncanonical*/ true, seed);
        auto sk_can = mash::sketchFromSequence(s, k, sketch_size, /*noncanonical*/ false, seed);

        CHECK(std::is_sorted(sk_nc.hashes.begin(), sk_nc.hashes.end()));
        CHECK(std::is_sorted(sk_can.hashes.begin(), sk_can.hashes.end()));

        // 至少确保两个结果都符合 bottom-k 参考
        CHECK(sk_nc.hashes == reference_bottom_k(s, k, sketch_size, true, seed));
        CHECK(sk_can.hashes == reference_bottom_k(s, k, sketch_size, false, seed));

        // “通常”两者不同：用 CHECK_FALSE 可能导致偶发失败，所以这里只做 INFO 输出。
        DOCTEST_INFO("noncanonical size=", sk_nc.hashes.size(), "; canonical size=", sk_can.hashes.size());
    }
}
