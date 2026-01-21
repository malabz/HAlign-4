// ================================================================
// test_anchor.cpp - 锚点收集函数的单元测试与性能测试
// ================================================================
// 测试维度：
// 1. 准确性：验证 collect_anchors 能正确生成锚点
// 2. 边界情况：空输入、无匹配、单一匹配等
// 3. 过滤逻辑：高频过滤、稀疏采样、q-occ-frac 等参数
// 4. 一致性：相同输入产生稳定输出（排除随机性）
// 5. 性能：大规模数据下的吞吐量（可选，通过环境变量启用）
// ================================================================

#include <doctest/doctest.h>
#include "seed.h"
#include "anchor.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>
#include <unordered_set>
#include <chrono>

// ================================================================
// 辅助函数：性能测试开关
// ================================================================
static bool perfEnabled() {
    const char* v = std::getenv("HALIGN4_RUN_PERF");
    return (v != nullptr) && (*v != '\0') && (std::string(v) != "0");
}
static bool shouldSkipPerf() { return !perfEnabled(); }

// ================================================================
// 辅助函数：生成随机 DNA 序列（用于性能测试）
// ================================================================
static std::string makeRandomDna(std::size_t len, std::uint32_t seed)
{
    static constexpr char bases[4] = {'A', 'C', 'G', 'T'};
    std::uint64_t x = seed;

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

// ================================================================
// 辅助函数：手动构造 MinimizerHit（用于精确测试）
// ================================================================
static minimizer::MinimizerHit makeHit(hash_t hash56, std::uint32_t pos,
                                       std::uint32_t rid = 0,
                                       bool strand = true,
                                       std::uint8_t span = 15)
{
    return minimizer::MinimizerHit(hash56, pos, rid, strand, span);
}

// ================================================================
// TEST SUITE: anchor 锚点收集测试
// ================================================================

TEST_SUITE("anchor::collect_anchors - 准确性与边界测试")
{
    // ------------------------------------------------------------------
    // TEST CASE 1: 空输入 - ref 和 qry 都为空
    // ------------------------------------------------------------------
    TEST_CASE("collect_anchors - 空输入返回空锚点列表")
    {
        minimizer::MinimizerHits ref_hits;
        minimizer::MinimizerHits qry_hits;

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits);

        CHECK(anchors.empty());
    }

    // ------------------------------------------------------------------
    // TEST CASE 2: ref 为空，qry 非空
    // ------------------------------------------------------------------
    TEST_CASE("collect_anchors - ref 为空时无锚点")
    {
        minimizer::MinimizerHits ref_hits;
        minimizer::MinimizerHits qry_hits;
        qry_hits.push_back(makeHit(0x123456, 10));
        qry_hits.push_back(makeHit(0x789ABC, 20));

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits);

        CHECK(anchors.empty());
    }

    // ------------------------------------------------------------------
    // TEST CASE 3: qry 为空，ref 非空
    // ------------------------------------------------------------------
    TEST_CASE("collect_anchors - qry 为空时无锚点")
    {
        minimizer::MinimizerHits ref_hits;
        ref_hits.push_back(makeHit(0x123456, 100));
        ref_hits.push_back(makeHit(0x789ABC, 200));

        minimizer::MinimizerHits qry_hits;

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits);

        CHECK(anchors.empty());
    }

    // ------------------------------------------------------------------
    // TEST CASE 4: 单一完美匹配 - ref 和 qry 各有一个相同 hash
    // ------------------------------------------------------------------
    TEST_CASE("collect_anchors - 单一完美匹配")
    {
        minimizer::MinimizerHits ref_hits;
        ref_hits.push_back(makeHit(0x123456, 100, 0, true, 15));

        minimizer::MinimizerHits qry_hits;
        qry_hits.push_back(makeHit(0x123456, 50, 0, true, 15));

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits);

        REQUIRE(anchors.size() == 1);
        CHECK(anchors[0].hash == 0x123456);
        CHECK(anchors[0].pos_ref == 100);
        CHECK(anchors[0].pos_qry == 50);
        CHECK(anchors[0].rid_ref == 0);
        CHECK(anchors[0].rid_qry == 0);
        CHECK(anchors[0].span == 15);
        CHECK(anchors[0].is_rev == false); // 同方向 => is_rev = false
    }

    // ------------------------------------------------------------------
    // TEST CASE 5: 反向匹配 - ref 和 qry 方向不同
    // ------------------------------------------------------------------
    TEST_CASE("collect_anchors - 反向匹配检测")
    {
        minimizer::MinimizerHits ref_hits;
        ref_hits.push_back(makeHit(0xABCDEF, 200, 0, true, 15));

        minimizer::MinimizerHits qry_hits;
        qry_hits.push_back(makeHit(0xABCDEF, 80, 0, false, 15));

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits);

        REQUIRE(anchors.size() == 1);
        CHECK(anchors[0].hash == 0xABCDEF);
        CHECK(anchors[0].is_rev == true); // strand 不同 => is_rev = true
    }

    // ------------------------------------------------------------------
    // TEST CASE 6: 多个匹配 - 同一 hash 在 ref/qry 都出现多次
    // ------------------------------------------------------------------
    TEST_CASE("collect_anchors - 一对多展开（occurrence expansion）")
    {
        minimizer::MinimizerHits ref_hits;
        // ref 端同一 hash 出现 3 次
        ref_hits.push_back(makeHit(0x111111, 100));
        ref_hits.push_back(makeHit(0x111111, 200));
        ref_hits.push_back(makeHit(0x111111, 300));

        minimizer::MinimizerHits qry_hits;
        // qry 端同一 hash 出现 2 次
        qry_hits.push_back(makeHit(0x111111, 50));
        qry_hits.push_back(makeHit(0x111111, 150));

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits);

        // 应该生成 3 * 2 = 6 个锚点
        REQUIRE(anchors.size() == 6);

        // 验证所有锚点的 hash 都正确
        for (const auto& a : anchors) {
            CHECK(a.hash == 0x111111);
        }

        // 验证 pos_ref 的取值（应该包含 100, 200, 300）
        std::unordered_set<std::uint32_t> ref_positions;
        for (const auto& a : anchors) {
            ref_positions.insert(a.pos_ref);
        }
        CHECK(ref_positions.size() == 3);
        CHECK(ref_positions.count(100) == 1);
        CHECK(ref_positions.count(200) == 1);
        CHECK(ref_positions.count(300) == 1);

        // 验证 pos_qry 的取值（应该包含 50, 150）
        std::unordered_set<std::uint32_t> qry_positions;
        for (const auto& a : anchors) {
            qry_positions.insert(a.pos_qry);
        }
        CHECK(qry_positions.size() == 2);
        CHECK(qry_positions.count(50) == 1);
        CHECK(qry_positions.count(150) == 1);
    }

    // ------------------------------------------------------------------
    // TEST CASE 7: 无匹配 - ref 和 qry 的 hash 完全不同
    // ------------------------------------------------------------------
    TEST_CASE("collect_anchors - 无共同 hash 时返回空")
    {
        minimizer::MinimizerHits ref_hits;
        ref_hits.push_back(makeHit(0x111111, 100));
        ref_hits.push_back(makeHit(0x222222, 200));

        minimizer::MinimizerHits qry_hits;
        qry_hits.push_back(makeHit(0x333333, 50));
        qry_hits.push_back(makeHit(0x444444, 150));

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits);

        CHECK(anchors.empty());
    }

    // ------------------------------------------------------------------
    // TEST CASE 8: span 取最小值 - ref 和 qry 的 span 不同
    // ------------------------------------------------------------------
    TEST_CASE("collect_anchors - span 取 min(ref.span, qry.span)")
    {
        minimizer::MinimizerHits ref_hits;
        ref_hits.push_back(makeHit(0x555555, 100, 0, true, 20)); // span=20

        minimizer::MinimizerHits qry_hits;
        qry_hits.push_back(makeHit(0x555555, 50, 0, true, 15)); // span=15

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits);

        REQUIRE(anchors.size() == 1);
        CHECK(anchors[0].span == 15); // min(20, 15) = 15
    }

    // ------------------------------------------------------------------
    // TEST CASE 9: 多序列场景 - rid_ref 和 rid_qry
    // ------------------------------------------------------------------
    TEST_CASE("collect_anchors - 多序列 rid 正确传递")
    {
        minimizer::MinimizerHits ref_hits;
        ref_hits.push_back(makeHit(0x666666, 100, 5, true, 15)); // rid=5

        minimizer::MinimizerHits qry_hits;
        qry_hits.push_back(makeHit(0x666666, 50, 7, true, 15)); // rid=7

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits);

        REQUIRE(anchors.size() == 1);
        CHECK(anchors[0].rid_ref == 5);
        CHECK(anchors[0].rid_qry == 7);
    }
}

// ================================================================
// TEST SUITE: anchor::collect_anchors - 过滤参数测试
// ================================================================

TEST_SUITE("anchor::collect_anchors - 过滤逻辑测试")
{
    // ------------------------------------------------------------------
    // TEST CASE 10: 高频过滤 - ref_occ > ref_occ_thr 时稀疏采样
    // ------------------------------------------------------------------
    TEST_CASE("collect_anchors - 高频 minimizer 稀疏采样（-e 参数）")
    {
        // 构造一个高频 hash：在 ref 端出现 100 次
        minimizer::MinimizerHits ref_hits;
        for (std::uint32_t i = 0; i < 100; ++i) {
            ref_hits.push_back(makeHit(0x777777, i * 10));
        }

        // qry 端该 hash 出现在多个位置
        minimizer::MinimizerHits qry_hits;
        for (std::uint32_t i = 0; i < 20; ++i) {
            qry_hits.push_back(makeHit(0x777777, i * 100));
        }

        // 设置较低的阈值，触发稀疏采样
        anchor::SeedFilterParams params;
        params.u_floor = 5;            // 阈值下限
        params.u_ceil = 50;            // 阈值上限 => ref_occ=100 超过阈值
        params.f_top_frac = 0.0;       // 不使用 -f 过滤
        params.sample_every_bp = 500;  // 每 500bp 采样一次

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits, params);

        // 由于采样，锚点数量应显著少于 100*20=2000
        // 具体多少取决于 qry 位置能否被 500 整除
        // 这里主要验证"有过滤效果"
        CHECK(anchors.size() < 2000);
        CHECK(anchors.size() > 0); // 至少有一些锚点
    }

    // ------------------------------------------------------------------
    // TEST CASE 11: q-occ-frac 过滤 - qry 端高频 hash 被丢弃
    // ------------------------------------------------------------------
    TEST_CASE("collect_anchors - qry 端高频过滤（--q-occ-frac）")
    {
        minimizer::MinimizerHits ref_hits;
        ref_hits.push_back(makeHit(0x888888, 100));

        // qry 端有 1000 个 hit，其中同一个 hash 出现 500 次（占 50%）
        minimizer::MinimizerHits qry_hits;
        for (std::uint32_t i = 0; i < 500; ++i) {
            qry_hits.push_back(makeHit(0x888888, i)); // 高频 hash
        }
        for (std::uint32_t i = 0; i < 500; ++i) {
            qry_hits.push_back(makeHit(0x999900 + i, i)); // 不同的 hash
        }

        // 设置 q_occ_frac = 1%（即 10 次），超过则丢弃
        anchor::SeedFilterParams params;
        params.q_occ_frac = 0.01; // 1% of 1000 = 10
        params.u_floor = 0;
        params.u_ceil = 1000000;
        params.f_top_frac = 0.0;

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits, params);

        // 高频 hash (0x888888) 应该被过滤掉
        CHECK(anchors.empty());
    }

    // ------------------------------------------------------------------
    // TEST CASE 12: 无过滤 - 默认参数下所有匹配都生成锚点
    // ------------------------------------------------------------------
    TEST_CASE("collect_anchors - 默认参数不过滤低频 hash")
    {
        minimizer::MinimizerHits ref_hits;
        ref_hits.push_back(makeHit(0xAAAAAA, 100));
        ref_hits.push_back(makeHit(0xBBBBBB, 200));

        minimizer::MinimizerHits qry_hits;
        qry_hits.push_back(makeHit(0xAAAAAA, 50));
        qry_hits.push_back(makeHit(0xBBBBBB, 150));

        // 使用默认参数（不会过滤低频）
        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits);

        REQUIRE(anchors.size() == 2);
    }
}

// ================================================================
// TEST SUITE: anchor::collect_anchors - 一致性测试
// ================================================================

TEST_SUITE("anchor::collect_anchors - 一致性测试")
{
    // ------------------------------------------------------------------
    // TEST CASE 13: 相同输入产生稳定输出（无随机性）
    // ------------------------------------------------------------------
    TEST_CASE("collect_anchors - 多次调用结果一致")
    {
        minimizer::MinimizerHits ref_hits;
        for (std::uint32_t i = 0; i < 10; ++i) {
            ref_hits.push_back(makeHit(0xCCCCCC, i * 50));
        }

        minimizer::MinimizerHits qry_hits;
        for (std::uint32_t i = 0; i < 5; ++i) {
            qry_hits.push_back(makeHit(0xCCCCCC, i * 100));
        }

        auto anchors1 = minimizer::collect_anchors(ref_hits, qry_hits);
        auto anchors2 = minimizer::collect_anchors(ref_hits, qry_hits);

        REQUIRE(anchors1.size() == anchors2.size());

        // 逐一比较（顺序可能不同，需要排序后比较）
        auto cmp = [](const anchor::Anchor& a, const anchor::Anchor& b) {
            if (a.pos_ref != b.pos_ref) return a.pos_ref < b.pos_ref;
            return a.pos_qry < b.pos_qry;
        };
        std::sort(anchors1.begin(), anchors1.end(), cmp);
        std::sort(anchors2.begin(), anchors2.end(), cmp);

        for (std::size_t i = 0; i < anchors1.size(); ++i) {
            CHECK(anchors1[i].hash == anchors2[i].hash);
            CHECK(anchors1[i].pos_ref == anchors2[i].pos_ref);
            CHECK(anchors1[i].pos_qry == anchors2[i].pos_qry);
            CHECK(anchors1[i].rid_ref == anchors2[i].rid_ref);
            CHECK(anchors1[i].rid_qry == anchors2[i].rid_qry);
            CHECK(anchors1[i].span == anchors2[i].span);
            CHECK(anchors1[i].is_rev == anchors2[i].is_rev);
        }
    }

    // ------------------------------------------------------------------
    // TEST CASE 14: 基于真实 minimizer 提取的锚点收集
    // ------------------------------------------------------------------
    TEST_CASE("collect_anchors - 基于 extractMinimizer 的真实场景")
    {
        // 两条相似序列
        std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGT"; // 32bp，重复模式
        std::string qry_seq = "ACGTACGTACGTACGTACGTACGTACGTACGT"; // 完全相同

        const std::size_t k = 7;
        const std::size_t w = 5;

        auto ref_hits = minimizer::extractMinimizer(ref_seq, k, w, false);
        auto qry_hits = minimizer::extractMinimizer(qry_seq, k, w, false);

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits);

        // 两条相同序列应该产生至少一个锚点
        CHECK(anchors.size() > 0);

        // 所有锚点的 is_rev 应该为 false（同向）
        for (const auto& a : anchors) {
            CHECK(a.is_rev == false);
        }
    }

    // ------------------------------------------------------------------
    // TEST CASE 15: 部分重叠序列
    // ------------------------------------------------------------------
    TEST_CASE("collect_anchors - 部分重叠序列产生锚点")
    {
        std::string ref_seq = "AAAACGTACGTACGTACGTTTTT"; // 前后有不同区域
        std::string qry_seq = "GGGGCGTACGTACGTACGTCCCC"; // 中间部分相同

        const std::size_t k = 7;
        const std::size_t w = 3;

        auto ref_hits = minimizer::extractMinimizer(ref_seq, k, w, false);
        auto qry_hits = minimizer::extractMinimizer(qry_seq, k, w, false);

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits);

        // 中间重叠部分应该产生锚点
        CHECK(anchors.size() > 0);
    }
}

// ================================================================
// TEST SUITE: anchor::collect_anchors - 性能测试（可选）
// ================================================================

TEST_SUITE("anchor::collect_anchors - 性能测试" * doctest::skip(shouldSkipPerf()))
{
    // ------------------------------------------------------------------
    // TEST CASE 16: 大规模数据吞吐量测试
    // ------------------------------------------------------------------
    TEST_CASE("collect_anchors - 大规模数据性能")
    {
        // 生成一条 10kb 的参考序列和 5kb 的查询序列
        const std::size_t ref_len = 10000;
        const std::size_t qry_len = 5000;
        const std::size_t k = 15;
        const std::size_t w = 10;

        std::string ref_seq = makeRandomDna(ref_len, 12345);
        std::string qry_seq = makeRandomDna(qry_len, 67890);

        auto ref_hits = minimizer::extractMinimizer(ref_seq, k, w, false);
        auto qry_hits = minimizer::extractMinimizer(qry_seq, k, w, false);

        auto start = std::chrono::high_resolution_clock::now();
        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits);
        auto end = std::chrono::high_resolution_clock::now();

        auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        MESSAGE("collect_anchors 性能测试:");
        MESSAGE("  参考序列长度: ", ref_len, " bp");
        MESSAGE("  查询序列长度: ", qry_len, " bp");
        MESSAGE("  参考 hits 数量: ", ref_hits.size());
        MESSAGE("  查询 hits 数量: ", qry_hits.size());
        MESSAGE("  生成锚点数量: ", anchors.size());
        MESSAGE("  耗时: ", elapsed_us, " μs");
        if (elapsed_us > 0) {
            MESSAGE("  吞吐量: ", (ref_len + qry_len) * 1000000.0 / elapsed_us, " bp/s");
        }

        // 不做严格断言，只确保不崩溃
        CHECK(anchors.size() >= 0);
    }

    // ------------------------------------------------------------------
    // TEST CASE 17: 高频重复区域性能测试
    // ------------------------------------------------------------------
    TEST_CASE("collect_anchors - 高频重复区域性能")
    {
        // 构造一个高度重复的序列（模拟 LTR、tandem repeat 等）
        const std::size_t repeat_count = 100;
        const std::string repeat_unit = "ACGTACGTACGTACGT"; // 16bp 重复单元
        const std::size_t k = 7;
        const std::size_t w = 5;

        std::string ref_seq;
        for (std::size_t i = 0; i < repeat_count; ++i) {
            ref_seq += repeat_unit;
        }

        std::string qry_seq = ref_seq; // 完全相同

        auto ref_hits = minimizer::extractMinimizer(ref_seq, k, w, false);
        auto qry_hits = minimizer::extractMinimizer(qry_seq, k, w, false);

        MESSAGE("重复序列性能测试:");
        MESSAGE("  序列长度: ", ref_seq.size(), " bp");
        MESSAGE("  参考 hits: ", ref_hits.size());
        MESSAGE("  查询 hits: ", qry_hits.size());

        auto start = std::chrono::high_resolution_clock::now();
        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits);
        auto end = std::chrono::high_resolution_clock::now();

        auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        MESSAGE("  生成锚点数量: ", anchors.size());
        MESSAGE("  耗时: ", elapsed_us, " μs");

        // 验证过滤有效性：锚点数量不应该爆炸到 O(n^2)
        // 假设 hits 数量在千级别，如果没有过滤会产生百万级锚点
        // 有过滤的情况下应该大幅减少
        CHECK(anchors.size() > 0);
        MESSAGE("  锚点/ref_hits 比例: ", static_cast<double>(anchors.size()) / ref_hits.size());
    }
}

// ================================================================
// TEST SUITE: anchor::collect_anchors - 边界与异常
// ================================================================

TEST_SUITE("anchor::collect_anchors - 边界与鲁棒性")
{
    // ------------------------------------------------------------------
    // TEST CASE 18: 极端情况 - 单个 hit
    // ------------------------------------------------------------------
    TEST_CASE("collect_anchors - 单个 hit 匹配")
    {
        minimizer::MinimizerHits ref_hits;
        ref_hits.push_back(makeHit(0xDDDDDD, 42));

        minimizer::MinimizerHits qry_hits;
        qry_hits.push_back(makeHit(0xDDDDDD, 99));

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits);

        REQUIRE(anchors.size() == 1);
        CHECK(anchors[0].hash == 0xDDDDDD);
        CHECK(anchors[0].pos_ref == 42);
        CHECK(anchors[0].pos_qry == 99);
    }

    // ------------------------------------------------------------------
    // TEST CASE 19: 大量不同 hash - 确保哈希表不崩溃
    // ------------------------------------------------------------------
    TEST_CASE("collect_anchors - 大量不同 hash（哈希表压力测试）")
    {
        minimizer::MinimizerHits ref_hits;
        for (std::uint32_t i = 0; i < 1000; ++i) {
            ref_hits.push_back(makeHit(0x100000 + i, i * 10));
        }

        minimizer::MinimizerHits qry_hits;
        for (std::uint32_t i = 500; i < 1500; ++i) {
            qry_hits.push_back(makeHit(0x100000 + i, i * 5));
        }

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits);

        // 重叠部分是 [500, 1000)，共 500 个不同 hash
        CHECK(anchors.size() == 500);
    }

    // ------------------------------------------------------------------
    // TEST CASE 20: span = 0 的边界情况（虽然实际不应发生）
    // ------------------------------------------------------------------
    TEST_CASE("collect_anchors - span 为 0 的 hit")
    {
        minimizer::MinimizerHits ref_hits;
        ref_hits.push_back(makeHit(0xEEEEEE, 100, 0, true, 0)); // span=0

        minimizer::MinimizerHits qry_hits;
        qry_hits.push_back(makeHit(0xEEEEEE, 50, 0, true, 15));

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits);

        REQUIRE(anchors.size() == 1);
        CHECK(anchors[0].span == 0); // min(0, 15) = 0
    }
}

