#include <doctest/doctest.h>

#include <filesystem>
#include <fstream>
#include <string>

#include "consensus.h"

namespace fs = std::filesystem;

static void writeTextFile(const fs::path& p, const std::string& content) {
    std::ofstream ofs(p, std::ios::binary);
    REQUIRE_MESSAGE(ofs.good(), "cannot write file: " << p.string());
    ofs << content;
}

static std::string readSingleFastaSequence(const fs::path& fasta) {
    std::ifstream ifs(fasta, std::ios::binary);
    REQUIRE_MESSAGE(ifs.good(), "cannot read fasta: " << fasta.string());

    std::string line, seq;
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') continue;
        seq += line;
    }
    return seq;
}

static fs::path makeTempDirSimple() {
    fs::path dir = fs::temp_directory_path() / "halign4_tests_simple";
    std::error_code ec;
    fs::remove_all(dir, ec);
    fs::create_directories(dir, ec);
    REQUIRE_MESSAGE(!ec, "cannot create temp dir: " << dir.string() << " (" << ec.message() << ")");
    return dir;
}

// ---------- 测试用例 ----------
// 测试说明：
// 1) 构造一个包含 3 条对齐序列的小 FASTA（长度为 5），其中某些位置有 gap('-')；
// 2) 调用 consensus::generateConsensusSequence（批次大小和线程数可控）生成共识；
// 3) 验证返回的共识字符串与预期一致，验证输出 FASTA 中的序列也一致；
// 4) 验证 counts.json 存在且非空（简单检查，不逐字段比对）

TEST_SUITE("generate_consensus") {

TEST_CASE("generateConsensusSequence - correctness (gap majority is ignored)") {
    auto dir = makeTempDirSimple();
    fs::path in_fa  = dir / "aligned.fasta";
    fs::path out_fa = dir / "consensus.fasta";
    fs::path out_js = dir / "counts.json";

    // s1 A C G T -
    // s2 A C - T -
    // s3 A C G T -
    //
    // 你的 pickConsensusChar 只在 A/C/G/T/U 中挑多数，不会输出 '-'.
    // 第 5 列全是 '-' => A/C/G/T/U 都为 0 => 返回 'A'.
    // 所以期望共识是 "ACGTA"，不是 "ACGT-'.
    const std::string aligned =
        ">s1\nACGT-\n"
        ">s2\nAC-T-\n"
        ">s3\nACGT-\n";

    writeTextFile(in_fa, aligned);

    std::string cons = consensus::generateConsensusSequence(
        in_fa, out_fa, out_js,
        0, 4, 2
    );

    CHECK(cons == "ACGTA");
    CHECK(readSingleFastaSequence(out_fa) == "ACGTA");

    // JSON 先只检查：文件存在且非空（避免 cereal 字段名不一致导致读失败）
    std::error_code ec;
    CHECK(fs::exists(out_js, ec));
    if (!ec && fs::exists(out_js, ec)) {
        CHECK(fs::file_size(out_js, ec) > 0);
    }
}

TEST_CASE("generateConsensusSequence - tie breaks to A (A > C > G > T > U)") {
    auto dir = makeTempDirSimple();
    fs::path in_fa  = dir / "aligned_tie.fasta";
    fs::path out_fa = dir / "consensus_tie.fasta";
    fs::path out_js = dir / "counts_tie.json";

    // 两条序列在该位点 A=1, C=1 => tie，按你的逻辑不更新 => 返回 'A'
    const std::string aligned =
        ">s1\nA\n"
        ">s2\nC\n";

    writeTextFile(in_fa, aligned);

    std::string cons = consensus::generateConsensusSequence(
        in_fa, out_fa, out_js,
        0, 2, 2
    );

    CHECK(cons == "A");
    CHECK(readSingleFastaSequence(out_fa) == "A");
}

TEST_CASE("generateConsensusSequence - all gaps => consensus becomes A") {
    auto dir = makeTempDirSimple();
    fs::path in_fa  = dir / "aligned_allgap.fasta";
    fs::path out_fa = dir / "consensus_allgap.fasta";
    fs::path out_js = dir / "counts_allgap.json";

    // 所有位点全是 '-'，A/C/G/T/U 都是 0 => 每一列都返回 'A'
    const std::string aligned =
        ">s1\n---\n"
        ">s2\n---\n"
        ">s3\n---\n";

    writeTextFile(in_fa, aligned);

    std::string cons = consensus::generateConsensusSequence(
        in_fa, out_fa, out_js,
        0, 4, 2
    );

    CHECK(cons == "AAA");
    CHECK(readSingleFastaSequence(out_fa) == "AAA");
}

// ========== 新增测试用例 ==========

TEST_CASE("generateConsensusSequence - single sequence returns same") {
    auto dir = makeTempDirSimple();
    fs::path in_fa  = dir / "single.fasta";
    fs::path out_fa = dir / "single_cons.fasta";
    fs::path out_js = dir / "single_counts.json";

    const std::string seq = ">s1\nACGTACGT\n";
    writeTextFile(in_fa, seq);

    std::string cons = consensus::generateConsensusSequence(
        in_fa, out_fa, out_js,
        0, 1, 1
    );

    CHECK(cons == "ACGTACGT");
    CHECK(readSingleFastaSequence(out_fa) == "ACGTACGT");
}

TEST_CASE("generateConsensusSequence - U (uracil) handling") {
    auto dir = makeTempDirSimple();
    fs::path in_fa  = dir / "u_vs_t.fasta";
    fs::path out_fa = dir / "u_vs_t_cons.fasta";
    fs::path out_js = dir / "u_vs_t_counts.json";

    // 第一位 U=2, T=1 => 共识 U
    const std::string aligned =
        ">s1\nUAAAA\n"
        ">s2\nTAAAA\n"
        ">s3\nUAAAA\n";

    writeTextFile(in_fa, aligned);

    std::string cons = consensus::generateConsensusSequence(
        in_fa, out_fa, out_js,
        0, 2, 2
    );

    CHECK(cons.size() == 5);
    CHECK(cons[0] == 'U');
}

TEST_CASE("generateConsensusSequence - seq_limit affects result") {
    auto dir = makeTempDirSimple();
    fs::path in_fa  = dir / "limit.fasta";
    fs::path out_fa = dir / "limit_cons.fasta";
    fs::path out_js = dir / "limit_counts.json";

    // 构造 5 条序列：前 2 条 A，后 3 条 C => 全体多数为 C，但若只处理前 2 条则为 A
    std::string content;
    content += ">s0\nA\n";
    content += ">s1\nA\n";
    content += ">s2\nC\n";
    content += ">s3\nC\n";
    content += ">s4\nC\n";
    writeTextFile(in_fa, content);

    // 不限制 => 5 条 => 共识 C
    std::string cons_all = consensus::generateConsensusSequence(
        in_fa, out_fa, out_js,
        0, 1, 2
    );
    CHECK(cons_all == "C");

    // 限制为 2 => 只考虑前两条 => 共识 A
    std::string cons_lim = consensus::generateConsensusSequence(
        in_fa, out_fa, out_js,
        2, 1, 2
    );
    CHECK(cons_lim == "A");
}

TEST_CASE("generateConsensusSequence - empty input throws") {
    auto dir = makeTempDirSimple();
    fs::path in_fa  = dir / "empty.fasta";
    fs::path out_fa = dir / "empty_cons.fasta";
    fs::path out_js = dir / "empty_counts.json";

    // 写一个空文件
    writeTextFile(in_fa, "");

    CHECK_THROWS_AS(
        consensus::generateConsensusSequence(in_fa, out_fa, out_js, 0, 1, 1),
        std::runtime_error
    );
}

// ========== 新增测试用例结束 ==========

} // end TEST_SUITE("generate_consensus")

