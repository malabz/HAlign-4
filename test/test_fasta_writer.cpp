#include <doctest/doctest.h>

#include <filesystem>
#include <fstream>
#include <string>
#include <chrono>
#include <cstdlib>

#include "utils.h"

namespace fs = std::filesystem;

// perf gating：只有设置 HALIGN4_RUN_PERF=1 才会执行性能用例
static bool perfEnabled() {
    const char* v = std::getenv("HALIGN4_RUN_PERF");
    return (v != nullptr) && (*v != '\0') && (std::string(v) != "0");
}
static bool shouldSkipPerf() { return !perfEnabled(); }

static fs::path makeTempDir(std::string_view name) {
    fs::path base = fs::current_path();
    fs::path dir = base / std::string(name);
    std::error_code ec;
    fs::remove_all(dir, ec);
    fs::create_directories(dir, ec);
    REQUIRE_MESSAGE(!ec, "cannot create temp dir: " << dir.string() << " (" << ec.message() << ")");
    return dir;
}

static std::string slurp(const fs::path& p) {
    std::ifstream ifs(p, std::ios::binary);
    REQUIRE_MESSAGE(ifs.good(), "cannot open file: " << p.string());
    std::string s;
    ifs.seekg(0, std::ios::end);
    s.resize(static_cast<std::size_t>(ifs.tellg()));
    ifs.seekg(0, std::ios::beg);
    if (!s.empty()) ifs.read(s.data(), static_cast<std::streamsize>(s.size()));
    return s;
}

static double toMiB(double bytes) { return bytes / (1024.0 * 1024.0); }

TEST_SUITE("write_fasta")
{
    TEST_CASE("FastaWriter - buffered flush on flush() and destructor")
    {
        auto dir = makeTempDir("halign4_tests_fasta_writer");
        fs::path out = dir / "out.fasta";

        // 使用一个很大的阈值，确保 write() 过程中不会自动触发 flushBuffer_。
        {
            seq_io::FastaWriter w(out, /*line_width=*/4, /*buffer_threshold_bytes=*/1ULL << 30);

            seq_io::SeqRecord r1{"id1", "", "ACGTACGT", ""};
            seq_io::SeqRecord r2{"id2", "desc", "TT", ""};

            w.write(r1);
            w.write(r2);

            // 未显式 flush() 时，数据应该仍在内存缓冲区里；文件可能还没落盘。
            // 注意：不同标准库/平台对文件大小可见性有差异，所以这里用“==0 或非常小”做弱断言。
            std::error_code ec;
            auto sz = fs::exists(out, ec) ? fs::file_size(out, ec) : 0ULL;
            CHECK_MESSAGE(!ec, "file_size failed: " << ec.message());
            CHECK_MESSAGE(sz == 0ULL, "expected buffered writer not to write on disk before flush (size=" << sz << ")");

            // flush 后应落盘
            w.flush();
        }

        const std::string got = slurp(out);
        const std::string expected =
            ">id1\n"
            "ACGT\n"
            "ACGT\n"
            ">id2 desc\n"
            "TT\n";

        CHECK(got == expected);

        std::error_code ec;
        fs::remove_all(dir, ec);
    }

    TEST_CASE("FastaWriter - disable buffer (threshold=0) writes immediately")
    {
        auto dir = makeTempDir("halign4_tests_fasta_writer_nobuf");
        fs::path out = dir / "out.fasta";

        seq_io::FastaWriter w(out, /*line_width=*/80, /*buffer_threshold_bytes=*/0);
        seq_io::SeqRecord r{"id", "", "AAAA", ""};
        w.write(r);

        // 关闭额外缓冲时，write() 会直接写入到 std::ofstream 的用户态缓冲区，
        // 但文件系统层面的“文件大小可见性”可能仍依赖 flush()/close()。
        // 因此这里先 flush，再检查大小与内容，保证跨平台稳定。
        w.flush();

        std::error_code ec;
        const auto sz = fs::file_size(out, ec);
        CHECK_MESSAGE(!ec, "file_size failed: " << ec.message());
        CHECK(sz > 0);

        const std::string got = slurp(out);
        const std::string expected = ">id\nAAAA\n";
        CHECK(got == expected);

        fs::remove_all(dir, ec);
    }

    TEST_CASE("perf: write time n=10000 len=30000 (FastaWriter)")
    {
        if (shouldSkipPerf()) {
            doctest::skip(true);
        }

        constexpr std::size_t N   = 10000;
        constexpr std::size_t LEN = 30000;

        auto dir = makeTempDir("halign4_tests_write_fasta_perf");
        fs::path out_buf   = dir / "out_buffered.fasta";
        fs::path out_nobuf = dir / "out_nobuf.fasta";

        // 准备测试数据：复用同一条序列内容，避免把随机数生成也算进写出时间里
        seq_io::SeqRecord rec;
        rec.id = "s";
        rec.desc.clear();
        rec.seq.assign(LEN, 'A');

        // 估算写出字节数（用于吞吐率）：每条记录大约为
        // - header: ">s\n" -> 3
        // - sequence: LEN + "\n" -> LEN+1
        const double approx_bytes = static_cast<double>(N) * static_cast<double>(3 + LEN + 1);

        // 1) 默认 writer（内部额外缓冲阈值约 8MiB）
        {
            seq_io::FastaWriter w(out_buf);
            auto t0 = std::chrono::steady_clock::now();
            for (std::size_t i = 0; i < N; ++i) {
                rec.id = "s" + std::to_string(i);
                w.write(rec);
            }
            w.flush();
            auto t1 = std::chrono::steady_clock::now();
            const double sec = std::chrono::duration<double>(t1 - t0).count();
            MESSAGE("write_buffered: n=" << N << " len=" << LEN
                    << " time_s=" << sec
                    << " approx_throughput_MiBps=" << (sec > 0.0 ? (toMiB(approx_bytes) / sec) : 0.0));
        }

        // 2) 关闭额外缓冲（阈值=0）
        {
            seq_io::FastaWriter w(out_nobuf, /*line_width=*/80, /*buffer_threshold_bytes=*/0);
            auto t0 = std::chrono::steady_clock::now();
            for (std::size_t i = 0; i < N; ++i) {
                rec.id = "s" + std::to_string(i);
                w.write(rec);
            }
            w.flush();
            auto t1 = std::chrono::steady_clock::now();
            const double sec = std::chrono::duration<double>(t1 - t0).count();
            MESSAGE("write_nobuf: n=" << N << " len=" << LEN
                    << " time_s=" << sec
                    << " approx_throughput_MiBps=" << (sec > 0.0 ? (toMiB(approx_bytes) / sec) : 0.0));
        }

        // 基本正确性 smoke：输出文件应非空（不做逐条解析，避免 perf 用例耗时过长）
        std::error_code ec;
        CHECK(fs::file_size(out_buf, ec) > 0);
        CHECK(fs::file_size(out_nobuf, ec) > 0);

        fs::remove_all(dir, ec);
    }
}
