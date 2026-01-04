#include "utils.h"
#include <cstdio>
#include <cstdlib> // for malloc/free

#if __has_include(<zlib.h>)
    #include <zlib.h>
    #define HALIGN4_HAVE_ZLIB 1
#else
    #define HALIGN4_HAVE_ZLIB 0
#endif

#include "kseq.h"

#if HALIGN4_HAVE_ZLIB
    KSEQ_INIT(gzFile, gzread)
#else
    static int fileRead(std::FILE* fp, void* buf, int len)
    {
        const std::size_t n = std::fread(buf, 1, static_cast<std::size_t>(len), fp);
        return static_cast<int>(n);
    }
    KSEQ_INIT(std::FILE*, fileRead)
#endif

namespace seq_io
{
    struct KseqReader::Impl
    {

#if HALIGN4_HAVE_ZLIB
        gzFile fp{nullptr};
#else
        std::FILE* fp{nullptr};
#endif
        kseq_t* seq{nullptr};
        FilePath file_path; // store the input path for diagnostics and error messages

        // 新增：用于加大 stdio / zlib 缓冲的持久化缓冲区（非 zlib 情况）
        // 这样可以减少 fread/gzread 的系统调用次数，提高大文件读取吞吐量
        char* io_buf{nullptr};
        std::size_t io_buf_size{1 << 20}; // 默认 1 MiB，可根据需要调整
    };

    static std::runtime_error makeIoError(const std::string& msg, const FilePath& p)
    {
        return std::runtime_error(msg + ": " + p.string());
    }

    static void assignKstring(std::string& dst, const kstring_t& ks)
    {
        if (ks.s && ks.l > 0) dst.assign(ks.s, ks.l);
        else dst.clear();
    }

    KseqReader::KseqReader(const FilePath& file_path)
        : impl_(std::make_unique<Impl>())
    {
        impl_->file_path = file_path;

#if HALIGN4_HAVE_ZLIB
        // gzopen/gzread 是 kseq 常用组合。
        impl_->fp = gzopen(file_path.string().c_str(), "rb");
        if (!impl_->fp) {
            throw makeIoError("failed to open input", file_path);
        }

        // 提升 zlib 的内部缓冲区大小以加速大文件读取（单位：字节，接口为 int）
        // 如果 io_buf_size 超过 int 上限，gzbuffer 的行为未定义；这里选用 1MiB，安全且有效。
        gzbuffer(impl_->fp, static_cast<int>(impl_->io_buf_size));

        impl_->seq = kseq_init(impl_->fp);
#else
        impl_->fp = std::fopen(file_path.string().c_str(), "rb");
        if (!impl_->fp) {
            throw makeIoError("failed to open input", file_path);
        }

        // 为 stdio 分配一个较大的缓冲区并设置为全缓冲，减少 fread 的系统调用次数
        impl_->io_buf = static_cast<char*>(std::malloc(impl_->io_buf_size));
        if (impl_->io_buf) {
            if (setvbuf(impl_->fp, impl_->io_buf, _IOFBF, static_cast<size_t>(impl_->io_buf_size)) != 0) {
                // 如果设置失败，释放并继续（仍能工作，只是没有额外加速）
                std::free(impl_->io_buf);
                impl_->io_buf = nullptr;
            }
        }

        impl_->seq = kseq_init(impl_->fp);
#endif

        if (!impl_->seq) {
            throw makeIoError("failed to init kseq", file_path);
        }
    }

    KseqReader::~KseqReader()
    {
        if (!impl_) return;

        if (impl_->seq) {
            kseq_destroy(impl_->seq);
            impl_->seq = nullptr;
        }

#if HALIGN4_HAVE_ZLIB
        if (impl_->fp) {
            gzclose(impl_->fp);
            impl_->fp = nullptr;
        }
#else
        if (impl_->fp) {
            std::fclose(impl_->fp);
            impl_->fp = nullptr;
        }
        // 释放为 stdio 分配的缓冲区（必须在 fclose 之后）
        if (impl_->io_buf) {
            std::free(impl_->io_buf);
            impl_->io_buf = nullptr;
        }
#endif
    }

    KseqReader::KseqReader(KseqReader&& other) noexcept = default;
    KseqReader& KseqReader::operator=(KseqReader&& other) noexcept = default;

    bool KseqReader::next(SeqRecord& rec)
    {
        if (!impl_ || !impl_->seq) {
            throw std::runtime_error("KseqReader is not initialized");
        }

        const int ret = kseq_read(impl_->seq);

        if (ret >= 0) {
            assignKstring(rec.id,   impl_->seq->name);
            assignKstring(rec.desc, impl_->seq->comment);
            assignKstring(rec.seq,  impl_->seq->seq);
            return true;
        }

        if (ret == -1) {
            return false; // EOF
        }

        // 其它负值视为解析错误（如 FASTQ 质量串截断等场景）
        throw std::runtime_error("kseq_read() failed with code " + std::to_string(ret) +
                                 " for file: " + impl_->file_path.string());
    }

    // ------------------------- FastaWriter -------------------------

    FastaWriter::FastaWriter(const FilePath& file_path, std::size_t line_width)
        : out_(file_path, std::ios::binary), line_width_(line_width)
    {
        if (!out_) {
            throw makeIoError("failed to open output", file_path);
        }
        if (line_width_ == 0) line_width_ = 80;
    }

    void FastaWriter::writeWrapped(std::ofstream& out, std::string_view s, std::size_t width)
    {
        for (std::size_t i = 0; i < s.size(); i += width) {
            const std::size_t n = (i + width <= s.size()) ? width : (s.size() - i);
            out.write(s.data() + static_cast<std::streamoff>(i), static_cast<std::streamsize>(n));
            out.put('\n');
        }
    }

    // cpp
    void FastaWriter::write(const SeqRecord& rec)
    {
        if (!out_) throw std::runtime_error("FastaWriter output stream is not ready");

        // 构造并一次性写入 header
        std::string header;
        header.reserve(1 + rec.id.size() + (rec.desc.empty() ? 0 : 1 + rec.desc.size()) + 1);
        header.push_back('>');
        header.append(rec.id);
        if (!rec.desc.empty()) {
            header.push_back(' ');
            header.append(rec.desc);
        }
        header.push_back('\n');
        out_.write(header.data(), static_cast<std::streamsize>(header.size()));

        // 构造带换行的序列缓冲区并一次性写出，减少多次 write/put 调用
        const std::size_t L = rec.seq.size();
        const std::size_t width = (line_width_ == 0) ? 80 : line_width_;

        if (L == 0) {
            out_.put('\n');
            return;
        }

        std::string seqbuf;
        seqbuf.reserve(L + (L / width) + 1);
        for (std::size_t i = 0; i < L; i += width) {
            const std::size_t n = (i + width <= L) ? width : (L - i);
            seqbuf.append(rec.seq.data() + static_cast<std::size_t>(i), n);
            seqbuf.push_back('\n');
        }

        out_.write(seqbuf.data(), static_cast<std::streamsize>(seqbuf.size()));
    }

    void FastaWriter::flush()
    {
        out_.flush();
    }

} // namespace seq_io
