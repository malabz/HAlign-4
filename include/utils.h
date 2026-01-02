//
// Created by 30451 on 2026/1/2.
//

#ifndef HALIGN4_UTILS_H
#define HALIGN4_UTILS_H

namespace file_io
{
namespace fs = std::filesystem;

inline std::string formatFsError(std::string_view msg,
                                 const fs::path& p,
                                 const std::error_code& ec) {
    std::string s;
    s.reserve(msg.size() + p.string().size() + 64);
    s.append(msg);
    s.append(": ");
    s.append(p.string());
    if (ec) {
        s.append(" (");
        s.append(ec.message());
        s.append(")");
    }
    return s;
}

inline void requireExists(const fs::path& p, std::string_view what) {
    std::error_code ec;
    const bool ok = fs::exists(p, ec); // non-throwing overload clears/sets ec accordingly :contentReference[oaicite:1]{index=1}
    if (ec || !ok) {
        throw std::runtime_error(formatFsError(std::string(what) + " does not exist", p, ec));
    }
}

inline void requireRegularFile(const fs::path& p, std::string_view what) {
    std::error_code ec;
    const bool ok = fs::is_regular_file(p, ec); // non-throwing overload :contentReference[oaicite:2]{index=2}
    if (ec || !ok) {
        throw std::runtime_error(formatFsError(std::string(what) + " is not a regular file", p, ec));
    }
}

inline void requireDirectory(const fs::path& p, std::string_view what) {
    std::error_code ec;
    const bool ok = fs::is_directory(p, ec); // non-throwing overload :contentReference[oaicite:3]{index=3}
    if (ec || !ok) {
        throw std::runtime_error(formatFsError(std::string(what) + " is not a directory", p, ec));
    }
}

inline void ensureDirectoryExists(const fs::path& p, std::string_view what = "directory") {
    std::error_code ec;
    if (fs::exists(p, ec)) {
        if (ec) {
            throw std::runtime_error(formatFsError(std::string(what) + " status check failed", p, ec));
        }
        requireDirectory(p, what);
        return;
    }
    if (ec) {
        throw std::runtime_error(formatFsError(std::string(what) + " status check failed", p, ec));
    }

    // create_directories has error_code overload (non-throwing) :contentReference[oaicite:4]{index=4}
    fs::create_directories(p, ec);
    if (ec) {
        throw std::runtime_error(formatFsError(std::string("failed to create ") + std::string(what), p, ec));
    }
}

inline bool isEmpty(const fs::path& p) {
    std::error_code ec;
    const bool empty = fs::is_empty(p, ec); // checks empty file or directory :contentReference[oaicite:5]{index=5}
    if (ec) {
        throw std::runtime_error(formatFsError("failed to check emptiness", p, ec));
    }
    return empty;
}

// 核心：准备工作目录（不存在则创建；存在则要求为目录；可选要求为空）
inline void prepareEmptydir(const fs::path& workdir, bool must_be_empty) {
    if (workdir.empty()) {
        throw std::runtime_error("workdir is empty");
    }

    ensureDirectoryExists(workdir, "workdir");

    if (must_be_empty && !isEmpty(workdir)) {
        throw std::runtime_error("workdir must be empty: " + workdir.string());
    }
}

// 可选：确保输出文件的父目录存在
inline void ensureParentDirExists(const fs::path& out_file) {
    if (out_file.empty()) return;
    auto parent = out_file.parent_path();
    if (parent.empty()) return;
    ensureDirectoryExists(parent, "output parent dir");
}

// 可选：危险操作（如你未来需要清空 workdir）
// remove_all 的语义与返回值见标准库文档 :contentReference[oaicite:6]{index=6}
inline void removeAll(const fs::path& p) {
    std::error_code ec;
    fs::remove_all(p, ec);
    if (ec) {
        throw std::runtime_error(formatFsError("remove_all failed", p, ec));
    }
}
}


#endif //HALIGN4_UTILS_H