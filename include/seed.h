#ifndef HALIGN4_SEED_H
#define HALIGN4_SEED_H

#include <cstdint>
#include <type_traits>
#include <string_view>
#include <stdexcept>
#include <utility>
#include <vector>
#include <string>
#include <cstddef>
#include <algorithm>
#include "hash.h"

// ================================================================
// 抽象 seed 接口（不绑定具体 minimizer/syncmer/strobemer 实现）
//
// 目标：
// - 你后续可能支持多种 seed（minimizer 及各种变体，例如 syncmer/strobemer）
// - 同时希望在不同阶段选择不同存储：只存 hash 或存 hash+位置
// - 这里仅提供“最抽象”的接口与工具（traits/比较器），不提供任何具体实现。
// ================================================================
namespace seed
{
    // ------------------------------------------------------------
    // SeedKind：建议放在“容器/批次”层面共享（节省内存）
    // ------------------------------------------------------------
    enum class SeedKind : std::uint8_t
    {
        minimizer = 0,
        syncmer   = 1,
        strobemer = 2,
    };

    inline constexpr const char* seedKindToString(SeedKind k) noexcept
    {
        switch (k) {
        case SeedKind::minimizer: return "minimizer";
        case SeedKind::syncmer:   return "syncmer";
        case SeedKind::strobemer: return "strobemer";
        default:                  return "unknown";
        }
    }


    // ------------------------------------------------------------
    // CRTP 基类：hash+位置 seed（hit）
    // ------------------------------------------------------------
    // 派生类需要提供：
    //   hash_t hash() const noexcept;     // hash 视角
    //   uint32_t pos() const noexcept;      // 位置（0-based）
    //   uint32_t rid() const noexcept;      // 序列 id（多序列场景）
    //   bool strand() const noexcept;       // 方向
    //   uint32_t span() const noexcept;     // 覆盖范围（minimizer 一般等于 k；strobemer 可更大）
    template <typename Derived>
    struct SeedHitBase
    {
        constexpr hash_t hash() const noexcept { return static_cast<const Derived&>(*this).hash(); }
        constexpr std::uint32_t pos() const noexcept { return static_cast<const Derived&>(*this).pos(); }
        constexpr std::uint32_t rid() const noexcept { return static_cast<const Derived&>(*this).rid(); }
        constexpr bool strand() const noexcept { return static_cast<const Derived&>(*this).strand(); }
        constexpr std::uint32_t span() const noexcept { return static_cast<const Derived&>(*this).span(); }

        // 一个通用排序：先 hash，再 rid/pos/strand
        friend constexpr bool operator<(const SeedHitBase& a, const SeedHitBase& b) noexcept
        {
            if (a.hash() != b.hash()) return a.hash() < b.hash();
            if (a.rid() != b.rid()) return a.rid() < b.rid();
            if (a.pos() != b.pos()) return a.pos() < b.pos();
            return a.strand() < b.strand();
        }

        friend constexpr bool operator==(const SeedHitBase& a, const SeedHitBase& b) noexcept
        {
            return a.hash() == b.hash() && a.rid() == b.rid() && a.pos() == b.pos() && a.strand() == b.strand() && a.span() == b.span();
        }
    };

    // ------------------------------------------------------------
    // 统一访问接口（free-functions）
    // ------------------------------------------------------------

    // 默认 trait：没有位置信息（你仍然可以手动特化覆盖）
    template <typename T>
    struct has_position : std::false_type {};

    // 自动推导：如果 T 继承了 SeedHitBase<T>，则认为它“有位置”。
    // 这样后续你每增加一个 Hit 类型，只要继承 SeedHitBase 就自动生效，不需要再写特化。
    template <typename T>
    struct has_position_auto : std::is_base_of<SeedHitBase<T>, T> {};

    template <typename T>
    inline constexpr bool has_position_v = has_position<T>::value || has_position_auto<T>::value;

    // hash-only：要求类型提供 hash()
    template <typename SeedT>
    inline constexpr hash_t hash_value(const SeedT& s) noexcept
    {
        return s.hash();
    }

    // hit：要求类型提供 hash()/pos()/rid()/strand()/span()
    template <typename HitT>
    inline constexpr std::uint32_t get_pos(const HitT& h) noexcept { return h.pos(); }
    template <typename HitT>
    inline constexpr std::uint32_t get_rid(const HitT& h) noexcept { return h.rid(); }
    template <typename HitT>
    inline constexpr bool get_strand(const HitT& h) noexcept { return h.strand(); }
    template <typename HitT>
    inline constexpr std::uint32_t get_span(const HitT& h) noexcept { return h.span(); }

    // 只看 hash 的比較器（適用於 sort+unique / Jaccard / containment 等）
    struct HashOnlyLess
    {
        template <typename M>
        constexpr bool operator()(const M& a, const M& b) const noexcept
        {
            return hash_value(a) < hash_value(b);
        }
    };

    struct HashOnlyEqual
    {
        template <typename M>
        constexpr bool operator()(const M& a, const M& b) const noexcept
        {
            return hash_value(a) == hash_value(b);
        }
    };



    struct Anchor
    {
        hash_t hash{};              // 一般就是 seed hash（期望 ref/query 相同）
        std::uint32_t rid_ref{};    // ref 序列 id
        std::uint32_t pos_ref{};    // ref 上位置（0-based）
        std::uint32_t rid_qry{};    // query 序列 id（很多场景固定 0）
        std::uint32_t pos_qry{};    // query 上位置（0-based, forward 坐标系）
        std::uint32_t span{};       // 覆盖长度（可用 min(ref.span, qry.span)）
        bool is_rev{};              // ref/query 是否为"相反链"
        // 可选：预先缓存对角线，chaining 常用
        // int32_t diag{}; // (int32_t)pos_ref - (int32_t)pos_qry
    };

    // =====================================================================
    // collect_anchors - 收集锚点（参考 minimap2 实现）
    // ---------------------------------------------------------------------
    // 功能：
    // 从 ref_hits 和 qry_hits 中收集所有共享相同 hash 的锚点对
    //
    // 算法流程（minimap2 风格）：
    // 1. 对 ref_hits 按 hash 排序，构建 hash -> position 索引
    // 2. 遍历 qry_hits，通过 hash 查找匹配的 ref_hits
    // 3. 对于每个匹配，生成一个 Anchor
    //
    // 参数：
    //   ref_hits : 参考序列的种子命中列表
    //   qry_hits : 查询序列的种子命中列表
    //
    // 返回：
    //   std::vector<Anchor> - 锚点列表（未排序）
    //
    // 复杂度：
    //   时间：O(R log R + Q * avg_matches)
    //   空间：O(R + A)
    //
    // 注意：模板声明在这里，显式实例化在 chain.cpp 中
    // =====================================================================
    template <typename HitT>
    std::vector<Anchor> collect_anchors(const std::vector<HitT>& ref_hits, const std::vector<HitT>& qry_hits);

    // =====================================================================
    // sortAnchorsByDiagonal - 按对角线排序锚点（用于链化算法）
    // ---------------------------------------------------------------------
    // 排序规则：
    // 1. 首先按 rid_ref 排序
    // 2. 其次按对角线 (pos_ref - pos_qry) 排序
    // 3. 最后按 pos_ref 排序
    // =====================================================================
    void sortAnchorsByDiagonal(std::vector<Anchor>& anchors);

    // =====================================================================
    // sortAnchorsByPosition - 按位置排序锚点
    // ---------------------------------------------------------------------
    // 排序规则：
    // 1. 首先按 rid_ref 排序
    // 2. 其次按 pos_ref 排序
    // 3. 最后按 pos_qry 排序
    // =====================================================================
    void sortAnchorsByPosition(std::vector<Anchor>& anchors);

    // =====================================================================
    // filterHighFrequencyAnchors - 过滤高频锚点（参考 minimap2）
    // ---------------------------------------------------------------------
    // 功能：
    // 移除出现次数超过阈值的锚点（重复区域噪声过滤）
    //
    // 参数：
    //   anchors : 输入/输出锚点列表（原地修改）
    //   max_occ : 单个 hash 允许的最大出现次数（默认 500）
    // =====================================================================
    void filterHighFrequencyAnchors(std::vector<Anchor>& anchors, std::size_t max_occ = 500);


} // namespace seed

// ================================================================
namespace minimizer
{


    // ------------------------- minimizer hit（带位置） -------------------------
    // 用于 chaining/定位：需要知道落在哪条序列(rid)、什么位置(pos)、方向(strand)以及覆盖长度(span)。
    //
    // 内存优化：采用 minimap2 的 mm128_t 打包方式，只用 16 字节存下所有信息：
    //   x = (hash << 8) | span
    //     - bit[0..7]   : span (8 bit)
    //     - bit[8..63]  : hash (56 bit)
    //   y = (rid_with_strand << 32) | pos
    //     - bit[0..31]  : pos (32 bit)
    //     - bit[32..62] : rid (31 bit)
    //     - bit[63]     : strand (1 bit)
    //
    // 备注：
    // - 56bit hash 对绝大多数场景足够；如果你希望保留 64bit hash，需要改变打包布局。
    // - span 用 8bit，适合 k<=255；若以后需要更大 k，可调整布局。
    struct MinimizerHit : public seed::SeedHitBase<MinimizerHit>
    {
        hash_t x{0};
        hash_t y{0};

        constexpr MinimizerHit() = default;
        constexpr MinimizerHit(hash_t x_, hash_t y_) : x(x_), y(y_) {}

        // 打包/解包 x
        static constexpr hash_t pack_x(hash_t hash56, std::uint8_t span) noexcept
        {
            return (hash56 << 8) | static_cast<hash_t>(span);
        }
        static constexpr std::uint8_t span_from_x(hash_t x) noexcept
        {
            return static_cast<std::uint8_t>(x & 0xffULL);
        }
        static constexpr hash_t hash_from_x(hash_t x) noexcept
        {
            return (x >> 8);
        }

        // 打包/解包 y
        static constexpr hash_t pack_y(std::uint32_t pos, std::uint32_t rid, bool strand) noexcept
        {
            const std::uint32_t rid_with_strand = (rid & 0x7fffffffU) | (strand ? 0x80000000U : 0U);
            return (static_cast<hash_t>(rid_with_strand) << 32) | static_cast<hash_t>(pos);
        }
        static constexpr std::uint32_t pos_from_y(hash_t y) noexcept
        {
            return static_cast<std::uint32_t>(y & 0xffffffffULL);
        }
        static constexpr std::uint32_t rid_with_strand_from_y(hash_t y) noexcept
        {
            return static_cast<std::uint32_t>((y >> 32) & 0xffffffffULL);
        }
        static constexpr std::uint32_t rid_from_y(hash_t y) noexcept
        {
            return rid_with_strand_from_y(y) & 0x7fffffffU;
        }
        static constexpr bool strand_from_y(hash_t y) noexcept
        {
            return (rid_with_strand_from_y(y) & 0x80000000U) != 0U;
        }

        // 便捷构造：输入“完整语义字段”，内部自动打包
        constexpr MinimizerHit(hash_t hash56, std::uint32_t pos, std::uint32_t rid, bool strand, std::uint8_t span) noexcept
            : x(pack_x(hash56, span)), y(pack_y(pos, rid, strand))
        {
        }

        // SeedHitBase required API（仍然提供 hash/pos/rid/strand/span）
        constexpr hash_t hash() const noexcept { return hash_from_x(x); }
        constexpr std::uint32_t pos() const noexcept { return pos_from_y(y); }
        constexpr std::uint32_t rid() const noexcept { return rid_from_y(y); }
        constexpr bool strand() const noexcept { return strand_from_y(y); }
        constexpr std::uint32_t span() const noexcept { return span_from_x(x); }
    };

    static_assert(sizeof(MinimizerHit) == 16, "MinimizerHit should be 16 bytes when packed as (x,y)");
    using MinimizerHits = std::vector<MinimizerHit>;

    // =====================================================================
    // extractMinimizer
    // ---------------------------------------------------------------------
    // 从一条输入序列中提取 minimizer hit 列表（hash+位置）。
    //
    // 参数：
    //   seq       : 输入序列
    //   k         : k-mer 大小
    //   w         : 窗口大小（以 k-mer 为单位）
    //
    // 返回：
    //   minimizer hit 列表（按扫描顺序）。
    //
    // 注意：这里只提供声明；实现请放到对应的 .cpp 文件中。
    // =====================================================================
    MinimizerHits extractMinimizer(const std::string& seq,
                                   std::size_t k,
                                   std::size_t w,
                                   bool non_canonical);


    // =============================================================
    // nt4_table
    // -------------------------------------------------------------
    // 把输入字符映射为 2-bit 编码（0..3），其余字符为 invalid(4)。
    // - 'A'/'a' -> 0
    // - 'C'/'c' -> 1
    // - 'G'/'g' -> 2
    // - 'T'/'t' -> 3
    // - 'U'/'u' -> 3   (RNA 的 U 当作 T)
    // - 其它（包括 'N'/'n'、'-' 等） -> 4
    //
    // 设计目的：
    // - 让高性能实现（rolling k-mer/minimap2风格）可以 O(1) 查表，避免 switch 分支。
    // - 表放头文件里，多个 .cpp 复用，不用重复写 256 项初始化。
    // =============================================================
    inline constexpr std::uint8_t nt4_table[256] = {
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,3,3,4,4,
        4,4,4,4,4,4,4,4,4,4,0,4,1,4,4,4,
        2,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,
        3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4
    };

} // namespace minimizer


#endif //HALIGN4_SEED_H
