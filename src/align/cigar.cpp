#include "align.h"

#include <cassert>
#include <cstdint>
#include <stdexcept>

namespace cigar {

// ------------------------------------------------------------------
// 约定：SAM CIGAR 在 BAM 中的压缩编码方式为：
//   cigar_unit = (len << 4) | op
// 其中 op 为 4-bit 操作码（0..15），len 为 28-bit（0..2^28-1）。
//
// 本文件提供：
// - cigarToInt : (char op, len) -> CigarUnit
// - intToCigar : CigarUnit -> (char op, len)
//
// 注意：
// - 这里的编码/解码只做“表示层”的转换，不对 len 做业务意义校验（例如 len==0 是否允许）。
// - 对于未知操作符，我们选择抛出异常（比静默映射成 '?' 更安全）。
// ------------------------------------------------------------------

namespace {

// SAM/BAM 的 CIGAR 操作码（4-bit）。
// 这套编号是行业约定（htslib/BAM spec），必须固定，否则与外部工具不兼容。
constexpr std::uint32_t kOpM  = 0;  // M: alignment match (match or mismatch)
constexpr std::uint32_t kOpI  = 1;  // I: insertion to the reference
constexpr std::uint32_t kOpD  = 2;  // D: deletion from the reference
constexpr std::uint32_t kOpN  = 3;  // N: skipped region from the reference
constexpr std::uint32_t kOpS  = 4;  // S: soft clipping (clipped sequences present in SEQ)
constexpr std::uint32_t kOpH  = 5;  // H: hard clipping (clipped sequences NOT present in SEQ)
constexpr std::uint32_t kOpP  = 6;  // P: padding (silent deletion from padded reference)
constexpr std::uint32_t kOpEq = 7;  // =: sequence match
constexpr std::uint32_t kOpX  = 8;  // X: sequence mismatch

constexpr std::uint32_t kOpMask = 0xFu;
constexpr std::uint32_t kLenBits = 28u;
constexpr std::uint32_t kMaxLen = (1u << kLenBits) - 1u;

// 把字符 op 映射为 4-bit opCode。
// 说明：用 switch 是最直接且零开销（编译器可生成跳表/比较链）。
constexpr std::uint32_t opCharToCode(const char operation)
{
    switch (operation) {
    case 'M': return kOpM;
    case 'I': return kOpI;
    case 'D': return kOpD;
    case 'N': return kOpN;
    case 'S': return kOpS;
    case 'H': return kOpH;
    case 'P': return kOpP;
    case '=': return kOpEq;
    case 'X': return kOpX;
    default:  return 0xFFFFFFFFu;
    }
}

// 把 4-bit opCode 映射回字符 op。
constexpr char opCodeToChar(const std::uint32_t op_code)
{
    switch (op_code) {
    case kOpM:  return 'M';
    case kOpI:  return 'I';
    case kOpD:  return 'D';
    case kOpN:  return 'N';
    case kOpS:  return 'S';
    case kOpH:  return 'H';
    case kOpP:  return 'P';
    case kOpEq: return '=';
    case kOpX:  return 'X';
    default:    return '?';
    }
}

} // namespace

CigarUnit cigarToInt(const char operation, const std::uint32_t len)
{
    // 关键正确性：len 只能占用 28bit。超出会破坏 opCode 位域。
    // 这里选择抛异常，避免 silent overflow。
    if (len > kMaxLen) {
        throw std::invalid_argument("cigarToInt: len exceeds 28-bit limit");
    }

    const std::uint32_t op_code = opCharToCode(operation);
    if (op_code == 0xFFFFFFFFu) {
        throw std::invalid_argument("cigarToInt: unknown CIGAR operation");
    }

    // (len << 4) | op
    return static_cast<CigarUnit>((len << 4) | (op_code & kOpMask));
}

void intToCigar(const CigarUnit cigar_unit, char& operation, std::uint32_t& len)
{
    // 低 4 位是 opCode，高 28 位是长度。
    const std::uint32_t op_code = static_cast<std::uint32_t>(cigar_unit) & kOpMask;
    len = static_cast<std::uint32_t>(cigar_unit) >> 4;

    operation = opCodeToChar(op_code);

    // 对于未知 opcode：这里用断言 + 保守返回 '?'。
    // - release 模式下不中断（避免解析外部脏数据时直接崩溃）
    // - debug 模式下尽早暴露问题
    assert(operation != '?');
}

// ------------------------------------------------------------------
// 函数：hasInsertion
// 功能：检测 CIGAR 序列中是否存在插入操作（'I'）
//
// 实现说明：
// 1. 性能优化：直接对 CigarUnit 的低 4 位进行位运算检查，避免完整解码
// 2. 短路优化：找到第一个 'I' 操作即返回 true，无需遍历整个序列
// 3. 复杂度：O(N)，最坏情况遍历所有操作；最好情况 O(1)（第一个就是 'I'）
// 4. 正确性：依赖 kOpI = 1 的编码约定（SAM/BAM 规范固定值）
//
// 参数：cigar - CIGAR 操作序列（压缩形式，每个元素为 CigarUnit）
// 返回：true 表示存在至少一个插入操作，false 表示不存在
// ------------------------------------------------------------------
bool hasInsertion(const Cigar_t& cigar)
{
    // 遍历所有 CIGAR 操作，检查低 4 位是否为 kOpI（插入操作码 = 1）
    for (const auto cigar_unit : cigar) {
        // 提取低 4 位的操作码
        const std::uint32_t op_code = static_cast<std::uint32_t>(cigar_unit) & kOpMask;

        // 如果是插入操作，立即返回 true（短路优化）
        if (op_code == kOpI) {
            return true;
        }
    }

    // 遍历完所有操作都没有找到插入，返回 false
    return false;
}

// ------------------------------------------------------------------
// 函数：cigarToString
// 功能：将 CIGAR 从压缩格式转换为 SAM 标准字符串格式
//
// 实现说明：
// 1. 性能优化：预分配字符串空间（cigar.size() * 5），避免多次内存重新分配
//    - 估算：平均每个 CIGAR 操作占用约 3-5 字符（如 "100M" = 4 字符）
//    - 对于大多数情况，reserve(cigar.size() * 5) 足够，避免了 string 的自动扩容
// 2. 使用 std::to_string 转换整数，编译器会优化为高效实现
// 3. 复杂度：O(N)，N 为 CIGAR 操作数量
// 4. 字符串拼接：使用 += 操作符，在预分配空间足够时为 O(1) 追加
//
// 参数：cigar - CIGAR 操作序列（压缩形式）
// 返回：SAM 格式的 CIGAR 字符串，例如 "100M5I95M"
//
// 示例：
//   输入：[cigarToInt('M', 100), cigarToInt('I', 5), cigarToInt('M', 95)]
//   输出："100M5I95M"
// ------------------------------------------------------------------
std::string cigarToString(const Cigar_t& cigar)
{
    // 预分配字符串空间，减少内存重新分配次数
    // 估算：每个 CIGAR 操作平均占用约 5 个字符（例如 "1000M" = 5 字符）
    std::string result;
    result.reserve(cigar.size() * 5);

    // 遍历所有 CIGAR 操作，逐个解码并拼接到结果字符串
    for (const auto cigar_unit : cigar) {
        char op;
        std::uint32_t len;
        intToCigar(cigar_unit, op, len);

        // 格式：<长度><操作符>，例如 "100M"
        result += std::to_string(len);
        result += op;
    }

    return result;
}

// ------------------------------------------------------------------
// 函数：alignQueryToRef
// 功能：根据 CIGAR 操作对 query 序列插入 gap 字符（'-'），使其与参考序列对齐
//
// **关键特性：保留原有 gap 字符**
// - 输入 query 中已存在的 gap 字符（'-'）会被当作普通碱基处理（不删除、不特殊对待）
// - 只根据 CIGAR 操作在适当位置插入新的 gap
// - 例如：query = "A-CG-T"，CIGAR = "2M1D2M2M" 时：
//   * 2M 消耗 "A-"（包括原有 gap），拷贝到输出
//   * 1D 插入新 gap "-"
//   * 2M 消耗 "CG"，拷贝到输出
//   * 2M 消耗 "-T"（包括原有 gap），拷贝到输出
//   * 结果："A--CG-T"（长度 7 = 6 + 1个新gap）
//
// 实现说明：
// 1. 预计算对齐后的总长度（遍历 CIGAR，累加所有操作对应的对齐序列长度）
// 2. 使用"从后往前"填充策略，避免频繁插入导致的 O(N²) 复杂度
// 3. 算法核心：
//    - 先将 query 移动到临时变量（避免拷贝，O(1) 复杂度）
//    - 将 query resize 到最终长度，预填充为 gap 字符 '-'（O(N) 复杂度，单次分配）
//    - 从后往前遍历 CIGAR，根据操作类型决定是拷贝原字符（包括原有 gap）还是保持新 gap
// 4. 复杂度：O(M + N)，M 为 CIGAR 操作数，N 为 query 长度（包括原有 gap）
//
// CIGAR 操作处理逻辑：
// - M/=/X: query 和 ref 都消耗，拷贝原 query 字符（**包括原有 gap**）
// - I: query 相对 ref 的插入，只消耗 query，拷贝原 query 字符（**包括原有 gap**）
// - D: query 相对 ref 的缺失，只消耗 ref，保持新插入的 gap 字符（已预填充）
// - S: soft clip，拷贝原 query 字符（**包括原有 gap**）
// - H: hard clip，不处理（序列已不存在）
//
// 性能优化：
// - 只遍历 CIGAR 两次（第一次计算长度 O(M)，第二次填充 O(M+N)）
// - 单次内存分配（assign），避免多次 resize/reserve
// - 从后往前填充，利用局部性原理，cache 友好
// - 移动语义减少内存拷贝
// ------------------------------------------------------------------
void alignQueryToRef(std::string& query, const Cigar_t& cigar)
{
    // 边界情况：空序列或空 CIGAR，直接返回
    if (query.empty() || cigar.empty()) {
        return;
    }
    // ------------------------------------------------------------------
    // 步骤1：计算对齐后的总长度
    // 说明：遍历 CIGAR，累加所有会消耗 ref 位置的操作长度
    //       - M/=/X/D: 消耗 ref，累加到 aligned_length
    //       - I/S/H: 不消耗 ref，也需要累加（因为 query 中存在）
    // ------------------------------------------------------------------
    std::size_t aligned_length = 0;
    std::size_t query_consumed = 0;  // 记录 CIGAR 会消耗的 query 长度（用于校验）

    for (const auto c : cigar) {
        char op;
        std::uint32_t len;
        intToCigar(c, op, len);

        // 计算对齐后的长度（包括 gap）
        if (op == 'M' || op == '=' || op == 'X' || op == 'D' || op == 'I' || op == 'S') {
            aligned_length += len;
        }
        // H 操作不占用空间（已从序列中移除）

        // 记录会消耗 query 的操作（用于后续校验）
        if (op == 'M' || op == '=' || op == 'X' || op == 'I' || op == 'S') {
            query_consumed += len;
        }
    }

    // 边界情况：空 CIGAR（所有操作都是 H），直接返回
    if (aligned_length == 0) {
        return;
    }

    // ------------------------------------------------------------------
    // 校验：CIGAR 消耗的 query 长度必须与原始 query 长度一致
    // 说明：
    // - 原始 query 包含已存在的 gap 字符（'-'），这些 gap 也会被 CIGAR 消耗
    // - 例如：query = "A-CG"，CIGAR = "4M" 是合法的（消耗 4 个字符，包括 1 个原有 gap）
    // - 如果长度不匹配，说明 CIGAR 与 query 不一致（可能是外部调用错误）
    // ------------------------------------------------------------------
    #ifdef _DEBUG
    assert(query_consumed == query.size() &&
           "CIGAR 操作消耗的 query 长度与原始长度不匹配（包括原有 gap）");
    #endif

    // ------------------------------------------------------------------
    // 步骤2：从后往前构建对齐序列（避免频繁插入导致的 O(N²) 复杂度）
    // 说明：
    // 1. 先将 query 移动到临时变量（避免拷贝，利用移动语义）
    // 2. 将 query resize 到最终长度，预填充为 gap 字符 '-'（新插入的 gap）
    // 3. 从后往前遍历 CIGAR，根据操作类型决定是拷贝原字符（包括原有 gap）还是保持新 gap
    // 4. aligned_pos：对齐序列的写入位置（从后往前递减）
    // 5. query_pos：原始 query 的读取位置（从后往前递减）
    // 注意：原有的 gap 字符（'-'）会被当作普通字符拷贝
    // ------------------------------------------------------------------
    const std::string original_query = std::move(query);  // 移动语义，避免拷贝
    query.assign(aligned_length, '-');                    // 预分配空间，填充新 gap

    std::size_t aligned_pos = aligned_length;             // 对齐序列的写入位置（从末尾开始）
    std::size_t query_pos = original_query.size();        // 原始序列的读取位置（从末尾开始）

    // ------------------------------------------------------------------
    // 步骤3：从后往前遍历 CIGAR，填充对齐序列
    // 复杂度：O(M + N)，M 为 CIGAR 操作数，N 为 query 长度（包括原有 gap）
    // 说明：原有的 gap 字符（'-'）会被当作普通字符处理，与 A/C/G/T 等碱基无区别
    // ------------------------------------------------------------------
    for (auto it = cigar.rbegin(); it != cigar.rend(); ++it) {
        char op;
        std::uint32_t len;
        intToCigar(*it, op, len);

        if (op == 'M' || op == '=' || op == 'X') {
            // Match/Mismatch：拷贝原 query 字符（包括原有 gap）
            // 说明：这些操作消耗 query 和 ref，将原序列字符拷贝到对齐序列
            //       如果原序列包含 gap 字符（'-'），也会被拷贝
            #ifdef _DEBUG
            assert(query_pos >= len && "CIGAR M/=/X 操作超出 query 长度");
            assert(aligned_pos >= len && "对齐位置越界");
            #endif

            for (std::uint32_t i = 0; i < len; ++i) {
                query[--aligned_pos] = original_query[--query_pos];
            }
        } else if (op == 'I') {
            // Insertion：拷贝原 query 字符（包括原有 gap）
            // 说明：query 相对 ref 的插入，只消耗 query，将原序列字符拷贝到对齐序列
            //       如果原序列包含 gap 字符（'-'），也会被拷贝
            #ifdef _DEBUG
            assert(query_pos >= len && "CIGAR I 操作超出 query 长度");
            assert(aligned_pos >= len && "对齐位置越界");
            #endif

            for (std::uint32_t i = 0; i < len; ++i) {
                query[--aligned_pos] = original_query[--query_pos];
            }
        } else if (op == 'D') {
            // Deletion：保持新插入的 gap 字符
            // 说明：query 相对 ref 的缺失，不消耗 query，跳过（已在 assign 时预填充 '-'）
            //       这里插入的是新 gap，与原序列中已存在的 gap 无关
            #ifdef _DEBUG
            assert(aligned_pos >= len && "CIGAR D 操作超出对齐序列长度");
            #endif

            aligned_pos -= len;  // 跳过新插入的 gap 字符（已在 assign 时填充）
        } else if (op == 'S') {
            // Soft Clipping：拷贝原 query 字符（包括原有 gap）
            // 说明：query 中存在但未比对的部分，将原序列字符拷贝到对齐序列
            //       如果原序列包含 gap 字符（'-'），也会被拷贝
            #ifdef _DEBUG
            assert(query_pos >= len && "CIGAR S 操作超出 query 长度");
            assert(aligned_pos >= len && "对齐位置越界");
            #endif

            for (std::uint32_t i = 0; i < len; ++i) {
                query[--aligned_pos] = original_query[--query_pos];
            }
        } else if (op == 'H') {
            // Hard Clipping：不处理
            // 说明：已从 query 中移除，不消耗 query，也不插入字符
            // 无操作
        } else if (op == 'N' || op == 'P') {
            // N (skipped region) / P (padding)：当作 deletion 处理
            // 说明：这些操作在典型的 alignment 中很少见，保守处理为新 gap
            #ifdef _DEBUG
            assert(aligned_pos >= len && "CIGAR N/P 操作超出对齐序列长度");
            #endif

            aligned_pos -= len;  // 保持新插入的 gap 字符
        } else {
            // 未知操作符：抛出异常
            throw std::runtime_error(
                std::string("alignQueryToRef: 不支持的 CIGAR 操作符 '") + op + "'"
            );
        }
    }

    // ------------------------------------------------------------------
    // 步骤4：验证对齐结果（仅在 Debug 模式下）
    // 说明：确保所有 query 字符（包括原有 gap）都已被处理，且对齐序列已完全填充
    // ------------------------------------------------------------------
    #ifdef _DEBUG
    assert(query_pos == 0 && "CIGAR 操作未完全消耗 query 序列（包括原有 gap）");
    assert(aligned_pos == 0 && "对齐序列未完全填充");
    #endif
}

} // namespace cigar

