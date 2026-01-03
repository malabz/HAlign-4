#include "consensus.h"

// ---------------- TopKLongestSelector 实现说明 ----------------
// TopKLongestSelector 用于在单次顺序扫描中选出长度最长的 K 条序列。
// 设计目标：
// - 空间复杂度 O(K)，只保留 K 条记录（使用堆存储）；
// - 时间复杂度：每个元素的插入/维护为 O(log K)；遍历 N 条记录总体为 O(N log K)；
// - 稳定性：对相同长度的序列，优先保留较早出现的序列（即出现更早的序列被视为更“好”）。
//
// 实现细节：
// - 使用一个最小堆（min-heap），heap_[0] 始终为当前已收集 K 条中的最差（即最短或最晚出现的）元素；
// - 当堆未满（size < K）时直接加入并上浮（siftUp）；
// - 当堆已满时，比较候选项与 heap_[0]：如果候选项更好（长度更长或相同长度但出现更早），则替换 heap_[0] 并下沉（siftDown）；
// - 比较规则（betterThan / worseThan）使用 (length, order) 的字典序：长度为主，出现顺序为次（更早更优）。

TopKLongestSelector::TopKLongestSelector(std::size_t k)
        : k_(k)
    {
        heap_.reserve(k_);
    }

    void TopKLongestSelector::reset(std::size_t k)
    {
        k_ = k;
        order_counter_ = 0; // order_counter_ 用于记录元素出现的先后顺序，以便在长度相同时保持稳定性
        heap_.clear();
        heap_.reserve(k_);
    }

    std::size_t TopKLongestSelector::size() const
    {
        return heap_.size();
    }

    std::size_t TopKLongestSelector::capacity() const
    {
        return k_;
    }

    bool TopKLongestSelector::empty() const
    {
        return heap_.empty();
    }

    // 比较函数：判断 a 是否比 b 更 "差"
    // 语义：在最小堆中，“更差”的元素会被置于堆顶（即应当被替换或弹出）。
    bool TopKLongestSelector::worseThan(const Item& a, const Item& b)
    {
        if (a.len != b.len) return a.len < b.len;      // 长度更短的被认为更差
        return a.order > b.order;                      // 同长度时，后来出现的（order 值更大）被认为更差
    }

    // 判断候选项 cand 是否比当前最坏项 worst 更好（用于决定是否替换堆顶）
    bool TopKLongestSelector::betterThan(const Item& cand, const Item& worst)
    {
        if (cand.len != worst.len) return cand.len > worst.len; // 更长则更好
        // 同长度时，我们偏向保留更早出现的（即 order 更小的项优先）
        return cand.order < worst.order;
    }

    // 上浮操作：当新元素插入到堆尾时调用，若其比父节点更差则与父节点交换直到堆序恢复
    // 说明：这里实现的是一个 min-heap 的上浮（维护堆顶为最差元素），使得堆顶始终为当前最差
    void TopKLongestSelector::siftUp(std::size_t idx)
    {
        while (idx > 0) {
            const std::size_t parent = (idx - 1) / 2;
            // min-heap：如果当前更“差”（更小），就上浮（与父节点交换）
            if (worseThan(heap_[idx], heap_[parent])) {
                std::swap(heap_[idx], heap_[parent]);
                idx = parent;
            } else {
                break;
            }
        }
    }

    // 下沉操作：当堆顶被替换后调用，把新的堆顶下沉到合适位置
    // 通过比较左右孩子找到更"差"的孩子并与当前交换，直到堆序恢复
    void TopKLongestSelector::siftDown(std::size_t idx)
    {
        const std::size_t n = heap_.size();
        while (true) {
            const std::size_t left = idx * 2 + 1;
            if (left >= n) break; // 无孩子，停止

            const std::size_t right = left + 1;

            // 选择更“差”的孩子（min-heap 中更差的孩子应当成为交换对象）
            std::size_t worst_child = left;
            if (right < n && worseThan(heap_[right], heap_[left])) {
                worst_child = right;
            }

            // 如果孩子比当前更“差”，则下沉（交换）
            if (worseThan(heap_[worst_child], heap_[idx])) {
                std::swap(heap_[idx], heap_[worst_child]);
                idx = worst_child;
            } else {
                break;
            }
        }
    }

    // 考虑一条新记录：将其包装为 Item 并尝试加入堆中
    // 逻辑：
    // - 如果 k_ == 0，直接返回（不保存任何记录）
    // - 如果堆未满，则直接 push 并上浮
    // - 否则比较候选与堆顶（最差）；若候选更好，则替换堆顶并下沉
    void TopKLongestSelector::consider(seq_io::SeqRecord rec)
    {
        if (k_ == 0) return;

        Item cand;
        cand.len = rec.seq.size();
        cand.order = order_counter_++;
        cand.rec = std::move(rec);

        if (heap_.size() < k_) {
            heap_.push_back(std::move(cand));
            siftUp(heap_.size() - 1);
            return;
        }

        // heap_[0] 是当前“最差”的那条（堆顶）
        if (betterThan(cand, heap_[0])) {
            heap_[0] = std::move(cand);
            siftDown(0);
        }
    }

    // 将保留的堆内容按长度降序输出（同长度按出现顺序升序，保证稳定性）
    // 返回值：vector<SeqRecord>，其中包含已选择的记录（按长度从长到短排序）
    std::vector<seq_io::SeqRecord> TopKLongestSelector::takeSortedDesc()
    {
        // 为避免复制大量序列数据，move 出 heap_ 的内容到一个临时 vector 上进行排序
        std::vector<Item> items = std::move(heap_);
        heap_.clear();
        heap_.shrink_to_fit();
        heap_.reserve(k_);

        // 自定义排序：长度降序；若长度相同则按原始出现顺序升序（更早的排在前面）
        std::sort(items.begin(), items.end(),
                  [](const Item& a, const Item& b) {
                      if (a.len != b.len) return a.len > b.len;
                      return a.order < b.order;
                  });

        // 将 SeqRecord 从 items 移出到返回向量中
        std::vector<seq_io::SeqRecord> out;
        out.reserve(items.size());
        for (auto& it : items) {
            out.push_back(std::move(it.rec));
        }
        return out;
    }
