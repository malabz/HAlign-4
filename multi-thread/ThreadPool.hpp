#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <vector>
#include <atomic>

class ThreadPool {
public:
    explicit ThreadPool(size_t num_threads);
    ~ThreadPool();

    void enqueue(std::function<void()> task);

    void wait_for_all();

private:
    std::vector<std::thread> workers;          // 工作线程集合
    std::queue<std::function<void()>> tasks;   // 任务队列（用队列存储，能够先进先出）

    std::mutex queue_mutex;                    // 任务队列互斥锁（用于线程同步）
    std::condition_variable condition;         // 条件变量（用于线程同步）
    std::atomic<bool> stop;                    // 线程池停止标志
    std::atomic<int> active_workers;
};
