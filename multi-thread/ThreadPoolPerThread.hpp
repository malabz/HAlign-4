#pragma once
#include <vector>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <atomic>

class ThreadPoolPerThread {
public:
    explicit ThreadPoolPerThread(size_t num_threads);
    ~ThreadPoolPerThread();

    // 将任务投递给 tid 所属线程
    void enqueue(int tid, std::function<void()> task);

    // 等待所有任务完成
    void wait_for_all();

private:
    struct Worker {
        std::queue<std::function<void()>> tasks;
        std::mutex mutex;
        std::condition_variable cond;
        bool stop = false;
        std::thread thread;
        int active_tasks = 0;
    };

    std::vector<Worker> workers;
    std::atomic<size_t> total_active;
};
