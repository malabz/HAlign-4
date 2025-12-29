#include "ThreadPool.hpp"

ThreadPool::ThreadPool(size_t num_threads) : stop(false), active_workers(0) {
    for (size_t i = 0; i < num_threads; ++i) {
        workers.emplace_back([this]() {
            while (true) {
                std::function<void()> task;

                {
                    std::unique_lock<std::mutex> lock(this->queue_mutex);
                    this->condition.wait(lock, [this]() {
                        return this->stop || !this->tasks.empty();
                    });
                    if (this->stop && this->tasks.empty())
                        return;
                    task = std::move(this->tasks.front());
                    this->tasks.pop();
                    ++active_workers;
                }

                try {
                    task();
                } catch (...) {
                    std::lock_guard<std::mutex> lock(queue_mutex);
                    --active_workers;
                    condition.notify_all();
                    throw; // 或者处理异常
                }

                {
                    std::lock_guard<std::mutex> lock(queue_mutex);
                    --active_workers;
                    condition.notify_all();
                }
            }
        });
    }
}

void ThreadPool::enqueue(std::function<void()> task) {
    {
        std::lock_guard<std::mutex> lock(queue_mutex);
        tasks.emplace(std::move(task));
    }
    condition.notify_one();
}

void ThreadPool::wait_for_all() {
    std::unique_lock<std::mutex> lock(queue_mutex);
    condition.wait(lock, [this]() {
        return tasks.empty() && active_workers == 0;
    });
}

ThreadPool::~ThreadPool() {
    {
        std::lock_guard<std::mutex> lock(queue_mutex);
        stop = true;
    }
    condition.notify_all();
    for (std::thread &worker : workers)
        worker.join();
}
