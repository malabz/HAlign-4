#include "ThreadPoolPerThread.hpp"

ThreadPoolPerThread::ThreadPoolPerThread(size_t num_threads)
    : workers(num_threads), total_active(0) {
    for (size_t i = 0; i < num_threads; ++i) {
        workers[i].thread = std::thread([this, i]() {
            Worker& w = workers[i];
            while (true) {
                std::function<void()> task;
                {
                    std::unique_lock<std::mutex> lock(w.mutex);
                    w.cond.wait(lock, [&] {
                        return w.stop || !w.tasks.empty();
                    });

                    if (w.stop && w.tasks.empty())
                        break;

                    task = std::move(w.tasks.front());
                    w.tasks.pop();
                    ++w.active_tasks;
                    ++total_active;
                }

                task();

                {
                    std::lock_guard<std::mutex> lock(w.mutex);
                    --w.active_tasks;
                    --total_active;
                }
                w.cond.notify_all();
            }
        });
    }
}

void ThreadPoolPerThread::enqueue(int tid, std::function<void()> task) {
    Worker& w = workers[tid];
    {
        std::lock_guard<std::mutex> lock(w.mutex);
        w.tasks.emplace(std::move(task));
    }
    w.cond.notify_one();
}

void ThreadPoolPerThread::wait_for_all() {
    for (Worker& w : workers) {
        std::unique_lock<std::mutex> lock(w.mutex);
        w.cond.wait(lock, [&] {
            return w.tasks.empty() && w.active_tasks == 0;
        });
    }
}

ThreadPoolPerThread::~ThreadPoolPerThread() {
    for (auto& w : workers) {
        {
            std::lock_guard<std::mutex> lock(w.mutex);
            w.stop = true;
        }
        w.cond.notify_one();
    }

    for (auto& w : workers) {
        if (w.thread.joinable())
            w.thread.join();
    }
}
