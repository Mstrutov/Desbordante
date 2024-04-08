#pragma once

#include <atomic>
#include <barrier>
#include <cstddef>
#include <functional>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

#include "model/index.h"

namespace util {
class WorkerThreadPool {
public:
    // Both functions must be thread-safe. Worker must return false when no work is left.
    using Worker = std::function<void()>;

private:
    struct Completion {
        WorkerThreadPool& thread_pool;

        void operator()() noexcept {
            thread_pool.is_working_.clear();
        }
    };

    Worker work_;
    std::vector<std::jthread> worker_threads_;
    std::barrier<Completion> barrier_;
    std::atomic_flag is_working_;

    void Work() {
        while (true) {
            is_working_.wait(false);
            if (!work_) break;
            WorkUntilComplete();
        }
    }

    void Terminate() {
        SetWork(nullptr);
    }

public:
    WorkerThreadPool(std::size_t thread_num);

    // Previous tasks must be finished before calling this.
    template <typename WorkerType>
    void SetWork(WorkerType work) {
        work_ = std::move(work);
        is_working_.test_and_set();
        is_working_.notify_all();
    }

    template <typename FunctionType>
    void ExecSingle(FunctionType func) {
        SetWork([func, flag = std::make_shared<std::once_flag>()]() {
            std::call_once(*flag, func);
        });
    }

    template <typename FunctionType>
    void ExecIndex(FunctionType func, model::Index size) {
        auto work = [func, size, index = std::make_shared<std::atomic<model::Index>>()]() mutable {
            model::Index i;
            while ((i = (*index)++) < size) {
                func(i);
            }
            *index = size;
        };
        SetWork(work);
    }

    // Main thread must call this to finish.
    void WorkUntilComplete() {
        auto call_unchecked = [](Worker const& w) {
            if (!w) __builtin_unreachable();
            w();
        };
        call_unchecked(work_);
        barrier_.arrive_and_wait();
    }

    ~WorkerThreadPool() {
        Terminate();
    }
};
}  // namespace util
