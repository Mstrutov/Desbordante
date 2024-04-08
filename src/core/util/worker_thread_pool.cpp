#include "util/worker_thread_pool.h"

#include <cassert>

namespace util {
WorkerThreadPool::WorkerThreadPool(std::size_t thread_num)
    : barrier_(thread_num, Completion{*this}) {
    assert(thread_num > 1);
    worker_threads_.reserve(--thread_num);
    try {
        while (thread_num--) {
            worker_threads_.emplace_back(&WorkerThreadPool::Work, this);
        }
    } catch (std::system_error&) {
        Terminate();
        barrier_.wait(barrier_.arrive(++thread_num));
        throw;
    }
}
}  // namespace util
