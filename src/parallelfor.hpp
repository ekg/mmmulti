#pragma once

#include <thread>
#include <cstdint>
#include <functional>
#include <atomic>
#include "atomic_queue.h"

namespace parallelfor {

using std::uint32_t;
using std::uint64_t;
using std::uint8_t;

template <typename I>
void parallelfor(const I& begin,
                 const I& end,
                 const uint64_t& nthreads,
                 const std::function<void(I)>& func) {
    auto queue_ptr = new atomic_queue::AtomicQueue2<I, 2 << 16>;
    auto& queue = *queue_ptr;
    std::atomic<bool> work_todo;
    auto worker =
        [&queue,&work_todo,&func](void) {
            I i;
            while (work_todo.load()) {
                if (queue.try_pop(i)) {
                    func(i);
                } else {
                    std::this_thread::sleep_for(std::chrono::nanoseconds(1));
                }
            }
        };
    std::vector<std::thread> workers;
    workers.reserve(nthreads);
    work_todo.store(true);
    for (uint64_t t = 0; t < nthreads; ++t) {
        workers.emplace_back(worker);
    }
    I todo_i = begin;
    while (todo_i != end) {
        if (queue.try_push(todo_i)) {
            ++todo_i;
        } else {
            std::this_thread::sleep_for(std::chrono::nanoseconds(1));
        }
    }
    while (!queue.was_empty()) {
        std::this_thread::sleep_for(std::chrono::nanoseconds(1));
    }
    work_todo.store(false);
    for (uint64_t t = 0; t < nthreads; ++t) {
        workers[t].join();
    }
    delete queue_ptr;
}

}
