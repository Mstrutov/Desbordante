#include "timer.h"

namespace algos::fastod {

Timer::Timer(bool start) : is_started_(false) {
    if (start) {
        Start();
    }
}

void Timer::Start() {
    start_time_ = end_time_ = std::chrono::high_resolution_clock::now();
    is_started_ = true;
}

void Timer::Stop() {
    if (!is_started_) {
        return;
    }

    end_time_ = std::chrono::high_resolution_clock::now();
    is_started_ = false;
}

bool Timer::IsStarted() const {
    return is_started_;
}

double Timer::GetElapsedSeconds() const {
    std::chrono::nanoseconds elapsed_nanoseconds =
            is_started_ ? std::chrono::high_resolution_clock::now() - start_time_
                        : end_time_ - start_time_;

    return std::chrono::duration<double>(elapsed_nanoseconds).count();
}

}  // namespace algos::fastod
