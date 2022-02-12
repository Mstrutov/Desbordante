#include "primitive.h"

#include <cassert>

void algos::Primitive::AddProgress(double const val) noexcept {
    assert(val >= 0);
    std::scoped_lock lock(progress_mutex_);
    cur_phase_progress_ += val;
    assert(cur_phase_progress_ < 101);
}

void algos::Primitive::SetProgress(double const val) noexcept {
    assert(0 <= val && val < 101);
    std::scoped_lock lock(progress_mutex_);
    cur_phase_progress_ = val;
}

std::pair<uint8_t, double> algos::Primitive::GetProgress() const noexcept {
    std::scoped_lock lock(progress_mutex_);
    return std::make_pair(cur_phase_id_, cur_phase_progress_);
}

void algos::Primitive::ToNextProgressPhase() noexcept {
    /* Current phase done, ensure that this is displayed in the progress bar */
    SetProgress(kTotalProgressPercent);

    std::scoped_lock lock(progress_mutex_);
    ++cur_phase_id_;
    assert(cur_phase_id_ < phase_names_.size());
    cur_phase_progress_ = 0;
}
