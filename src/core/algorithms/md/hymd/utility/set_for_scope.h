#pragma once

#include "algorithms/md/decision_boundary.h"

namespace algos::hymd::utility {
class SetForScope {
    model::md::DecisionBoundary& ref_;
    model::md::DecisionBoundary const old_value_;

public:
    SetForScope(model::md::DecisionBoundary& ref,
                model::md::DecisionBoundary const new_value) noexcept
        : ref_(ref), old_value_(ref) {
        ref_ = new_value;
    }

    ~SetForScope() {
        ref_ = old_value_;
    }
};
}  // namespace algos::hymd::utility
