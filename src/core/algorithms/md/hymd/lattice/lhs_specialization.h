#pragma once

#include "algorithms/md/hymd/md_lhs.h"

namespace algos::hymd::lattice {
struct SpecializationData {
    MdLhs::iterator spec_before;
    LhsNode new_child;
};

struct LhsSpecialization {
    using Unspecialized = MdLhs const&;

    MdLhs const& old_lhs;
    SpecializationData specialization_data;

    Unspecialized ToUnspecialized() const {
        return old_lhs;
    }

    LhsSpecialization const& GetLhsSpecialization() const {
        return *this;
    }
};
}  // namespace algos::hymd::lattice
