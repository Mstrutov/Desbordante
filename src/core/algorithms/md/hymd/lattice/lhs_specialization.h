#pragma once

#include "algorithms/md/hymd/md_lhs.h"

namespace algos::hymd::lattice {
struct SpecializationData {
    MdLhs::iterator const spec_before;
    LhsNode const new_child;
};

struct LhsSpecialization {
    using Unspecialized = MdLhs const&;

    MdLhs const& old_lhs;
    SpecializationData const specialization_data;

    Unspecialized ToUnspecialized() const {
        return old_lhs;
    }

    LhsSpecialization const& GetLhsSpecialization() const {
        return *this;
    }
};
}  // namespace algos::hymd::lattice
