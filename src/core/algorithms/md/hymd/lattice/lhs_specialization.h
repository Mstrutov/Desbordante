#pragma once

#include "algorithms/md/hymd/md_element.h"
#include "algorithms/md/hymd/md_lhs.h"

namespace algos::hymd::lattice {
struct LhsSpecialization {
    using Unspecialized = MdLhs const&;

    MdLhs const& old_lhs;
    MdElement specialized;

    Unspecialized ToUnspecialized() const {
        return old_lhs;
    }

    LhsSpecialization const& GetLhsSpecialization() const {
        return *this;
    }
};
}  // namespace algos::hymd::lattice
