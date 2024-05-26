#pragma once

#include "algorithms/md/hymd/lattice/lhs_specialization.h"
#include "algorithms/md/hymd/lattice/md.h"
#include "algorithms/md/hymd/md_element.h"

namespace algos::hymd::lattice {
struct MdSpecialization {
    using Unspecialized = Md;

    LhsSpecialization const& lhs_specialization;
    MdElement rhs;

    Unspecialized ToUnspecialized() const noexcept {
        return {lhs_specialization.ToUnspecialized(), rhs};
    }

    LhsSpecialization const& GetLhsSpecialization() const noexcept {
        return lhs_specialization;
    }
};
}  // namespace algos::hymd::lattice
