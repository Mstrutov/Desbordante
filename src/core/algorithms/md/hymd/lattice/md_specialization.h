#pragma once

#include "algorithms/md/hymd/lattice/lhs_specialization.h"
#include "algorithms/md/hymd/lattice/md.h"
#include "algorithms/md/hymd/md_element.h"

namespace algos::hymd::lattice {
struct MdSpecialization {
    LhsSpecialization const& lhs_specialization;
    MdElement rhs;

    Md ToUnspecialized() const {
        return {lhs_specialization.ToUnspecialized(), rhs};
    }

    LhsSpecialization const& GetLhsSpecialization() const {
        return lhs_specialization;
    }
};
}  // namespace algos::hymd::lattice
