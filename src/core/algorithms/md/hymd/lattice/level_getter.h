#pragma once

#include <cstddef>
#include <vector>

#include "algorithms/md/hymd/lattice/full_lattice.h"
#include "algorithms/md/hymd/lattice/validation_info.h"

namespace algos::hymd::lattice {

class LevelGetter {
protected:
    std::size_t cur_level_ = 0;
    FullLattice* const lattice_;
    // Prevent lifetime issues.
    std::vector<lattice::MdLatticeNodeInfo> lattice_level_info_;

    virtual std::vector<ValidationInfo> GetCurrentMdsInternal(
            std::vector<lattice::MdLatticeNodeInfo>& level_mds) = 0;

public:
    LevelGetter(FullLattice* lattice) : lattice_(lattice) {}

    bool AreLevelsLeft() const noexcept {
        return cur_level_ <= lattice_->GetMaxLevel();
    }

    std::vector<ValidationInfo> GetCurrentMds() {
        lattice_level_info_ = lattice_->GetLevel(cur_level_);
        return GetCurrentMdsInternal(lattice_level_info_);
    }

    virtual ~LevelGetter() = default;
};

}  // namespace algos::hymd::lattice
