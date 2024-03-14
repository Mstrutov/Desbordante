#pragma once

#include <cstddef>
#include <vector>

#include "algorithms/md/hymd/lattice/md_lattice.h"
#include "algorithms/md/hymd/lattice/validation_info.h"

namespace algos::hymd::lattice {

class LevelGetter {
protected:
    std::size_t cur_level_ = 0;
    MdLattice* const lattice_;
    // Prevent lifetime issues.
    std::vector<MdLattice::MdVerificationMessenger> messengers_;

    virtual std::vector<ValidationInfo> GetCurrentMdsInternal(
            std::vector<MdLattice::MdVerificationMessenger>& level_mds) = 0;

public:
    LevelGetter(MdLattice* lattice) : lattice_(lattice) {}

    bool AreLevelsLeft() const noexcept {
        return cur_level_ <= lattice_->GetMaxLevel();
    }

    std::vector<ValidationInfo> GetCurrentMds() {
        messengers_ = lattice_->GetLevel(cur_level_);
        return GetCurrentMdsInternal(messengers_);
    }

    virtual ~LevelGetter() = default;
};

}  // namespace algos::hymd::lattice
