#pragma once

#include "algorithms/md/hymd/indexes/dictionary_compressor.h"
#include "algorithms/md/hymd/lattice/cardinality/min_picker_lattice.h"
#include "algorithms/md/hymd/lattice/full_lattice.h"
#include "algorithms/md/hymd/lattice/level_getter.h"
#include "algorithms/md/hymd/similarity_data.h"
#include "algorithms/md/hymd/specializer.h"
#include "algorithms/md/hymd/validator.h"

namespace algos::hymd {

class LatticeTraverser {
private:
    lattice::FullLattice* const lattice_;
    Recommendations recommendations_;

    std::unique_ptr<lattice::LevelGetter> const level_getter_;
    Validator const validator_;

    Specializer* const specializer_;

    void LowerAndSpecialize(Validator::Result& validation_result,
                            lattice::ValidationInfo& validation_info);

public:
    LatticeTraverser(lattice::FullLattice* lattice,
                     std::unique_ptr<lattice::LevelGetter> level_getter, Validator validator,
                     Specializer* specializer) noexcept
        : lattice_(lattice),
          level_getter_(std::move(level_getter)),
          validator_(validator),
          specializer_(specializer) {}

    bool TraverseLattice(bool traverse_all);

    Recommendations TakeRecommendations() noexcept {
        auto recommendations = std::move(recommendations_);
        recommendations_.clear();
        return recommendations;
    }
};

}  // namespace algos::hymd
