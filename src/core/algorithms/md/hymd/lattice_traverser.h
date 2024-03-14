#pragma once

#include "algorithms/md/hymd/indexes/dictionary_compressor.h"
#include "algorithms/md/hymd/lattice/cardinality/min_picker_lattice.h"
#include "algorithms/md/hymd/lattice/level_getter.h"
#include "algorithms/md/hymd/lattice/md_lattice.h"
#include "algorithms/md/hymd/similarity_data.h"
#include "algorithms/md/hymd/validator.h"
#include "util/worker_thread_pool.h"

namespace algos::hymd {

class LatticeTraverser {
private:
    lattice::MdLattice* const lattice_;
    Recommendations recommendations_;

    std::unique_ptr<lattice::LevelGetter> const level_getter_;
    Validator const validator_;

    util::WorkerThreadPool* pool_;

public:
    LatticeTraverser(lattice::MdLattice* lattice,
                     std::unique_ptr<lattice::LevelGetter> level_getter, Validator validator,
                     util::WorkerThreadPool* pool) noexcept
        : lattice_(lattice),
          level_getter_(std::move(level_getter)),
          validator_(validator),
          pool_(pool) {}

    bool TraverseLattice(bool traverse_all);

    Recommendations TakeRecommendations() noexcept {
        auto recommendations = std::move(recommendations_);
        recommendations_.clear();
        return recommendations;
    }
};

}  // namespace algos::hymd
