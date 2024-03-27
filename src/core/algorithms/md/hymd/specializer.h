#pragma once

#include <optional>
#include <vector>

#include "algorithms/md/decision_boundary.h"
#include "algorithms/md/hymd/column_match_info.h"
#include "algorithms/md/hymd/decision_boundary_vector.h"
#include "algorithms/md/hymd/invalidated_rhs.h"
#include "algorithms/md/hymd/lattice/full_lattice.h"
#include "algorithms/md/hymd/similarity_vector.h"
#include "model/index.h"

namespace algos::hymd {

class Specializer {
private:
    std::vector<ColumnMatchInfo> const* const column_matches_info_;
    lattice::FullLattice* const lattice_;
    bool const prune_nondisjoint_;

    [[nodiscard]] std::optional<model::md::DecisionBoundary> SpecializeOneLhs(
            model::Index col_match_index, model::md::DecisionBoundary lhs_bound) const;

public:
    Specializer(std::vector<ColumnMatchInfo> const& column_matches_info,
                lattice::FullLattice* lattice, bool prune_nondisjoint)
        : column_matches_info_(&column_matches_info),
          lattice_(lattice),
          prune_nondisjoint_(prune_nondisjoint) {}

    void SpecializeFor(SimilarityVector const& sim, DecisionBoundaryVector& lhs_bounds,
                       model::Index rhs_index, model::md::DecisionBoundary old_rhs_bound);

    void SpecializeInvalidated(DecisionBoundaryVector& lhs_bounds,
                               InvalidatedRhss const& invalidated);
};

}  // namespace algos::hymd
