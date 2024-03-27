#pragma once

#include <cassert>
#include <cstddef>
#include <vector>

#include "algorithms/md/decision_boundary.h"
#include "algorithms/md/hymd/decision_boundary_vector.h"
#include "algorithms/md/hymd/lattice/md_lattice.h"
#include "algorithms/md/hymd/lattice/md_lattice_node_info.h"
#include "algorithms/md/hymd/lattice/single_level_func.h"
#include "algorithms/md/hymd/lattice/support_lattice.h"
#include "model/index.h"
#include "util/erase_if_replace.h"

namespace algos::hymd::lattice {

class FullLattice {
private:
    MdLattice md_lattice_;
    SupportLattice support_lattice_;

public:
    std::size_t GetColMatchNumber() const noexcept {
        return md_lattice_.GetColMatchNumber();
    }

    [[nodiscard]] bool HasGeneralization(DecisionBoundaryVector const& lhs_bounds,
                                         model::md::DecisionBoundary const rhs_bound,
                                         model::Index const rhs_index) const {
        return md_lattice_.HasGeneralization(lhs_bounds, rhs_bound, rhs_index);
    }

    [[nodiscard]] std::size_t GetMaxLevel() const {
        return md_lattice_.GetMaxLevel();
    }

    std::vector<MdLatticeNodeInfo> GetLevel(std::size_t level) {
        // TODO: traverse both simultaneously.
        std::vector<MdLatticeNodeInfo> mds = md_lattice_.GetLevel(level);
        util::EraseIfReplace(mds, [this](MdLatticeNodeInfo const& node_info) {
            return support_lattice_.IsUnsupported(node_info.lhs_bounds);
        });
        return mds;
    }

    std::vector<MdLatticeNodeInfo> GetAll() {
        std::vector<MdLatticeNodeInfo> mds = md_lattice_.GetAll();
        util::EraseIfReplace(mds, [this](MdLatticeNodeInfo const& node_info) {
            return support_lattice_.IsUnsupported(node_info.lhs_bounds);
        });
        return mds;
    }

    std::vector<model::md::DecisionBoundary> GetRhsInterestingnessBounds(
            DecisionBoundaryVector const& lhs_bounds,
            std::vector<model::Index> const& indices) const {
        return md_lattice_.GetRhsInterestingnessBounds(lhs_bounds, indices);
    }

    void AddIfMinimal(DecisionBoundaryVector const& lhs_bounds,
                      model::md::DecisionBoundary const rhs_bound, model::Index const rhs_index) {
        assert(!support_lattice_.IsUnsupported(lhs_bounds));
        // TODO: Use info about what LHS was specialized from.
        md_lattice_.AddIfMinimal(lhs_bounds, rhs_bound, rhs_index);
    }

    void AddIfMinimalAndNotUnsupported(DecisionBoundaryVector const& lhs_bounds,
                                       model::md::DecisionBoundary const rhs_bound,
                                       model::Index const rhs_index) {
        // TODO: traverse both simultaneously.
        if (support_lattice_.IsUnsupported(lhs_bounds)) return;
        AddIfMinimal(lhs_bounds, rhs_bound, rhs_index);
    }

    std::vector<MdLatticeNodeInfo> FindViolated(SimilarityVector const& similarity_vector) {
        return md_lattice_.FindViolated(similarity_vector);
    }

    void MarkUnsupported(DecisionBoundaryVector const& lhs_bounds) {
        support_lattice_.MarkUnsupported(lhs_bounds);
    }

    bool IsUnsupported(DecisionBoundaryVector const& lhs_bounds) {
        return support_lattice_.IsUnsupported(lhs_bounds);
    }

    FullLattice(std::size_t column_matches_size, SingleLevelFunc single_level_func)
        : md_lattice_(column_matches_size, std::move(single_level_func)),
          support_lattice_(column_matches_size) {}
};

}  // namespace algos::hymd::lattice
