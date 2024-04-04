#pragma once

#include "algorithms/md/hymd/decision_boundary_vector.h"
#include "algorithms/md/hymd/lattice/md.h"
#include "algorithms/md/hymd/lattice/md_specialization.h"
#include "algorithms/md/hymd/lattice/node_base.h"

namespace algos::hymd::lattice {
class MdNode : public NodeBase<MdNode> {
public:
    using Specialization = MdSpecialization;
    using Unspecialized = Md;

    DecisionBoundaryVector rhs_bounds;

    static MdLhs const& GetLhs(Unspecialized const& md) noexcept {
        return md.lhs;
    }

    bool ContainsGeneralizationOf(Unspecialized const& md) const noexcept {
        return rhs_bounds[md.rhs.index] >= md.rhs.decision_boundary;
    }

    MdNode* AddOneUnchecked(model::Index child_array_index, model::md::DecisionBoundary bound) {
        return AddOneUncheckedBase(child_array_index, bound, rhs_bounds.size());
    }

    MdNode(std::size_t attributes_num, std::size_t children_number)
        : NodeBase<MdNode>(children_number), rhs_bounds(attributes_num) {}

    explicit MdNode(DecisionBoundaryVector rhs)
        : NodeBase<MdNode>(rhs.size()), rhs_bounds(std::move(rhs)) {}
};
}  // namespace algos::hymd::lattice
