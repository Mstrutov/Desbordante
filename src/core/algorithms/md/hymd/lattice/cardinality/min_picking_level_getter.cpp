#include "algorithms/md/hymd/lattice/cardinality/min_picking_level_getter.h"

#include "algorithms/md/hymd/lowest_bound.h"
#include "util/erase_if_replace.h"

namespace algos::hymd::lattice::cardinality {

std::vector<ValidationInfo> MinPickingLevelGetter::GetCurrentMdsInternal(
        std::vector<MdLattice::MdVerificationMessenger>& level_lattice_info) {
    min_picker_.NewBatch(level_lattice_info.size());
    for (MdLattice::MdVerificationMessenger& messenger : level_lattice_info) {
        std::size_t const column_match_number = messenger.GetRhs().size();
        boost::dynamic_bitset<> const& previously_picked_rhs =
                picked_.try_emplace(messenger.GetLhs(), column_match_number).first->second;
        boost::dynamic_bitset<> indices(column_match_number);
        DecisionBoundaryVector const& rhs = messenger.GetRhs();
        for (model::Index i = 0; i < column_match_number; ++i) {
            if (rhs[i] != kLowestBound) {
                indices.set(i);
            }
        }
        indices -= previously_picked_rhs;
        if (indices.none()) continue;
        min_picker_.AddGeneralizations(messenger, indices);
    }
    std::vector<ValidationInfo> collected = min_picker_.GetAll();
    if constexpr (MinPickerType::kNeedsEmptyRemoval) {
        if constexpr (kEraseEmptyKeepOrder) {
            std::erase_if(collected, [](ValidationInfo const& validation_info) {
                return validation_info.rhs_indices.none();
            });
        } else {
            util::EraseIfReplace(collected, [](ValidationInfo const& validation_info) {
                return validation_info.rhs_indices.none();
            });
        }
    }
    for (ValidationInfo const& validation_info : collected) {
        std::size_t const column_match_number = validation_info.messenger->GetRhs().size();
        boost::dynamic_bitset<>& validated_indices =
                picked_.try_emplace(validation_info.messenger->GetLhs(), column_match_number)
                        .first->second;
        assert((validated_indices & validation_info.rhs_indices).none());
        validated_indices |= validation_info.rhs_indices;
    }
    if (collected.empty()) {
        picked_.clear();
        ++cur_level_;
    }
    return collected;
}

}  // namespace algos::hymd::lattice::cardinality
