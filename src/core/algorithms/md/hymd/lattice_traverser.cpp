#include "algorithms/md/hymd/lattice_traverser.h"

#include "algorithms/md/hymd/lowest_bound.h"
#include "algorithms/md/hymd/utility/set_for_scope.h"
#include "model/index.h"

namespace algos::hymd {

bool LatticeTraverser::TraverseLattice(bool const traverse_all) {
    using model::Index;
    while (level_getter_->AreLevelsLeft()) {
        std::vector<lattice::ValidationInfo> validations = level_getter_->GetCurrentMds();
        if (validations.empty()) {
            continue;
        }

        std::size_t const validations_size = validations.size();
        std::vector<Validator::Result> results(validations_size);
        auto validate_at_index = [&](Index i) { results[i] = validator_.Validate(validations[i]); };
        pool_->ExecIndex(validate_at_index, validations_size);
        pool_->WorkUntilComplete();
        auto viol_func = [this, &results]() {
            for (Validator::Result& result : results) {
                for (std::vector<Recommendation> const& rhs_violations : result.recommendations) {
                    recommendations_.insert(rhs_violations.begin(), rhs_violations.end());
                };
            }
        };
        pool_->ExecSingle(viol_func);
        for (Index i = 0; i < validations_size; ++i) {
            Validator::Result const& result = results[i];
            lattice::MdLattice::MdVerificationMessenger& messenger = *validations[i].messenger;
            if (result.is_unsupported) {
                messenger.MarkUnsupported();
            } else {
                messenger.LowerAndSpecialize(result.invalidated);
            }
        }
        pool_->WorkUntilComplete();
        if (!traverse_all) return false;
    }
    return true;
}

}  // namespace algos::hymd
