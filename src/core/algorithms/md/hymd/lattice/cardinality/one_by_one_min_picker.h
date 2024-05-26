#pragma once

#include <unordered_set>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "algorithms/md/hymd/lattice/md_lattice_node_info.h"
#include "algorithms/md/hymd/lattice/validation_info.h"
#include "algorithms/md/hymd/md_lhs.h"
#include "model/index.h"

namespace algos::hymd::lattice::cardinality {

class OneByOnePicker {
private:
    enum class ComparisonResult {
        Specialization,
        Generalization,
        Incomparable,
    };

    std::vector<ValidationInfo> currently_picked_;

    static bool IsGeneralization(MdLhs::iterator gen_it, MdLhs::iterator spec_it,
                                 MdLhs::iterator gen_end, MdLhs::iterator spec_end);

    static bool IsGeneralizationIncr(MdLhs::iterator gen_it, MdLhs::iterator spec_it,
                                     MdLhs::iterator gen_end, MdLhs::iterator spec_end) {
        // Fewer classifiers in gen, gen generalizes spec.
        if (++gen_it == gen_end) return true;
        // More classifiers in gen, gen is incomparable with spec.
        if (++spec_it == spec_end) return false;
        return IsGeneralization(gen_it, spec_it, gen_end, spec_end);
    }

    static ComparisonResult CompareLhss(MdLhs const& cur, MdLhs const& prev);

public:
    static constexpr bool kNeedsEmptyRemoval = true;

    void NewBatch(std::size_t elements);
    void AddGeneralizations(MdLattice::MdVerificationMessenger& messenger,
                            boost::dynamic_bitset<>& considered_indices);
    std::vector<ValidationInfo> GetAll() noexcept;
};

}  // namespace algos::hymd::lattice::cardinality
