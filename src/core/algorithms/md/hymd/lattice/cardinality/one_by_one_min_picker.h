#pragma once

#include <unordered_set>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "algorithms/md/hymd/lattice/md_lattice_node_info.h"
#include "algorithms/md/hymd/lattice/validation_info.h"
#include "model/index.h"

namespace algos::hymd::lattice::cardinality {

class OneByOnePicker {
private:
    std::vector<ValidationInfo> currently_picked_;

public:
    static constexpr bool kNeedsEmptyRemoval = true;

    void NewBatch(std::size_t elements);
    void AddGeneralizations(MdLattice::MdVerificationMessenger& messenger,
                            boost::dynamic_bitset<>& considered_indices);
    std::vector<ValidationInfo> GetAll() noexcept;
};

}  // namespace algos::hymd::lattice::cardinality
