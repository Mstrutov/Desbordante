#pragma once

#include <cstddef>
#include <unordered_set>
#include <vector>

#include "algorithms/md/decision_boundary.h"
#include "algorithms/md/hymd/column_match_info.h"
#include "algorithms/md/hymd/indexes/records_info.h"
#include "algorithms/md/hymd/lattice/md_lattice.h"
#include "algorithms/md/hymd/lattice/validation_info.h"
#include "algorithms/md/hymd/recommendation.h"
#include "algorithms/md/hymd/table_identifiers.h"
#include "algorithms/md/hymd/utility/invalidated_rhss.h"
#include "model/index.h"

namespace algos::hymd {

class Validator {
public:
    struct Result {
        std::vector<std::vector<Recommendation>> recommendations;
        utility::InvalidatedRhss invalidated;
        bool is_unsupported;
    };

private:
    template <typename PairProvider>
    class SetPairProcessor;
    class OneCardPairProvider;
    class MultiCardPairProvider;

    indexes::RecordsInfo const* const records_info_;
    std::vector<ColumnMatchInfo> const* const column_matches_info_;
    std::size_t const min_support_;
    lattice::MdLattice* const lattice_;

    [[nodiscard]] bool Supported(std::size_t support) const noexcept {
        return support >= min_support_;
    }

    [[nodiscard]] model::Index GetLeftPliIndex(model::Index const column_match_index) const {
        return (*column_matches_info_)[column_match_index].left_column_index;
    }

    indexes::DictionaryCompressor const& GetLeftCompressor() const noexcept {
        return records_info_->GetLeftCompressor();
    }

    indexes::DictionaryCompressor const& GetRightCompressor() const noexcept {
        return records_info_->GetRightCompressor();
    }

    std::size_t GetLeftValueNum(model::Index const col_match_index) const {
        return GetLeftCompressor().GetPli(GetLeftPliIndex(col_match_index)).GetClusters().size();
    }

    std::size_t GetTotalPairsNum() const {
        return GetLeftCompressor().GetNumberOfRecords() * GetRightCompressor().GetNumberOfRecords();
    }

    [[nodiscard]] indexes::RecSet const* GetSimilarRecords(ValueIdentifier value_id,
                                                           model::md::DecisionBoundary lhs_bound,
                                                           model::Index column_match_index) const;

public:
    Validator(indexes::RecordsInfo const* records_info,
              std::vector<ColumnMatchInfo> const& column_matches_info, std::size_t min_support,
              lattice::MdLattice* lattice)
        : records_info_(records_info),
          column_matches_info_(&column_matches_info),
          min_support_(min_support),
          lattice_(lattice) {}

    [[nodiscard]] Result Validate(lattice::ValidationInfo& validation_info) const;
};

}  // namespace algos::hymd
