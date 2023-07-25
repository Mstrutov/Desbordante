#pragma once

#include <functional>
#include <set>

#include <boost/thread/mutex.hpp>

#include "algorithms/pli_based_fd_algorithm.h"
#include "model/column_layout_relation_data.h"
#include "model/vertical.h"
#include "util/config/max_lhs/type.h"
#include "util/config/thread_number/type.h"

namespace algos {

class FastFDs : public PliBasedFDAlgorithm {
public:
    FastFDs(bool request_prepared_data = false);

private:
    using OrderingComparator = std::function<bool(Column const&, Column const&)>;
    using DiffSet = Vertical;

    void RegisterOptions();

    void ResetStateFd() final;
    unsigned long long ExecuteInternal() final;

    // Computes all difference sets of `relation_` by complementing agree sets
    void GenDiffSets();

    /* Computes minimal difference sets
     * of `relation_` modulo `col`
     */
    std::vector<DiffSet> GetDiffSetsMod(Column const& col) const;
    /* Returns initial ordering,
     * the total ordering of { schema_->GetColumns() \ `attribute` } according to `diff_sets`
     */
    std::set<Column, OrderingComparator> GetInitOrdering(std::vector<DiffSet> const& diff_sets,
                                                         Column const& attribute) const;
    /* Returns next ordering,
     * the total ordering of { B in schema_->GetColumns() | B > `attribute` (in `cur_ordering`) }
     * according to `diff_sets`
     */
    std::set<Column, OrderingComparator> GetNextOrdering(
        std::vector<DiffSet> const& diff_sets, Column const& attribute,
        std::set<Column, OrderingComparator> const& cur_ordering) const;
    void FindCovers(Column const& attribute, std::vector<DiffSet> const& diff_sets_mod,
                    std::vector<DiffSet> const& cur_diff_sets, Vertical const& path,
                    std::set<Column, OrderingComparator> const& ordering);
    /* Returns true if `cover` is the minimal cover of `diff_sets_mod`,
     * false otherwise
     */
    bool CoverMinimal(Vertical const& cover, std::vector<DiffSet> const& diff_sets_mod) const;
    /* Returns true if `candidate` covers `sets`,
     * false otherwise
     */
    bool IsCover(Vertical const& candidate, std::vector<Vertical> const& sets) const;
    /* Returns true if `l_col` > `r_col`,
     * false otherwise.
     * `l_col` > `r_col` iff
     * `l_col` covers more sets in `diff_sets` than `r_col` or
     * `l_col` and `r_col` cover the same number of sets but
     * `l_col` index less than `r_col` index
     */
    bool OrderingComp(std::vector<DiffSet> const& diff_sets,
                      Column const& l_col, Column const& r_col) const;
    bool ColumnContainsOnlyEqualValues(Column const& column) const;

    RelationalSchema const* schema_;
    std::vector<DiffSet> diff_sets_;
    util::config::ThreadNumType threads_num_;
    util::config::MaxLhsType max_lhs_;
    double percent_per_col_;
};

}  // namespace algos
