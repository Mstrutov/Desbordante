#pragma once

#include <memory>
#include <unordered_map>
#include <unordered_set>

#include <boost/dynamic_bitset.hpp>

#include "config/tabular_data/input_table_type.h"
#include "fd/fd_algorithm.h"
#include "model/table/column.h"
#include "model/table/relational_schema.h"
#include "model/table/vertical.h"
#include "search_tree.h"

namespace algos {

class Aid : public FDAlgorithm {
private:
    using Cluster = std::vector<size_t>;

    config::InputTable input_table_;

    std::unique_ptr<RelationalSchema> schema_{};
    std::vector<std::vector<size_t>> tuples_;

    size_t number_of_attributes_{};
    size_t number_of_tuples_{};

    std::unordered_set<boost::dynamic_bitset<>> neg_cover_{};

    constexpr static double const growth_threshold_ = 0.01;
    constexpr static size_t const window_size_ = 10;
    constexpr static size_t const prime_ = 10619863;

    std::vector<double> prev_ratios_;
    double sum_{};

    std::vector<std::unordered_map<size_t, Cluster>> clusters_;
    std::vector<std::vector<size_t>> indices_in_clusters_;

    boost::dynamic_bitset<> constant_columns_;

    void RegisterOptions();

    void ResetStateFd() final;

    void LoadDataInternal() final;
    unsigned long long ExecuteInternal() final;

    void BuildClusters();
    void CreateNegativeCover();
    void InvertNegativeCover();

    void HandleTuple(size_t tuple_num, size_t index);
    void HandleInvalidFd(boost::dynamic_bitset<> const& neg_cover_el, SearchTree& pos_cover_tree,
                         size_t rhs);
    size_t GenerateSecondClusterIndex(size_t index_in_cluster, size_t iteration_num) const;
    bool IsNegativeCoverGrowthSmall(size_t iteration_num, double curr_ratio);

    void HandleConstantColumns(boost::dynamic_bitset<>& attributes);
    void RegisterFDs(size_t rhs, std::vector<boost::dynamic_bitset<>> const& list_of_lhs);

    static boost::dynamic_bitset<> ChangeAttributesOrder(
            boost::dynamic_bitset<> const& initial_bitset, std::vector<size_t> const& new_order);
    std::vector<size_t> GetAttributesSortedByFrequency(
            std::vector<boost::dynamic_bitset<>> const& neg_cover_vector) const;

    boost::dynamic_bitset<> BuildAgreeSet(size_t t1, size_t t2);

public:
    Aid();
};

}  // namespace algos
