#pragma once

#include <array>
#include <memory>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <boost/dynamic_bitset.hpp>
#include <boost/unordered_set.hpp>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>

#include "config/tabular_data/input_table/option.h"
#include "custom_random/type.h"
#include "config/equal_nulls/option.h"
#include "custom_random/option.h"
#include "util/custom_random.h"

#include "fd/fd_algorithm.h"
#include "model/table/column.h"
#include "model/table/relational_schema.h"
#include "model/table/vertical.h"

#include "search_tree.h"
#include "mlfq.h"

namespace algos {

class EulerFD : public FDAlgorithm {
    using Bitset = boost::dynamic_bitset<>;
    using RandomStrategy = Cluster::RandomStrategy;

    // random strategy for unit tests
    config::CustomRandomFlagType custom_random_opt_;
    RandomStrategy rand_function_;
    std::unique_ptr<CustomRandom> random_ {};
    constexpr static std::size_t random_upper_bound_ = 3047102;

    // data from load data
    size_t number_of_attributes_ {};
    size_t number_of_tuples_ {};
    config::InputTable input_table_;
    std::unique_ptr<RelationalSchema> schema_ {};
    std::vector<std::vector<size_t>> tuples_;

    config::EqNullsType is_null_equal_null_ {};

    // thresholds to checking criterion of euler fd cycles
    constexpr static double pos_cover_growth_treshold_ = 0.01;
    constexpr static double neg_cover_growth_treshold_ = 0.01;
    constexpr static size_t window_ = 3;
    std::array<double, window_> last_ncover_ratios_;
    std::array<double, window_> last_pcover_ratios_;

    // clusters from partition
    std::vector<Cluster> clusters_;
    Bitset constant_columns_;

    // mlfq sampling data
    bool is_first_sample_ = true;
    constexpr static const size_t queues_number_ = 5;
    MLFQ mlfq_;
    constexpr static double initial_effective_treshold_ = 0.01;
    double effective_treshold_ = initial_effective_treshold_;

    // invalid fds storages
    std::unordered_set<Bitset> invalids_;
    std::unordered_set<Bitset> new_invalids_;
    size_t fd_num_ = 0;
    size_t old_invalid_size_ = 0;

    // covers of fds
    std::vector<SearchTreeEulerFD> negative_cover_;
    std::vector<SearchTreeEulerFD> positive_cover_;

    // data for sorting columns by number of 1
    std::vector<size_t> attribute_indexes_;
    std::vector<size_t> attribute_frequences_;

    void LoadDataInternal() final;
    unsigned long long ExecuteInternal() final;
    void ResetStateFd() final;
    void MakeExecuteOptsAvailable() final;

    [[nodiscard]] Bitset BuildAgreeSet(size_t t1, size_t t2) const;

    void InitCovers();
    void BuildPartition();

    double SamlingInCluster(Cluster *cluster);
    void Sampling();
    size_t GenerateResults();

    [[nodiscard]] std::vector<size_t> GetAttributesSortedByFrequency(
        const std::vector<Bitset> &neg_cover_vector);
    [[nodiscard]] static Bitset ChangeAttributesOrder(
         const Bitset &initial_bitset, const std::vector<size_t> &new_order);

    [[nodiscard]] std::vector<Bitset> CreateNegativeCover(size_t rhs,
      const std::vector<Bitset> &neg_cover_vector);
    size_t Invert(size_t rhs, const std::vector<Bitset> &neg);

    static void AddInvalidAtTree(SearchTreeEulerFD &tree, const Bitset &invalid);
    static std::unordered_set<Bitset> RemoveGeneralizations(
        SearchTreeEulerFD &tree, const Bitset &invalid);

    bool IsNCoverGrowthSmall() const;
    bool IsPCoverGrowthSmall() const;

    void HandleInvalid();

    void SaveAnswer();

public:
    EulerFD();
};
}
