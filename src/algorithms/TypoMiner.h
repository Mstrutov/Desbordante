#pragma once

#include "Primitive.h"
#include "Pyro.h"

namespace algos {

template<typename PreciseAlgo, typename ApproxAlgo = Pyro>
class TypoMiner : public Primitive {
private:
    std::unique_ptr<FDAlgorithm> precise_algo_;
    std::unique_ptr<FDAlgorithm> approx_algo_;
    std::vector<FD> approx_fds_;
    std::shared_ptr<ColumnLayoutRelationData> relation_;

    static bool FDLess(FD const& l, FD const& r);
    static auto MakeTuplesByIndicesComparator(std::map<int, unsigned> const& frequency_map);

    std::map<int, unsigned> CreateFrequencyMap(Column const& cluster_col,
                                               util::PLI::Cluster const& cluster) const;

    void ClarifyTypos(FD const& fd) {
        /* TODO(polyntsov): implement this */
        (void)fd;
    }

public:
    using Config = FDAlgorithm::Config;

    struct SquashedElement {
        int tuple_index;  /* Tuple index */
        unsigned amount;  /* The number of tuples equal to the given and
                           * following immediately after the given */
    };

    TypoMiner(Config const& config);

    unsigned long long Execute() override;

    std::vector<util::PLI::Cluster> FindClustersWithTypos(FD const& typos_fd,
                                                          bool const sort_clusters = true);
    /* Returns squashed representation of a cluster with respect to given fd.
     * Check description of SquashedElement. Two tuples are considered equal if they are
     * equal in rhs attribute of given fd (they are automatically equal in lhs too if clusters
     * were retrieved from FindClustersWithTypos()).
     */
    std::vector<SquashedElement> SquashCluster(FD const& squash_on,
                                               util::PLI::Cluster const& cluster);
    /* Sorts given cluster in ascending order by uniqueness of the values in rhs of sort_on fd.
     * More strictly:
     * t1 tuple is considered less than t2 tuple iff t1[sort_on.rhs] is less frequent value
     * in rhs column than t2[sort_on.rhs].
     */
    void SortCluster(FD const& sort_on, util::PLI::Cluster& cluster) const;

    /* Returns vector of approximate fds only (there are no precise fds) */
    std::vector<FD> const& GetApproxFDs() const noexcept {
        return approx_fds_;
    }
    std::string GetApproxFDsAsJson() const {
        return FDAlgorithm::FDsToJson(approx_fds_);
    }
};

template <typename PreciseAlgo, typename ApproxAlgo>
TypoMiner<PreciseAlgo, ApproxAlgo>::TypoMiner(Config const& config)
    : Primitive(config.data, config.separator, config.has_header,
                {"Precise fd algorithm execution", "Approximate fd algoritm execution",
                 "Extracting fds with non-zero error"}) {
    static_assert(std::is_base_of_v<PliBasedFDAlgorithm, ApproxAlgo>,
                  "Approximate algorithm must be relation based");
    if (config.GetSpecialParam<double>("error") == 0.0) {
        throw std::invalid_argument("Typos mining with error = 0 is meaningless");
    }
    Config precise_config = config;
    precise_config.special_params["error"] = 0.0;
    relation_ = ColumnLayoutRelationData::CreateFrom(input_generator_, config.is_null_equal_null);
    if constexpr (std::is_base_of_v<PliBasedFDAlgorithm, PreciseAlgo>) {
        precise_algo_ = std::make_unique<PreciseAlgo>(relation_, precise_config);
    } else {
        precise_algo_ = std::make_unique<PreciseAlgo>(precise_config);
    }
    approx_algo_ = std::make_unique<ApproxAlgo>(relation_, config);
}

template<typename PreciseAlgo, typename ApproxAlgo>
unsigned long long TypoMiner<PreciseAlgo, ApproxAlgo>::Execute() {
    auto const start_time = std::chrono::system_clock::now();

    precise_algo_->Execute();
    approx_algo_->Execute();

    std::list<FD>& precise_fds = precise_algo_->FdList();
    std::list<FD>& approx_fds = approx_algo_->FdList();

    precise_fds.sort(FDLess);
    approx_fds.sort(FDLess);

    std::set_difference(approx_fds.begin(), approx_fds.end(),
                        precise_fds.begin(),
                        precise_fds.end(),
                        std::back_inserter(approx_fds_), FDLess);

    auto const elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now() - start_time);

    return elapsed_milliseconds.count();
}

template <typename PreciseAlgo, typename ApproxAlgo>
std::vector<util::PLI::Cluster> TypoMiner<PreciseAlgo, ApproxAlgo>::FindClustersWithTypos(
    FD const& typos_fd, bool const sort_clusters) {
    std::vector<util::PLI::Cluster> clusters;
    std::shared_ptr<util::PLI const> intersection_pli;
    std::vector<Column const*> const lhs_columns = typos_fd.GetLhs().GetColumns();
    std::vector<int> const& probing_table =
        relation_->GetColumnData(typos_fd.GetRhs().GetIndex()).GetProbingTable();

    assert(lhs_columns.size() != 0);

    for (Column const* col : lhs_columns) {
        ColumnData const& col_data = relation_->GetColumnData(col->GetIndex());
        std::shared_ptr<util::PLI const> pli = col_data.GetPliOwnership();

        if (intersection_pli == nullptr) {
            intersection_pli = pli;
        } else {
            intersection_pli = intersection_pli->Intersect(pli.get());
        }
    }

    for (util::PLI::Cluster const& cluster : intersection_pli->GetIndex()) {
        int cluster_rhs_value = -1;

        /* Check if fd has wrong rhs values in this cluster */
        for (int const tuple_index : cluster) {
            int const probing_table_value = probing_table[tuple_index];

            if (cluster_rhs_value == -1) {
                cluster_rhs_value = probing_table_value;
            } else if (cluster_rhs_value != probing_table_value) {
                cluster_rhs_value = -1;
                break;
            }
        }

        if (cluster_rhs_value == -1 || (cluster_rhs_value == 0 && cluster.size() != 1)) {
            /* Actually intersection_pli is owned only by this method most of the time
             * (when lhs_columns.size() != 1), so we can here move cluster when
             * lhs_columns.size() != 1. But that leads to more cumbersome and complex code.
             * So I decided to leave it as it is until we know for sure that this place causes
             * performance problems.
             */
            clusters.push_back(cluster);

            if (sort_clusters) {
                std::map<int, unsigned> const frequency_map =
                    CreateFrequencyMap(typos_fd.GetRhs(), clusters.back());
                std::stable_sort(clusters.back().begin(), clusters.back().end(),
                                 MakeTuplesByIndicesComparator(frequency_map));
            }
        }
    }

    return clusters;
}

template <typename PreciseAlgo, typename ApproxAlgo>
std::vector<typename TypoMiner<PreciseAlgo, ApproxAlgo>::SquashedElement>
TypoMiner<PreciseAlgo, ApproxAlgo>::SquashCluster(FD const& squash_on,
                                                   util::PLI::Cluster const& cluster) {
    std::vector<SquashedElement> squashed;
    std::vector<int> const& probing_table =
        relation_->GetColumnData(squash_on.GetRhs().GetIndex()).GetProbingTable();

    if(cluster.empty()) {
        return squashed;
    }

    auto prev = cluster.begin();
    squashed.push_back({ .tuple_index = *prev, .amount = 1});

    for (auto it = std::next(cluster.cbegin()); it != cluster.cend(); ++it) {
        if (probing_table[*it] == probing_table[*prev]) {
            squashed.back().amount++;
        } else {
            squashed.push_back({ .tuple_index = *it, .amount = 1 });
        }
        prev = it;
    }

    return squashed;
}

template <typename PreciseAlgo, typename ApproxAlgo>
void TypoMiner<PreciseAlgo, ApproxAlgo>::SortCluster(FD const& sort_on,
                                                      util::PLI::Cluster& cluster) const {
    std::map<int, unsigned> const frequency_map = CreateFrequencyMap(sort_on.GetRhs(), cluster);

    std::stable_sort(cluster.begin(), cluster.end(), MakeTuplesByIndicesComparator(frequency_map));
}

template <typename PreciseAlgo, typename ApproxAlgo>
std::map<int, unsigned> TypoMiner<PreciseAlgo, ApproxAlgo>::CreateFrequencyMap(
    Column const& cluster_col, util::PLI::Cluster const& cluster) const {
    std::map<int, unsigned> frequencies;
    std::map<int, unsigned> frequency_map;
    std::vector<int> const& probing_table =
        relation_->GetColumnData(cluster_col.GetIndex()).GetProbingTable();

    for (int const tuple_index : cluster) {
        int const probing_table_value = probing_table[tuple_index];

        if (probing_table_value != 0) {
            frequencies[probing_table_value]++;
        }
    }

    for (int const tuple_index : cluster) {
        int const probing_table_value = probing_table[tuple_index];
        int const value = (probing_table_value == 0) ? 1 : frequencies[probing_table_value];
        frequency_map[tuple_index] = value;
    }

    return frequency_map;
}

template <typename PreciseAlgo, typename ApproxAlgo>
bool TypoMiner<PreciseAlgo, ApproxAlgo>::FDLess(FD const& l, FD const& r) {
    if (l.GetLhs() < r.GetLhs()) {
        return true;
    } else if (l.GetLhs() == r.GetLhs()) {
        return l.GetRhs() < r.GetRhs();
    }

    return false;
}

template <typename PreciseAlgo, typename ApproxAlgo>
auto TypoMiner<PreciseAlgo, ApproxAlgo>::MakeTuplesByIndicesComparator(
    std::map<int, unsigned> const& frequency_map) {
    return [&frequency_map](int const l, int const r) {
        return frequency_map.at(l) < frequency_map.at(r);
    };
}

}  // namespace algos

