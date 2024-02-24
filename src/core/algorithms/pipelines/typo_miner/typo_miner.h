#pragma once

#include "algorithms/algorithm.h"
#include "algorithms/create_algorithm.h"
#include "algorithms/fd/pyro/pyro.h"
#include "config/equal_nulls/type.h"
#include "config/names.h"
#include "config/tabular_data/input_table_type.h"
#include "model/table/column_layout_typed_relation_data.h"
#include "model/table/idataset_stream.h"
#include "parser/csv_parser/csv_parser.h"
#include "types.h"

namespace algos {

class TypoMiner : public Algorithm {
private:
    std::unique_ptr<FDAlgorithm> precise_algo_;
    std::unique_ptr<FDAlgorithm> approx_algo_;
    std::vector<FD> approx_fds_;
    std::unique_ptr<model::ColumnLayoutTypedRelationData> typed_relation_;
    /* Config members */
    double radius_; /* Maximal distance between two values to consider one of them a typo */
    double ratio_;  /* Maximal fraction of deviations per cluster to flag the cluster as
                     * containing typos */

    config::InputTable input_table_;
    config::EqNullsType is_null_equal_null_;
    std::shared_ptr<ColumnLayoutRelationData> relation_;

    PliBasedFDAlgorithm::ColumnLayoutRelationDataManager const relation_manager_;

    PliBasedFDAlgorithm::ColumnLayoutRelationDataManager MakeRelationManager() {
        return {&input_table_, &is_null_equal_null_, &relation_};
    }

    void ResetState() final;

    void LoadDataInternal() final;
    unsigned long long ExecuteInternal() final;

    static bool FDLess(FD const& l, FD const& r);
    static auto MakeTuplesByIndicesComparator(std::map<int, unsigned> const& frequency_map);

    std::map<int, unsigned> CreateFrequencyMap(Column const& cluster_col,
                                               model::PLI::Cluster const& cluster) const;
    unsigned GetMostFrequentValueIndex(Column const& cluster_col,
                                       model::PLI::Cluster const& cluster) const;

    bool ValuesAreClose(std::byte const* l, std::byte const* r, model::Type const& type) const {
        assert(type.IsMetrizable());
        return static_cast<model::IMetrizableType const&>(type).Dist(l, r) < radius_;
    }

    template <typename GetPrecise, typename GetApprox>
    TypoMiner(GetPrecise get_precise, GetApprox get_approx);

    void RegisterOptions();
    void MakeExecuteOptsAvailable() final;
    void AddSpecificNeededOptions(
            std::unordered_set<std::string_view>& previous_options) const final;
    bool SetExternalOption(std::string_view option_name, boost::any const& value) final;
    int TrySetOption(std::string_view option_name, boost::any const& value_precise,
                     boost::any const& value_approx);

public:
    using TyposVec = std::vector<model::PLI::Cluster::value_type>;
    using ClusterTyposPair = std::pair<model::PLI::Cluster, TyposVec>;

    struct SquashedElement {
        int tuple_index; /* Tuple index */
        unsigned amount; /* The number of tuples equal to the given and
                          * following immediately after the given */
    };

    explicit TypoMiner(AlgorithmType precise, AlgorithmType approx = AlgorithmType::pyro);

    std::vector<model::PLI::Cluster> FindClustersWithTypos(FD const& typos_fd,
                                                           bool const sort_clusters = true) const;
    /* Returns squashed representation of a cluster with respect to given fd.
     * Check description of SquashedElement. Two tuples are considered equal if they are
     * equal in rhs attribute of given fd (they are automatically equal in lhs too if clusters
     * were retrieved from FindClustersWithTypos()).
     */
    std::vector<SquashedElement> SquashCluster(FD const& squash_on,
                                               model::PLI::Cluster const& cluster) const;
    /* Sorts given cluster in ascending order by uniqueness of the values in rhs of sort_on fd.
     * More strictly:
     * t1 tuple is considered less than t2 tuple iff t1[sort_on.rhs] is less frequent value
     * in rhs column than t2[sort_on.rhs].
     */
    void SortCluster(FD const& sort_on, model::PLI::Cluster& cluster) const;

    /* Sorts given cluster in ascending order of the indices */
    void RestoreLineOrder(model::PLI::Cluster& cluster) const;
    /* Sorts given squashed cluster in ascending order of the indices */
    void RestoreLineOrder(std::vector<TypoMiner::SquashedElement>& squashed_cluster) const;

    /* Finds lines in a cluster that has typos on typos_fd.GetRhs() column values. A value is said
     * to contain a typo if it differs from the most frequent value in the cluster by less
     * than radius_ and the fraction of such values (values_num / cluster_size) is less than ratio_.
     * NOTE: The cluster argument is assumed to be consistent with the rhs of typos_fd, i.e.
     * indices in the cluster represent value indices in the typos_fd.GetRhs() column. Otherwise,
     * the behavior of this method is undefined.
     * Most likely you want to pass to this method as arguments pure approximate FD some_fd and
     * a cluster retrieved from the FindClustersWithTypos(some_fd).
     */
    TyposVec FindLinesWithTypos(FD const& typos_fd, model::PLI::Cluster const& cluster) const;
    TyposVec FindLinesWithTypos(FD const& typos_fd, model::PLI::Cluster const& cluster,
                                double new_radius, double new_ratio);
    std::vector<ClusterTyposPair> FindClustersAndLinesWithTypos(
            FD const& typos_fd, bool const sort_clusters = true) const;

    /* Returns vector of approximate fds only (there are no precise fds) */
    std::vector<FD> const& GetApproxFDs() const noexcept {
        return approx_fds_;
    }

    double GetRadius() const noexcept {
        return radius_;
    }

    double GetRatio() const noexcept {
        return ratio_;
    }

    double SetRadius(double radius) {
        SetOption(config::names::kRadius, radius);
        return radius_;
    }

    double SetRatio(double ratio) {
        SetOption(config::names::kRatio, ratio);
        return ratio_;
    }

    ColumnLayoutRelationData const& GetRelationData() const noexcept {
        assert(relation_ != nullptr);
        return *relation_;
    }

    std::string GetApproxFDsAsJson() const {
        return FDAlgorithm::FDsToJson(approx_fds_);
    }
};

}  // namespace algos
