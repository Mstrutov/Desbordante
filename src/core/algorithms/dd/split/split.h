#pragma once

#include <cstddef>
#include <filesystem>
#include <list>
#include <memory>
#include <utility>
#include <vector>

#include "algorithms/algorithm.h"
#include "algorithms/dd/dd.h"
#include "config/tabular_data/input_table_type.h"
#include "enums.h"
#include "model/table/column_layout_relation_data.h"
#include "model/table/column_layout_typed_relation_data.h"

namespace algos::dd {

using df = model::DF;
using dd = model::DD;

class Split : public Algorithm {
private:
    config::InputTable input_table_;

    std::shared_ptr<ColumnLayoutRelationData> relation_;
    std::shared_ptr<model::ColumnLayoutTypedRelationData> typed_relation_;
    unsigned num_rows_;
    unsigned num_columns_;

    bool has_dif_table_;

    config::InputTable difference_table_;
    std::unique_ptr<model::ColumnLayoutTypedRelationData> difference_typed_relation_;

    Reduce reduce_method_ = Reduce::Hybrid;  // currently, the fastest method

    std::vector<std::pair<double, double>> min_max_dif_;
    std::vector<std::vector<std::vector<double>>> distances_;
    std::list<std::pair<std::size_t, std::size_t>> tuple_pairs_;
    std::list<dd> results_;

    void RegisterOptions();
    void SetLimits();
    void ParseDifferenceTable();

    void ResetState() final {
        results_.clear();
    }

    double CalculateDistance(std::size_t column_index,
                             std::pair<std::size_t, std::size_t> tuple_pair);
    bool CheckDF(df const& dep, std::pair<std::size_t, std::size_t> tuple_pair);
    bool VerifyDD(dd const& dep);
    void CalculateAllDistances();
    bool IsFeasible(df const& d);
    std::list<df> SearchSpace(std::list<size_t>& indices);
    std::list<df> SearchSpace(std::size_t index);
    bool Subsume(df const& df1, df const& df2);
    std::list<dd> NegativePruningReduce(df const& rhs, std::list<df> const& search, unsigned& cnt);
    std::list<dd> HybridPruningReduce(df const& rhs, std::list<df> const& search, unsigned& cnt);
    std::list<dd> InstanceExclusionReduce(
            std::list<std::pair<std::size_t, std::size_t>> const& tuple_pairs,
            std::list<df> const& search, df const& rhs, unsigned& cnt);
    void CalculateTuplePairs();
    unsigned ReduceDDs(auto start_time);
    unsigned RemoveRedundantDDs();
    unsigned RemoveTransitiveDDs();
    model::DDString DDToDDString(dd const& dd);
    void PrintResults();

protected:
    void LoadDataInternal() override;
    void MakeExecuteOptsAvailable() override;
    unsigned long long ExecuteInternal() override;

public:
    Split();
    std::list<dd> GetResults();
    std::vector<std::pair<double, double>> GetMinMaxDif();
    std::list<model::DDString> GetDDStringList();
};

}  // namespace algos::dd
