#include "pfdtane.h"

#include <chrono>
#include <iomanip>
#include <list>
#include <memory>

#include <easylogging++.h>

#include "config/error/option.h"
#include "config/error_measure/option.h"
#include "config/error_measure/type.h"
#include "config/max_lhs/option.h"
#include "config/names_and_descriptions.h"
#include "config/option.h"
#include "enums.h"
#include "fd/tane/lattice_level.h"
#include "fd/tane/lattice_vertex.h"
#include "model/table/column_data.h"
#include "model/table/column_layout_relation_data.h"
#include "model/table/relational_schema.h"

namespace algos {
using boost::dynamic_bitset;
using Cluster = model::PositionListIndex::Cluster;

PFDTane::PFDTane() : TaneCommon() {
    RegisterOptions();
}

config::ErrorType PFDTane::CalculateZeroAryFdError(ColumnData const* rhs,
                                                   ColumnLayoutRelationData const*) {
    std::size_t max = 1;
    model::PositionListIndex const* x_pli = rhs->GetPositionListIndex();
    for (Cluster const& x_cluster : x_pli->GetIndex()) {
        max = std::max(max, x_cluster.size());
    }
    return 1.0 - static_cast<double>(max) / x_pli->getRelationSize();
}

config::ErrorType PFDTane::CalculateFdError(model::PositionListIndex const* x_pli,
                                            model::PositionListIndex const* xa_pli,
                                            ErrorMeasure measure, ColumnLayoutRelationData const*) {
    std::deque<Cluster> xa_index = xa_pli->GetIndex();
    std::shared_ptr<Cluster const> probing_table = x_pli->CalculateAndGetProbingTable();
    std::sort(xa_index.begin(), xa_index.end(),
              [&probing_table](Cluster const& a, Cluster const& b) {
                  return probing_table->at(a.front()) < probing_table->at(b.front());
              });
    double sum = 0.0;
    std::size_t cluster_rows_count = 0;
    std::deque<Cluster> const& x_index = x_pli->GetIndex();
    auto xa_cluster_it = xa_index.begin();

    for (Cluster const& x_cluster : x_index) {
        std::size_t max = 1;
        for (int x_row : x_cluster) {
            if (xa_cluster_it == xa_index.end()) {
                break;
            }
            if (x_row == xa_cluster_it->at(0)) {
                max = std::max(max, xa_cluster_it->size());
                xa_cluster_it++;
            }
        }
        sum += measure == +ErrorMeasure::per_tuple ? static_cast<double>(max)
                                                   : static_cast<double>(max) / x_cluster.size();
        cluster_rows_count += x_cluster.size();
    }
    unsigned int unique_rows =
            static_cast<unsigned int>(x_pli->getRelationSize() - cluster_rows_count);
    double probability = static_cast<double>(sum + unique_rows) /
                         (measure == +ErrorMeasure::per_tuple ? x_pli->getRelationSize()
                                                              : x_index.size() + unique_rows);
    return 1.0 - probability;
}

config::ErrorType PFDTane::CalculateZeroAryFdError(ColumnData const* rhs) {
    return CalculateZeroAryFdError(rhs, relation_.get());
}

config::ErrorType PFDTane::CalculateFdError(model::PositionListIndex const* lhs_pli,
                                            model::PositionListIndex const* joint_pli) {
    return CalculateFdError(lhs_pli, joint_pli, error_measure_, relation_.get());
}

void PFDTane::RegisterOptions() {
    RegisterOption(config::ErrorOpt(&max_ucc_error_));
    RegisterOption(config::ErrorMeasureOpt(&error_measure_));
    RegisterOption(config::MaxLhsOpt(&max_lhs_));
}

void PFDTane::MakeExecuteOptsAvailable() {
    MakeOptionsAvailable({config::ErrorOpt.GetName(), config::ErrorMeasureOpt.GetName(),
                          config::MaxLhsOpt.GetName()});
}

}  // namespace algos
