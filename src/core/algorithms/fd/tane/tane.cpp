#include "tane.h"

#include <chrono>
#include <iomanip>
#include <list>
#include <memory>

#include <easylogging++.h>

#include "config/error/option.h"
#include "config/max_lhs/option.h"
#include "lattice_level.h"
#include "lattice_vertex.h"
#include "model/table/column_data.h"
#include "model/table/column_layout_relation_data.h"
#include "model/table/relational_schema.h"

namespace algos {

using boost::dynamic_bitset;

Tane::Tane() : TaneCommon() {
    RegisterOptions();
}

void Tane::RegisterOptions() {
    RegisterOption(config::ErrorOpt(&max_ucc_error_));
    RegisterOption(config::MaxLhsOpt(&max_lhs_));
}

void Tane::MakeExecuteOptsAvailable() {
    MakeOptionsAvailable({config::MaxLhsOpt.GetName(), config::ErrorOpt.GetName()});
}

config::ErrorType Tane::CalculateZeroAryFdError(ColumnData const* rhs,
                                                ColumnLayoutRelationData const* relation_data) {
    return 1 - rhs->GetPositionListIndex()->GetNepAsLong() /
                       static_cast<double>(relation_data->GetNumTuplePairs());
}

config::ErrorType Tane::CalculateFdError(model::PositionListIndex const* lhs_pli,
                                         model::PositionListIndex const* joint_pli,
                                         ColumnLayoutRelationData const* relation_data) {
    return (double)(lhs_pli->GetNepAsLong() - joint_pli->GetNepAsLong()) /
           static_cast<double>(relation_data->GetNumTuplePairs());
}

config::ErrorType Tane::CalculateZeroAryFdError(ColumnData const* rhs) {
    return CalculateZeroAryFdError(rhs, relation_.get());
}

config::ErrorType Tane::CalculateFdError(model::PositionListIndex const* lhs_pli,
                                         model::PositionListIndex const* joint_pli) {
    return CalculateFdError(lhs_pli, joint_pli, relation_.get());
}

}  // namespace algos
