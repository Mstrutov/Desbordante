#pragma once

#include "algorithms/fd/pli_based_fd_algorithm.h"
#include "algorithms/fd/tane/lattice_level.h"
#include "algorithms/fd/tane/tane_common.h"
#include "config/error/type.h"
#include "config/error_measure/type.h"
#include "config/max_lhs/type.h"
#include "enums.h"
#include "model/table/position_list_index.h"
#include "model/table/relation_data.h"

namespace algos {

class PFDTane : public TaneCommon {
private:
    ErrorMeasure error_measure_ = +ErrorMeasure::per_tuple;
    void RegisterOptions();
    void MakeExecuteOptsAvailable() final;
    config::ErrorType CalculateZeroAryFdError(ColumnData const* rhs) override;
    config::ErrorType CalculateFdError(model::PositionListIndex const* lhs_pli,
                                       model::PositionListIndex const* joint_pli) override;

public:
    PFDTane();
    static config::ErrorType CalculateZeroAryFdError(ColumnData const* rhs,
                                                     ColumnLayoutRelationData const* relation_data);
    static config::ErrorType CalculateFdError(model::PositionListIndex const* x_pli,
                                              model::PositionListIndex const* xa_pli,
                                              ErrorMeasure error_measure,
                                              ColumnLayoutRelationData const* relation_data);
};

}  // namespace algos
