#pragma once

#include <string>

#include "algorithms/fd/pli_based_fd_algorithm.h"
#include "config/error/type.h"
#include "config/max_lhs/type.h"
#include "model/table/position_list_index.h"
#include "model/table/relation_data.h"
#include "tane_common.h"

namespace algos {

class Tane : public TaneCommon {
private:
    void RegisterOptions();
    void MakeExecuteOptsAvailable() final;
    config::ErrorType CalculateZeroAryFdError(ColumnData const* rhs) override;
    config::ErrorType CalculateFdError(model::PositionListIndex const* lhs_pli,
                                       model::PositionListIndex const* joint_pli) override;

public:
    Tane();
    static config::ErrorType CalculateZeroAryFdError(ColumnData const* rhs,
                                                     ColumnLayoutRelationData const* relation_data);
    static config::ErrorType CalculateFdError(model::PositionListIndex const* lhs_pli,
                                              model::PositionListIndex const* joint_pli,
                                              ColumnLayoutRelationData const* relation_data);
};

}  // namespace algos
