#pragma once

#include "algorithms/fd/pli_based_fd_algorithm.h"
#include "algorithms/fd/tane/lattice_level.h"
#include "config/error/type.h"
#include "config/error_measure/type.h"
#include "config/max_lhs/type.h"
#include "enums.h"
#include "model/table/position_list_index.h"
#include "model/table/relation_data.h"

namespace algos {

class PFDTane : public PliBasedFDAlgorithm {
private:
    config::ErrorType max_fd_error_;
    config::ErrorType max_ucc_error_;
    config::MaxLhsType max_lhs_;
    ErrorMeasure error_measure_ = +ErrorMeasure::per_tuple;

    void ResetStateFd() final;
    void RegisterOptions();
    void MakeExecuteOptsAvailable() final;
    void Prune(model::LatticeLevel* level);
    void ComputeDependencies(model::LatticeLevel* level);
    unsigned long long ExecuteInternal() final;

public:
    PFDTane();
    static double CalculateUccError(model::PositionListIndex const* pli,
                                    ColumnLayoutRelationData const* relation_data);

    void RegisterAndCountFd(Vertical const& lhs, Column const* rhs, double error,
                            RelationalSchema const* schema);
    static config::ErrorType CalculateZeroAryFdError(ColumnData const* rhs);
    static config::ErrorType CalculateFdError(model::PositionListIndex const* x_pli,
                                              model::PositionListIndex const* xa_pli,
                                              ErrorMeasure error_measure);
};

}  // namespace algos
