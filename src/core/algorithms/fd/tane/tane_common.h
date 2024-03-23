#pragma once

#include "algorithms/fd/pli_based_fd_algorithm.h"
#include "algorithms/fd/tane/lattice_level.h"
#include "config/error/type.h"
#include "config/max_lhs/type.h"
#include "model/table/position_list_index.h"
#include "model/table/relation_data.h"

namespace algos {

class TaneCommon : public PliBasedFDAlgorithm {
protected:
    config::ErrorType max_fd_error_;
    config::ErrorType max_ucc_error_;
    config::MaxLhsType max_lhs_;

private:
    void ResetStateFd() final;
    void Prune(model::LatticeLevel* level);
    void ComputeDependencies(model::LatticeLevel* level);
    unsigned long long ExecuteInternal() final;
    virtual config::ErrorType CalculateZeroAryFdError(ColumnData const* rhs) = 0;
    virtual config::ErrorType CalculateFdError(model::PositionListIndex const* lhs_pli,
                                               model::PositionListIndex const* joint_pli) = 0;

public:
    TaneCommon();
    static double CalculateUccError(model::PositionListIndex const* pli,
                                    ColumnLayoutRelationData const* relation_data);
    void RegisterAndCountFd(Vertical const& lhs, Column const* rhs, double error,
                            RelationalSchema const* schema);
};

}  // namespace algos
