#pragma once

#include "algorithms/fd/pli_based_fd_algorithm.h"
#include "algorithms/fd/tane/lattice_level.h"
#include "config/error/type.h"
#include "enums.h"
#include "model/table/column_data.h"
#include "model/table/column_layout_relation_data.h"
#include "model/table/position_list_index.h"

namespace algos {
namespace tane {

class TaneCommon : public PliBasedFDAlgorithm {
protected:
    config::ErrorType max_fd_error_;
    config::ErrorType max_ucc_error_;

private:
    void ResetStateFd() final;
    void Prune(model::LatticeLevel* level);
    void ComputeDependencies(model::LatticeLevel* level);
    unsigned long long ExecuteInternal() final;
    virtual config::ErrorType CalculateZeroAryFdError(ColumnData const* rhs) = 0;
    virtual config::ErrorType CalculateFdError(model::PositionListIndex const* lhs_pli,
                                               model::PositionListIndex const* joint_pli) = 0;

public:
    TaneCommon(std::optional<ColumnLayoutRelationDataManager> relation_manager = std::nullopt);
    static double CalculateUccError(model::PositionListIndex const* pli,
                                    ColumnLayoutRelationData const* relation_data);
    void RegisterAndCountFd(Vertical const& lhs, Column const* rhs, double error,
                            RelationalSchema const* schema);
};

}  // namespace tane

class Tane : public tane::TaneCommon {
private:
    void RegisterOptions();
    void MakeExecuteOptsAvailableFDInternal() override final;
    config::ErrorType CalculateZeroAryFdError(ColumnData const* rhs) override;
    config::ErrorType CalculateFdError(model::PositionListIndex const* lhs_pli,
                                       model::PositionListIndex const* joint_pli) override;

public:
    Tane(std::optional<ColumnLayoutRelationDataManager> relation_manager = std::nullopt);
};

class PFDTane : public tane::TaneCommon {
private:
    ErrorMeasure error_measure_ = +ErrorMeasure::per_tuple;
    void RegisterOptions();
    void MakeExecuteOptsAvailableFDInternal() final;
    config::ErrorType CalculateZeroAryFdError(ColumnData const* rhs) override;
    config::ErrorType CalculateFdError(model::PositionListIndex const* lhs_pli,
                                       model::PositionListIndex const* joint_pli) override;

public:
    PFDTane(std::optional<ColumnLayoutRelationDataManager> relation_manager = std::nullopt);
    static config::ErrorType CalculateZeroAryPFDError(
            ColumnData const* rhs, ColumnLayoutRelationData const* relation_data);
    static config::ErrorType CalculatePFDError(model::PositionListIndex const* x_pli,
                                               model::PositionListIndex const* xa_pli,
                                               ErrorMeasure error_measure,
                                               ColumnLayoutRelationData const* relation_data);
};

}  // namespace algos
