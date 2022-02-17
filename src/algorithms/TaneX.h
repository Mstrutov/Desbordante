#pragma once

#include <string>

#include "CSVParser.h"
#include "PliBasedFDAlgorithm.h"
#include "PositionListIndex.h"
#include "RelationData.h"

class Tane : public PliBasedFDAlgorithm {
private:
    /* Special config parameters */
    constexpr static const char* kMaxError = "error";

    unsigned long long ExecuteInternal() override;
public:
    //TODO: these consts should go in class (or struct) Configuration
    const double max_fd_error_ = 0.01;
    const double max_ucc_error_ = 0.01;
    const unsigned int max_lhs_ = -1;

    int count_of_fd_ = 0;
    int count_of_ucc_ = 0;
    long apriori_millis_ = 0;

    explicit Tane(Config const& config)
        : PliBasedFDAlgorithm(config, {kDefaultPhaseName}), max_fd_error_(GetSpecialParam<double>(kMaxError)),
          max_ucc_error_(GetSpecialParam<double>(kMaxError)),
          max_lhs_(config_.max_lhs) {}
    explicit Tane(std::shared_ptr<ColumnLayoutRelationData> relation, Config const& config)
        : PliBasedFDAlgorithm(std::move(relation), config, {kDefaultPhaseName}),
          max_fd_error_(GetSpecialParam<double>(kMaxError)),
          max_ucc_error_(GetSpecialParam<double>(kMaxError)),
          max_lhs_(config_.max_lhs) {}

    static double CalculateZeroAryFdError(ColumnData const* rhs,
                                          ColumnLayoutRelationData const* relation_data);
    static double CalculateFdError(util::PositionListIndex const* lhs_pli,
                                   util::PositionListIndex const* joint_pli,
                                   ColumnLayoutRelationData const* relation_data);
    static double CalculateUccError(util::PositionListIndex const* pli,
                                    ColumnLayoutRelationData const* relation_data);

    //static double round(double error) { return ((int)(error * 32768) + 1)/ 32768.0; }

    void RegisterFd(Vertical const& lhs, Column const* rhs,
                    double error, RelationalSchema const* schema);
    // void RegisterFd(Vertical const* lhs, Column const* rhs, double error, RelationalSchema const* schema);
    void RegisterUcc(Vertical const& key, double error, RelationalSchema const* schema);
};
