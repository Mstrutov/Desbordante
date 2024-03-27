#pragma once

#include <string>
#include <utility>
#include <vector>

#include "algorithm.h"
#include "config/tabular_data/input_table_type.h"
#include "frequency_handler.h"
#include "index_pair.h"
#include "model/table/column_index.h"
#include "model/table/column_layout_typed_relation_data.h"
#include "sample.h"

namespace algos {

class SFDAlgorithm : public Algorithm {
private:
    using TypedRelation = model::ColumnLayoutTypedRelationData;
    config::InputTable input_table_;
    std::unique_ptr<TypedRelation> typed_relation_;
    bool only_sfd_;
    long double minimum_cardinality_;
    long double max_diff_vals_proportion_;
    long double min_sfd_strength_measure_;
    long double min_skew_threshold_;
    long double min_structural_zeroes_amount_;
    long double max_false_positive_probability_;
    long double delta_;
    size_t max_amount_of_categories_;
    std::vector<bool> is_skewed_;
    std::vector<size_t> domains_;
    std::vector<model::ColumnIndex> soft_keys_{};
    std::vector<model::ColumnIndex> trivial_columns_{};
    std::vector<ColumnIndexPair> correlations_;
    std::vector<ColumnIndexPair> unrelated_;
    std::vector<ColumnIndexPair> sfd_;

    void RegisterOptions();
    void LoadDataInternal() override;
    void MakeExecuteOptsAvailable() override;
    void ResetState() override;

    unsigned long long ExecuteInternal() override;

    [[nodiscard]] size_t Category(model::ColumnIndex col_ind, std::string const &val, size_t domain,
                                  bool skew, FrequencyHandler const &handler) const;

    void FillContingencyTable(Sample &smp, std::vector<model::TypedColumnData> const &data,
                              model::ColumnIndex col_i, model::ColumnIndex col_k,
                              std::vector<std::vector<long double>> &n_i_j,
                              std::vector<long double> &n_i, std::vector<long double> &n_j,
                              FrequencyHandler const &handler) const;

    bool TooMuchStructuralZeros(std::vector<std::vector<long double>> const &n_i_j,
                                model::ColumnIndex col_i, model::ColumnIndex col_k) const;

    long double CalculateChiSquared(std::vector<std::vector<long double>> const &n_i_j,
                                    std::vector<long double> const &n_i,
                                    std::vector<long double> const &n_j, model::ColumnIndex col_i,
                                    model::ColumnIndex col_k, long double sample_size) const;

    void Init(model::ColumnIndex columns);

    std::pair<model::ColumnIndex, model::ColumnIndex> SwitchIndicesIfNecessary(
            model::ColumnIndex ind1, model::ColumnIndex ind2, size_t card1, size_t card2);

    bool DetectAndEmplaceSfd(Sample const &smp);

    void SkewHandling(model::ColumnIndex col_i, model::ColumnIndex col_k,
                      std::vector<model::TypedColumnData> const &data, Sample &smp,
                      FrequencyHandler const &handler);

    bool IsTrivialCase(model::ColumnIndex col_ind, std::vector<bool> &is_soft_or_trivial,
                       FrequencyHandler const &handler, size_t row_count);

    void CheckForCorrelation(model::ColumnIndex col_i, model::ColumnIndex col_k,
                             std::vector<model::TypedColumnData> const &data, Sample &smp,
                             FrequencyHandler const &handler);

    bool ChiSquaredTest(model::ColumnIndex col_i, model::ColumnIndex col_k, Sample &smp,
                        std::vector<std::vector<long double>> const &n_i_j,
                        std::vector<long double> const &n_i,
                        std::vector<long double> const &n_j) const;

public:
    std::vector<model::ColumnIndex> const &GetSoftKeys() const;
    std::vector<model::ColumnIndex> const &GetTrivialColumns() const;
    std::vector<ColumnIndexPair> const &GetCorrelations() const;
    std::vector<ColumnIndexPair> const &GetUnrelated() const;
    std::vector<ColumnIndexPair> const &GetSFD() const;
    SFDAlgorithm();
};
}  // namespace algos
