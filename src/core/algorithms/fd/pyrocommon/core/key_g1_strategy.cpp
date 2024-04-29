#include "key_g1_strategy.h"

#include <unordered_map>

#include "../model/pli_cache.h"
#include "search_space.h"

double KeyG1Strategy::CalculateKeyError(model::PositionListIndex* pli) const {
    return CalculateKeyError(pli->GetNepAsLong());
}

double KeyG1Strategy::CalculateKeyError(double num_violating_tuple_pairs) const {
    unsigned long long num_tuple_pairs =
            context_->GetColumnLayoutRelationData()->GetNumTuplePairs();
    if (num_tuple_pairs == 0) return 0;
    double g1 = num_violating_tuple_pairs / num_tuple_pairs;
    return Round(g1);
}

void KeyG1Strategy::EnsureInitialized(SearchSpace* search_space) const {
    if (search_space->is_initialized_) return;

    for (auto const& column : context_->GetSchema()->GetColumns()) {
        if (IsIrrelevantColumn(column.GetIndex())) continue;

        search_space->AddLaunchPad(CreateDependencyCandidate(static_cast<Vertical>(column)));
    }

    search_space->is_initialized_ = true;
}

double KeyG1Strategy::CalculateError(Vertical const& key_candidate) const {
    auto pli = context_->GetPliCache()->GetOrCreateFor(key_candidate, context_);
    auto pli_pointer = std::holds_alternative<model::PositionListIndex*>(pli)
                               ? std::get<model::PositionListIndex*>(pli)
                               : std::get<std::unique_ptr<model::PositionListIndex>>(pli).get();
    double error = CalculateKeyError(pli_pointer);
    calc_count_++;
    return error;
}

model::ConfidenceInterval KeyG1Strategy::CalculateKeyError(
        model::ConfidenceInterval const& num_violations) const {
    return model::ConfidenceInterval(CalculateKeyError(num_violations.GetMin()),
                                     CalculateKeyError(num_violations.GetMean()),
                                     CalculateKeyError(num_violations.GetMax()));
}

DependencyCandidate KeyG1Strategy::CreateDependencyCandidate(Vertical const& vertical) const {
    if (vertical.GetArity() == 1) {
        auto pli = context_->GetPliCache()->GetOrCreateFor(vertical, context_);
        auto pli_pointer = std::holds_alternative<model::PositionListIndex*>(pli)
                                   ? std::get<model::PositionListIndex*>(pli)
                                   : std::get<std::unique_ptr<model::PositionListIndex>>(pli).get();
        double key_error = CalculateKeyError(pli_pointer->GetNepAsLong());
        return DependencyCandidate(vertical, model::ConfidenceInterval(key_error), true);
    }

    if (context_->IsAgreeSetSamplesEmpty()) {
        return DependencyCandidate(vertical, model::ConfidenceInterval(0, .5, 1), false);
    }

    auto agree_set_sample = context_->GetAgreeSetSample(vertical);
    model::ConfidenceInterval estimated_equality_pairs =
            agree_set_sample
                    ->EstimateAgreements(vertical, context_->GetParameters().estimate_confidence)
                    .Multiply(context_->GetColumnLayoutRelationData()->GetNumTuplePairs());
    model::ConfidenceInterval key_error = CalculateKeyError(estimated_equality_pairs);
    return DependencyCandidate(vertical, key_error, false);
}

void KeyG1Strategy::RegisterDependency(Vertical const& vertical, double error,
                                       DependencyConsumer const& discovery_unit) const {
    discovery_unit.RegisterUcc(vertical, error, 0);  // TODO: calculate score?
}

std::unique_ptr<DependencyStrategy> KeyG1Strategy::CreateClone() {
    return std::make_unique<KeyG1Strategy>((max_dependency_error_ + min_non_dependency_error_) / 2,
                                           (max_dependency_error_ - min_non_dependency_error_) / 2);
}
