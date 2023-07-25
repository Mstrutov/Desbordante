#include "key_g1_strategy.h"

#include <unordered_map>

#include "pli_cache.h"
#include "search_space.h"

double KeyG1Strategy::CalculateKeyError(util::PositionListIndex* pli) const {
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

    for (auto& column : context_->GetSchema()->GetColumns()) {
        if (IsIrrelevantColumn(column->GetIndex())) continue;

        search_space->AddLaunchPad(CreateDependencyCandidate(static_cast<Vertical>(*column)));
    }

    search_space->is_initialized_ = true;
}

double KeyG1Strategy::CalculateError(Vertical const& key_candidate) const {
    auto pli = context_->GetPliCache()->GetOrCreateFor(key_candidate, context_);
    auto pli_pointer = std::holds_alternative<util::PositionListIndex*>(pli)
                           ? std::get<util::PositionListIndex*>(pli)
                           : std::get<std::unique_ptr<util::PositionListIndex>>(pli).get();
    double error = CalculateKeyError(pli_pointer);
    calc_count_++;
    return error;
}

util::ConfidenceInterval KeyG1Strategy::CalculateKeyError(util::ConfidenceInterval const& num_violations) const {
    return util::ConfidenceInterval(CalculateKeyError(num_violations.GetMin()),
                                    CalculateKeyError(num_violations.GetMean()),
                                    CalculateKeyError(num_violations.GetMax()));
}

DependencyCandidate KeyG1Strategy::CreateDependencyCandidate(Vertical const& vertical) const {
    if (vertical.GetArity() == 1) {
        auto pli = context_->GetPliCache()->GetOrCreateFor(vertical, context_);
        auto pli_pointer = std::holds_alternative<util::PositionListIndex*>(pli)
                               ? std::get<util::PositionListIndex*>(pli)
                               : std::get<std::unique_ptr<util::PositionListIndex>>(pli).get();
        double key_error = CalculateKeyError(pli_pointer->GetNepAsLong());
        return DependencyCandidate(vertical, util::ConfidenceInterval(key_error), true);
    }

    if (context_->IsAgreeSetSamplesEmpty()) {
        return DependencyCandidate(vertical, util::ConfidenceInterval(0, .5, 1), false);
    }

    auto agree_set_sample = context_->GetAgreeSetSample(vertical);
    util::ConfidenceInterval estimated_equality_pairs =
            agree_set_sample
                    ->EstimateAgreements(vertical, context_->GetParameters().estimate_confidence)
                    .Multiply(context_->GetColumnLayoutRelationData()->GetNumTuplePairs());
    util::ConfidenceInterval key_error = CalculateKeyError(estimated_equality_pairs);
    return DependencyCandidate(vertical, key_error, false);
}

void KeyG1Strategy::RegisterDependency(Vertical const& vertical, double error,
                                       DependencyConsumer const& discovery_unit) const {
    discovery_unit.RegisterUcc(vertical, error, 0); // TODO: calculate score?
}

std::unique_ptr<DependencyStrategy> KeyG1Strategy::CreateClone() {
    return std::make_unique<KeyG1Strategy>((max_dependency_error_ + min_non_dependency_error_) / 2,
                                           (max_dependency_error_ - min_non_dependency_error_) / 2);
}
