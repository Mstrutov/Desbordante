#include <unordered_map>
#include "KeyG1Strategy.h"
#include "SearchSpace.h"
#include "PLICache.h"

double KeyG1Strategy::calculateKeyError(std::shared_ptr<PositionListIndex> pli) {
    return calculateKeyError(pli->getNepAsLong());
}

double KeyG1Strategy::calculateKeyError(double numViolatingTuplePairs) {
    unsigned long long numTuplePairs = context_->relationData_->getNumTuplePairs();
    if (numTuplePairs == 0) return 0;
    double g1 = numViolatingTuplePairs / numTuplePairs;
    return round(g1);
}

void KeyG1Strategy::ensureInitialized(std::shared_ptr<SearchSpace> searchSpace) {
    if (searchSpace->isInitialized_) return;

    for (auto& column : context_->getSchema()->getColumns()) {
        if (DependencyStrategy::isIrrelevantColumn(column)) continue;

        searchSpace->addLaunchPad(createDependencyCandidate(std::make_shared<Vertical>(static_cast<Vertical>(*column))));
    }

    searchSpace->isInitialized_ = true;
}

double KeyG1Strategy::calculateError(std::shared_ptr<Vertical> keyCandidate) {
    auto pli = context_->pliCache_->getOrCreateFor(*keyCandidate, *context_);
    double error = calculateKeyError(pli);
    calcCount_++;
    return error;
}

ConfidenceInterval KeyG1Strategy::calculateKeyError(ConfidenceInterval const &estimatedQualityPairs) {
    return ConfidenceInterval(calculateKeyError(estimatedQualityPairs.getMin()),
            calculateKeyError(estimatedQualityPairs.getMean()),
            calculateKeyError(estimatedQualityPairs.getMax()));
}

DependencyCandidate KeyG1Strategy::createDependencyCandidate(std::shared_ptr<Vertical> vertical) {
    if (vertical->getArity() == 1) {
        auto pli = context_->pliCache_->getOrCreateFor(*vertical, *context_);
        double keyError = calculateKeyError(pli->getNepAsLong());
        return DependencyCandidate(vertical, ConfidenceInterval(keyError), true);
    }

    if (context_->agreeSetSamples_ == nullptr) {
        return DependencyCandidate(vertical, ConfidenceInterval(0, .5, 1), false);
    }

    std::shared_ptr<AgreeSetSample> agreeSetSample = context_->getAgreeSetSample(vertical);
    ConfidenceInterval estimatedEqualityPairs = agreeSetSample
            ->estimateAgreements(vertical, context_->configuration_.estimateConfidence)
            .multiply(context_->relationData_->getNumTuplePairs());
    ConfidenceInterval keyError = calculateKeyError(estimatedEqualityPairs);
    return DependencyCandidate(vertical, keyError, false);
}

void KeyG1Strategy::registerDependency(std::shared_ptr<Vertical> vertical, double error,
                                      DependencyConsumer const& discoveryUnit) {
    discoveryUnit.registerUcc(vertical, error, 0); // TODO: calculate score?
}
