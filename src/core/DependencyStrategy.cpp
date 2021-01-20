#include "DependencyStrategy.h"
#include "PLICache.h"

bool DependencyStrategy::shouldResample(std::shared_ptr<Vertical> vertical, double boostFactor) {

    if (context_->configuration_.sampleSize <= 0 || vertical->getArity() < 1) return false;

    // Do we have an exact sample already?
    std::shared_ptr<AgreeSetSample> currentSample = context_->getAgreeSetSample(vertical);
    if (currentSample->isExact()) return false;

    // Get an estimate of the number of equality pairs in the vertical
    std::shared_ptr<PositionListIndex> pli = context_->pliCache_->get(*vertical);
    double nep = pli != nullptr ?
                pli->getNepAsLong() :       // TODO: getNepAsDouble; long -> long long
                currentSample->estimateAgreements(vertical) * context_->relationData_->getNumTuplePairs();

    // Should the new sample be exact?
    if (nep <= context_->configuration_.sampleSize * boostFactor) return true;

    // Will we achieve an improved sampling ratio?
    double newSamplingRatio = context_->configuration_.sampleSize * boostFactor / nep;
    return newSamplingRatio >= 2 * currentSample->getSamplingRatio();
}