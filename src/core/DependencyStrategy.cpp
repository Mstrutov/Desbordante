#include "DependencyStrategy.h"
#include "PLICache.h"

bool DependencyStrategy::shouldResample(Vertical const& vertical, double boostFactor) {

    if (context_->getConfiiguration().sampleSize <= 0 || vertical.getArity() < 1) return false;

    // Do we have an exact sample already?
    std::shared_ptr<AgreeSetSample> currentSample = context_->getAgreeSetSample(vertical);
    if (currentSample->isExact()) return false;

    // Get an estimate of the number of equality pairs in the vertical
    std::shared_ptr<PositionListIndex> pli = context_->getPLICache()->get(vertical);
    double nep = pli != nullptr
            ? pli->getNepAsLong()
            : currentSample->estimateAgreements(vertical) * context_->getColumnLayoutRelationData()->getNumTuplePairs();

    // Should the new sample be exact?
    if (nep <= context_->getConfiiguration().sampleSize * boostFactor) return true;

    // Will we achieve an improved sampling ratio?
    double newSamplingRatio = context_->getConfiiguration().sampleSize * boostFactor / nep;
    return newSamplingRatio >= 2 * currentSample->getSamplingRatio();
}