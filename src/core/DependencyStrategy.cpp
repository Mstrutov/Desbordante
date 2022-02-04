#include "DependencyStrategy.h"
#include "PLICache.h"

bool DependencyStrategy::ShouldResample(Vertical const& vertical, double boost_factor) const {
    if (context_->GetConfiguration().sample_size <= 0 || vertical.GetArity() < 1) return false;

    // Do we have an exact sample already?
    auto current_sample = context_->GetAgreeSetSample(vertical);
    if (current_sample->IsExact()) return false;

    // Get an estimate of the number of equality pairs in the vertical
    util::PositionListIndex* pli = context_->GetPliCache()->Get(vertical);
    double nep = pli != nullptr
            ? pli->GetNepAsLong()
            : current_sample->EstimateAgreements(vertical) *
            context_->GetColumnLayoutRelationData()->GetNumTuplePairs();

    // Should the new sample be exact?
    if (nep <= context_->GetConfiguration().sample_size * boost_factor) return true;

    // Will we achieve an improved sampling ratio?
    double new_sampling_ratio = context_->GetConfiguration().sample_size * boost_factor / nep;
    return new_sampling_ratio >= 2 * current_sample->GetSamplingRatio();
}
