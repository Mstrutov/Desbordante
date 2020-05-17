#pragma once

#include <list>
#include <set>

#include <util/VerticalMap.h>
#include "core/ProfilingContext.h"
#include "core/DependencyStrategy.h"
#include "core/VerticalInfo.h"
#include "core/DependencyCandidate.h"
#include "model/Vertical.h"
#include "model/RelationalSchema.h"
class SearchSpace {
private:
    std::shared_ptr<ProfilingContext> context_;
    std::shared_ptr<DependencyStrategy> strategy_;
    VerticalMap<std::shared_ptr<VerticalInfo>> globalVisitees_;              //should be stored as unique_ptr to avoid huge memory chunk allocation problems?
    std::set<DependencyCandidate, std::function<bool (DependencyCandidate const&, DependencyCandidate const&)>> launchPads_;
    VerticalMap<std::shared_ptr<DependencyCandidate>> launchPadIndex_;
    std::list<DependencyCandidate> deferredLaunchPads_;
    std::unique_ptr<VerticalMap<std::shared_ptr<Vertical>>> scope_;
    double sampleBoost_;
    int recursionDepth_;
    bool isAscendRandomly_ = false;

public:
    bool isInitialized_ = false;
    int id_;

    // check if globalVisitees should be stored by shared_ptr, unique_ptr or value
    SearchSpace(int id, std::shared_ptr<DependencyStrategy> strategy, std::unique_ptr<VerticalMap<std::shared_ptr<Vertical>>> scope,
            VerticalMap<std::shared_ptr<VerticalInfo>> const & globalVisitees, std::shared_ptr<RelationalSchema> schema,
            std::function<bool (DependencyCandidate const&, DependencyCandidate const&)> const& dependencyCandidateComparator,
            int recursionDepth, double sampleBoost) :
                id_(id), strategy_(strategy), scope_(std::move(scope)), globalVisitees_(globalVisitees), recursionDepth_(recursionDepth),
                sampleBoost_(sampleBoost), launchPadIndex_(schema), launchPads_(dependencyCandidateComparator) {}

    SearchSpace(int id, std::shared_ptr<DependencyStrategy> strategy, std::shared_ptr<RelationalSchema> schema,
            std::function<bool (DependencyCandidate const&, DependencyCandidate const&)> const& dependencyCandidateComparator):
            SearchSpace(id, strategy, nullptr, static_cast<VerticalMap<std::shared_ptr<VerticalInfo>>>(schema), schema, dependencyCandidateComparator, 0, 1) {}

            //void ensureInitialized() {strategy_->ensureI}
    void addLaunchPad(DependencyCandidate const& launchPad);

};