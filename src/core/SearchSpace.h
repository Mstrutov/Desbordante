#pragma once

#include <list>
#include <set>
#include <memory>

#include <util/VerticalMap.h>
#include "core/ProfilingContext.h"
#include "core/DependencyStrategy.h"
#include "core/VerticalInfo.h"
#include "core/DependencyCandidate.h"
#include "model/Vertical.h"
#include "model/RelationalSchema.h"
class SearchSpace : public std::enable_shared_from_this<SearchSpace> {
private:
    std::shared_ptr<ProfilingContext> context_;
    std::shared_ptr<DependencyStrategy> strategy_;
    std::shared_ptr<VerticalMap<std::shared_ptr<VerticalInfo>>> globalVisitees_;              //should be stored as unique_ptr to avoid huge memory chunk allocation problems?
    std::set<DependencyCandidate, std::function<bool (DependencyCandidate const&, DependencyCandidate const&)>> launchPads_;
    VerticalMap<std::shared_ptr<DependencyCandidate>> launchPadIndex_;
    std::list<DependencyCandidate> deferredLaunchPads_;
    std::unique_ptr<VerticalMap<std::shared_ptr<Vertical>>> scope_;
    double sampleBoost_;
    int recursionDepth_;
    bool isAscendRandomly_ = false;

    void discover(std::shared_ptr<VerticalMap<std::shared_ptr<VerticalInfo>>> localVisitees);
    std::shared_ptr<DependencyCandidate> pollLaunchPad(std::shared_ptr<VerticalMap<std::shared_ptr<VerticalInfo>>> localVisitees);
    void escapeLaunchPad(std::shared_ptr<Vertical> lanchPad, std::list<std::shared_ptr<Vertical>>&& pruningSupersets, std::shared_ptr<VerticalMap<std::shared_ptr<VerticalInfo>>> localVisitees);
    void returnLaunchPad(DependencyCandidate const& launchPad, bool isDefer);
    bool ascend(DependencyCandidate const& launchPad, std::shared_ptr<VerticalMap<std::shared_ptr<VerticalInfo>>> localVisitees);
    void checkEstimate(std::shared_ptr<DependencyStrategy> strategy, DependencyCandidate const& traversalCandidate);
    void trickleDown(std::shared_ptr<Vertical> mainPeak, double mainPeakError, std::shared_ptr<VerticalMap<std::shared_ptr<VerticalInfo>>> localVisitees);
    std::shared_ptr<Vertical> trickleDownFrom(DependencyCandidate const& minDepCandidate,
            std::shared_ptr<DependencyStrategy> strategy, std::shared_ptr<VerticalMap<std::shared_ptr<VerticalInfo>>> allegedMinDeps,
            std::set<std::shared_ptr<Vertical>> const& allegedNonDeps, std::shared_ptr<VerticalMap<std::shared_ptr<VerticalInfo>>> localVisitees,
            std::shared_ptr<VerticalMap<std::shared_ptr<VerticalInfo>>> globalVisitees, double boostFactor);
    static void requireMinimalDependency(std::shared_ptr<DependencyStrategy> strategy, std::shared_ptr<Vertical> minDependency);
    static void getSubsetDeps() = delete; // no idea of return type
    static bool isImpliedByMinDep(std::shared_ptr<Vertical> vertical, std::shared_ptr<VerticalMap<std::shared_ptr<VerticalInfo>>> verticalInfos);
    static bool isKnownNonDependency(std::shared_ptr<Vertical> vertical, std::shared_ptr<VerticalMap<std::shared_ptr<VerticalInfo>>> verticalInfos);
    static std::string formatArityHistogram() = delete;
    static std::string formatArityHistogram(VerticalMap<int*>) = delete;

public:
    bool isInitialized_ = false;
    int id_;

    // check if globalVisitees should be stored by shared_ptr, unique_ptr or value
    SearchSpace(int id, std::shared_ptr<DependencyStrategy> strategy, std::unique_ptr<VerticalMap<std::shared_ptr<Vertical>>> scope,
            std::shared_ptr<VerticalMap<std::shared_ptr<VerticalInfo>>> globalVisitees, std::shared_ptr<RelationalSchema> schema,
            std::function<bool (DependencyCandidate const&, DependencyCandidate const&)> const& dependencyCandidateComparator,
            int recursionDepth, double sampleBoost) :
                id_(id), strategy_(strategy), scope_(std::move(scope)), globalVisitees_(globalVisitees), recursionDepth_(recursionDepth),
                sampleBoost_(sampleBoost), launchPadIndex_(schema), launchPads_(dependencyCandidateComparator) {}

    SearchSpace(int id, std::shared_ptr<DependencyStrategy> strategy, std::shared_ptr<RelationalSchema> schema,
            std::function<bool (DependencyCandidate const&, DependencyCandidate const&)> const& dependencyCandidateComparator):
            SearchSpace(id, strategy, nullptr, std::static_pointer_cast<VerticalMap<std::shared_ptr<VerticalInfo>>>(schema), schema, dependencyCandidateComparator, 0, 1) {}

    void ensureInitialized() { strategy_->ensureInitialized(shared_from_this()); }
    void discover() { discover(nullptr); }
    void addLaunchPad(DependencyCandidate const& launchPad);
    void setContext(std::shared_ptr<ProfilingContext> context)  {
        context_ = context;
        strategy_->context_ = context;
    }
    std::shared_ptr<ProfilingContext> getContext() { return context_; }
};