#pragma once

#include <list>
#include <set>
#include <memory>
#include <utility>

#include <util/VerticalMap.h>
#include "ProfilingContext.h"
#include "DependencyStrategy.h"
#include "VerticalInfo.h"
#include "DependencyCandidate.h"
#include "Vertical.h"
#include "RelationalSchema.h"
class SearchSpace : public std::enable_shared_from_this<SearchSpace> {
private:
    ProfilingContext* context_;
    std::unique_ptr<DependencyStrategy> strategy_;
    std::unique_ptr<VerticalMap<VerticalInfo>> globalVisitees_;
    std::set<DependencyCandidate, std::function<bool (DependencyCandidate const&, DependencyCandidate const&)>> launchPads_;
    std::unique_ptr<VerticalMap<DependencyCandidate>> launchPadIndex_;
    std::list<DependencyCandidate> deferredLaunchPads_;
    std::unique_ptr<VerticalMap<Vertical>> scope_;
    double sampleBoost_;
    int recursionDepth_;
    bool isAscendRandomly_ = false;

    int numNested = 0;


    void discover(std::unique_ptr<VerticalMap<VerticalInfo>> localVisitees);
    std::optional<DependencyCandidate> pollLaunchPad(VerticalMap<VerticalInfo>* localVisitees);
    void escapeLaunchPad(Vertical const& launchPad, std::vector<Vertical> pruningSupersets, VerticalMap<VerticalInfo>* localVisitees);
    void returnLaunchPad(DependencyCandidate const& launchPad, bool isDefer);

    bool ascend(DependencyCandidate const& launchPad, std::shared_ptr<VerticalMap<std::shared_ptr<VerticalInfo>>> localVisitees);
    void checkEstimate(std::shared_ptr<DependencyStrategy> strategy, DependencyCandidate const& traversalCandidate);
    void trickleDown(std::shared_ptr<Vertical> mainPeak, double mainPeakError, std::shared_ptr<VerticalMap<std::shared_ptr<VerticalInfo>>> localVisitees);
    std::shared_ptr<Vertical> trickleDownFrom(DependencyCandidate & minDepCandidate,
            std::shared_ptr<DependencyStrategy> strategy, std::shared_ptr<VerticalMap<std::shared_ptr<VerticalInfo>>> allegedMinDeps,
            std::unordered_set<Vertical> & allegedNonDeps, std::shared_ptr<VerticalMap<std::shared_ptr<VerticalInfo>>> localVisitees,
            std::shared_ptr<VerticalMap<std::shared_ptr<VerticalInfo>>> globalVisitees, double boostFactor);
    static void requireMinimalDependency(std::shared_ptr<DependencyStrategy> strategy, std::shared_ptr<Vertical> minDependency);
    static std::list<std::shared_ptr<Vertical>> getSubsetDeps(Vertical const& vertical, VerticalMap<VerticalInfo>* verticalInfos); // no idea of return type
    static bool isImpliedByMinDep(Vertical const& vertical, VerticalMap<VerticalInfo>* verticalInfos);
    static bool isKnownNonDependency(Vertical const& vertical, VerticalMap<VerticalInfo>* verticalInfos);
    static std::string formatArityHistogram() = delete;
    static std::string formatArityHistogram(VerticalMap<int*>) = delete;

public:
    unsigned long long nanosSmartConstructing = 0;
    unsigned long long pollingLaunchPads = 0;
    unsigned long long ascending = 0;
    unsigned long long tricklingDown = 0;
    unsigned long long tricklingDownPart = 0;
    unsigned long long tricklingDownFrom = 0;
    unsigned long long returningLaunchPad = 0;

    bool isInitialized_ = false;
    int id_;

    SearchSpace(int id, std::unique_ptr<DependencyStrategy> strategy, std::unique_ptr<VerticalMap<std::shared_ptr<Vertical>>> scope,
            std::unique_ptr<VerticalMap<std::shared_ptr<VerticalInfo>>> globalVisitees, RelationalSchema const* schema,
            std::function<bool (DependencyCandidate const&, DependencyCandidate const&)> const& dependencyCandidateComparator,
            int recursionDepth, double sampleBoost) :
                strategy_(std::move(strategy)), globalVisitees_(std::move(globalVisitees)),
                launchPads_(dependencyCandidateComparator),
                launchPadIndex_(std::make_unique<VerticalMap<std::shared_ptr<DependencyCandidate>>>(schema)),
                scope_(std::move(scope)), sampleBoost_(sampleBoost), recursionDepth_(recursionDepth), id_(id) {}

    SearchSpace(int id, std::unique_ptr<DependencyStrategy> strategy, RelationalSchema const* schema,
            std::function<bool (DependencyCandidate const&, DependencyCandidate const&)> const& dependencyCandidateComparator):
            SearchSpace(id, std::move(strategy), nullptr,
                        std::make_unique<VerticalMap<std::shared_ptr<VerticalInfo>>>(schema),
                        schema, dependencyCandidateComparator, 0, 1) {}

    void ensureInitialized();
    void discover() { discover(nullptr); }
    void addLaunchPad(DependencyCandidate const& launchPad);
    void setContext(ProfilingContext* context)  {
        context_ = context;
        strategy_->context_ = context;
    }
    ProfilingContext* getContext() { return context_; }
    unsigned int getErrorCalcCount() { return strategy_->calcCount_; }
    void printStats() const;
};