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
    using DependencyCandidateComp = std::function<bool (DependencyCandidate const&,
                                                        DependencyCandidate const&)>;
    ProfilingContext* context_;
    std::unique_ptr<DependencyStrategy> strategy_;
    std::unique_ptr<util::VerticalMap<VerticalInfo>> localVisitees_ = nullptr;
    std::unique_ptr<util::VerticalMap<VerticalInfo>> globalVisitees_;
    std::set<DependencyCandidate, DependencyCandidateComp> launchPads_;
    std::unique_ptr<util::VerticalMap<DependencyCandidate>> launchPadIndex_;
    std::list<DependencyCandidate> deferredLaunchPads_;
    std::unique_ptr<util::VerticalMap<Vertical>> scope_;
    double sampleBoost_;
    int recursionDepth_;
    bool isAscendRandomly_ = false;

    int numNested = 0;


    // void discover(std::unique_ptr<VerticalMap<VerticalInfo>> localVisitees);
    std::optional<DependencyCandidate> pollLaunchPad();
    void escapeLaunchPad(Vertical const& launchPad, std::vector<Vertical> pruningSupersets);
    void returnLaunchPad(DependencyCandidate const& launchPad, bool isDefer);

    bool ascend(DependencyCandidate const& launchPad);
    void checkEstimate(DependencyStrategy* strategy, DependencyCandidate const& traversalCandidate);
    void trickleDown(Vertical const& mainPeak, double mainPeakError);
    std::optional<Vertical> trickleDownFrom(
            DependencyCandidate minDepCandidate, DependencyStrategy* strategy,
            util::VerticalMap<VerticalInfo>* allegedMinDeps,
            std::unordered_set<Vertical> & allegedNonDeps,
            util::VerticalMap<VerticalInfo>* globalVisitees, double boostFactor);

    // CAREFUL: resets globalVisitees_, therefore SearchSpace could become invalidated
    std::unique_ptr<util::VerticalMap<VerticalInfo>> moveOutGlobalVisitees() {
        return std::move(globalVisitees_);
    }
    // CAREFUL: resets localVisitees_, therefore SearchSpace could become invalidated
    std::unique_ptr<util::VerticalMap<VerticalInfo>> moveOutLocalVisitees() {
        return std::move(localVisitees_);
    }
    void moveInLocalVisitees(std::unique_ptr<util::VerticalMap<VerticalInfo>> localVisitees) {
        localVisitees_ = std::move(localVisitees);
    }


    static void requireMinimalDependency(DependencyStrategy* strategy,
                                         Vertical const& minDependency);
    static std::vector<Vertical> getSubsetDeps(Vertical const& vertical,
                                               util::VerticalMap<VerticalInfo>* verticalInfos);
    static bool isImpliedByMinDep(Vertical const& vertical,
                                  util::VerticalMap<VerticalInfo>* verticalInfos);
    static bool isKnownNonDependency(Vertical const& vertical,
                                     util::VerticalMap<VerticalInfo>* verticalInfos);
    static std::string formatArityHistogram() = delete;
    static std::string formatArityHistogram(util::VerticalMap<int*>) = delete;

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

    SearchSpace(int id, std::unique_ptr<DependencyStrategy> strategy,
                std::unique_ptr<util::VerticalMap<Vertical>> scope,
                std::unique_ptr<util::VerticalMap<VerticalInfo>> globalVisitees,
                RelationalSchema const* schema,
                DependencyCandidateComp const& dependencyCandidateComparator,
                int recursionDepth, double sampleBoost)
        : strategy_(std::move(strategy)), globalVisitees_(std::move(globalVisitees)),
          launchPads_(dependencyCandidateComparator),
          launchPadIndex_(std::make_unique<util::VerticalMap<DependencyCandidate>>(schema)),
          scope_(std::move(scope)), sampleBoost_(sampleBoost), recursionDepth_(recursionDepth),
          id_(id) {}

    SearchSpace(int id, std::unique_ptr<DependencyStrategy> strategy,
                RelationalSchema const* schema,
                DependencyCandidateComp const& dependencyCandidateComparator)
        : SearchSpace(id, std::move(strategy), nullptr,
                      std::make_unique<util::VerticalMap<VerticalInfo>>(schema),
                      schema, dependencyCandidateComparator, 0, 1) {}

    void ensureInitialized();
    void discover();
    void addLaunchPad(DependencyCandidate const& launchPad);
    void setContext(ProfilingContext* context)  {
        context_ = context;
        strategy_->context_ = context;
    }
    ProfilingContext* getContext() { return context_; }
    unsigned int getErrorCalcCount() { return strategy_->calcCount_; }
    void printStats() const;
};
