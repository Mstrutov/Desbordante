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
    using DependencyCandidateComp =
        std::function<bool(DependencyCandidate const&, DependencyCandidate const&)>;
    ProfilingContext* context_;
    std::unique_ptr<DependencyStrategy> strategy_;
    std::unique_ptr<util::VerticalMap<VerticalInfo>> local_visitees_ = nullptr;
    std::unique_ptr<util::VerticalMap<VerticalInfo>> global_visitees_;
    std::set<DependencyCandidate, DependencyCandidateComp> launch_pads_;
    std::unique_ptr<util::VerticalMap<DependencyCandidate>> launch_pad_index_;
    std::list<DependencyCandidate> deferred_launch_pads_;
    std::unique_ptr<util::VerticalMap<Vertical>> scope_;
    double sample_boost_;
    int recursion_depth_;
    bool is_ascend_randomly_ = false;

    int num_nested_ = 0;

    // void Discover(std::unique_ptr<VerticalMap<VerticalInfo>> localVisitees);
    std::optional<DependencyCandidate> PollLaunchPad();
    void EscapeLaunchPad(Vertical const& hitting_set_candidate,
                         std::vector<Vertical> pruning_supersets);
    void ReturnLaunchPad(DependencyCandidate const& launch_pad, bool is_defer);

    bool Ascend(DependencyCandidate const& launch_pad);
    void CheckEstimate(DependencyStrategy* strategy,
                       DependencyCandidate const& traversal_candidate);
    void TrickleDown(Vertical const& main_peak, double main_peak_error);
    std::optional<Vertical> TrickleDownFrom(DependencyCandidate min_dep_candidate,
                                            DependencyStrategy* strategy,
                                            util::VerticalMap<VerticalInfo>* alleged_min_deps,
                                            std::unordered_set<Vertical>& alleged_non_deps,
                                            util::VerticalMap<VerticalInfo>* global_visitees,
                                            double boost_factor);

    // CAREFUL: resets globalVisitees_, therefore SearchSpace could become invalidated
    std::unique_ptr<util::VerticalMap<VerticalInfo>> MoveOutGlobalVisitees() {
        return std::move(global_visitees_);
    }
    // CAREFUL: resets localVisitees_, therefore SearchSpace could become invalidated
    std::unique_ptr<util::VerticalMap<VerticalInfo>> MoveOutLocalVisitees() {
        return std::move(local_visitees_);
    }
    void MoveInLocalVisitees(std::unique_ptr<util::VerticalMap<VerticalInfo>> local_visitees) {
        local_visitees_ = std::move(local_visitees);
    }

    static void RequireMinimalDependency(DependencyStrategy* strategy,
                                         Vertical const& min_dependency);
    static std::vector<Vertical> GetSubsetDeps(Vertical const& vertical,
                                               util::VerticalMap<VerticalInfo>* vertical_infos);
    static bool IsImpliedByMinDep(Vertical const& vertical,
                                  util::VerticalMap<VerticalInfo>* vertical_infos);
    static bool IsKnownNonDependency(Vertical const& vertical,
                                     util::VerticalMap<VerticalInfo>* vertical_infos);
    static std::string FormatArityHistogram() = delete;
    static std::string FormatArityHistogram(util::VerticalMap<int*>) = delete;

public:
    unsigned long long nanos_smart_constructing_ = 0;
    unsigned long long polling_launch_pads_ = 0;
    unsigned long long ascending_ = 0;
    unsigned long long trickling_down_ = 0;
    unsigned long long trickling_down_part_ = 0;
    unsigned long long trickling_down_from_ = 0;
    unsigned long long returning_launch_pad_ = 0;

    bool is_initialized_ = false;
    int id_;

    SearchSpace(int id, std::unique_ptr<DependencyStrategy> strategy,
                std::unique_ptr<util::VerticalMap<Vertical>> scope,
                std::unique_ptr<util::VerticalMap<VerticalInfo>> global_visitees,
                RelationalSchema const* schema,
                DependencyCandidateComp const& dependency_candidate_comparator,
                int recursion_depth, double sample_boost)
        : strategy_(std::move(strategy)), global_visitees_(std::move(global_visitees)),
          launch_pads_(dependency_candidate_comparator),
          launch_pad_index_(std::make_unique<util::VerticalMap<DependencyCandidate>>(schema)),
          scope_(std::move(scope)), sample_boost_(sample_boost), recursion_depth_(recursion_depth),
          id_(id) {}

    SearchSpace(int id, std::unique_ptr<DependencyStrategy> strategy,
                RelationalSchema const* schema,
                DependencyCandidateComp const& dependency_candidate_comparator)
        : SearchSpace(id, std::move(strategy), nullptr,
                      std::make_unique<util::VerticalMap<VerticalInfo>>(schema),
                      schema, dependency_candidate_comparator, 0, 1) {}

    void EnsureInitialized();
    void Discover();
    void AddLaunchPad(DependencyCandidate const& launch_pad);
    void SetContext(ProfilingContext* context) {
        context_ = context;
        strategy_->context_ = context;
    }
    ProfilingContext* GetContext() { return context_; }
    unsigned int GetErrorCalcCount() { return strategy_->calc_count_; }
    void PrintStats() const;
};
