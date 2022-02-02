#include "logging/easylogging++.h"
#include "SearchSpace.h"
#include <queue>
// TODO: extra careful with const& -> shared_ptr conversions via make_shared-smart pointer may delete the object - pass empty deleter [](*) {}

void SearchSpace::discover() {
    LOG(TRACE) << "Discovering in: " << static_cast<std::string>(*strategy_);
    while (true) {  // на второй итерации дропается
        auto now = std::chrono::system_clock::now();
        std::optional<DependencyCandidate> launchPad = pollLaunchPad();
        if (!launchPad.has_value()) break;

        if (localVisitees_ == nullptr) {
            localVisitees_ = std::make_unique<util::VerticalMap<VerticalInfo>>(context_->getSchema());
        }

        bool isDependencyFound = ascend(*launchPad);
        pollingLaunchPads += std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::system_clock::now() - now).count();
        returnLaunchPad(*launchPad, !isDependencyFound);
    }
}

std::optional<DependencyCandidate> SearchSpace::pollLaunchPad() {
    while (true) {
        if (launchPads_.empty()) {
            if (deferredLaunchPads_.empty()) return std::optional<DependencyCandidate>();

            launchPads_.insert(deferredLaunchPads_.begin(), deferredLaunchPads_.end());
            deferredLaunchPads_.clear();
        }

        auto launchPad = launchPads_.extract(launchPads_.begin()).value();

        // launchPads_.erase(launchPads_.begin());
        launchPadIndex_->remove(launchPad.vertical_);

        if (isImpliedByMinDep(launchPad.vertical_, globalVisitees_.get())
            || (localVisitees_ != nullptr && isImpliedByMinDep(launchPad.vertical_, localVisitees_.get()))) {
            launchPadIndex_->remove(launchPad.vertical_);
            LOG(TRACE) << "* Removing subset-pruned launch pad {" << launchPad.vertical_.toString() << '}';
            continue;
        }

        auto supersetEntries = globalVisitees_->getSupersetEntries(launchPad.vertical_);
        if (localVisitees_ != nullptr) {
            auto localSupersetEntries = localVisitees_->getSupersetEntries(launchPad.vertical_);
            auto endIterator = std::remove_if(localSupersetEntries.begin(), localSupersetEntries.end(),
                                              [](auto& entry) { return !entry.second->isPruningSubsets(); });

            std::for_each(localSupersetEntries.begin(), endIterator,
                          [&supersetEntries](auto &entry) { supersetEntries.push_back(entry); });
        }
        if (supersetEntries.empty()) return launchPad;
        LOG(TRACE) << boost::format{"* Escaping launchPad %1% from: %2%"}
            % launchPad.vertical_.toString() % "[UNIMPLEMENTED]";
        std::vector<Vertical> supersetVerticals;

        for (auto& entry : supersetEntries) {
            supersetVerticals.push_back(entry.first);
        }

        escapeLaunchPad(launchPad.vertical_, std::move(supersetVerticals));
    }
}

// this move looks legit IMO
void SearchSpace::escapeLaunchPad(Vertical const& launchPad,
                                  std::vector<Vertical> pruningSupersets) {
    std::transform(pruningSupersets.begin(), pruningSupersets.end(), pruningSupersets.begin(),
            [this](auto& superset) { return superset.invert().without(strategy_->getIrrelevantColumns()); } );

    std::function<bool (Vertical const&)> pruningFunction =
            [this, &launchPad] (auto const& hittingSetCandidate) -> bool {
        if (scope_ != nullptr && scope_->getAnySupersetEntry(hittingSetCandidate).second == nullptr) {
            return true;
        }

        auto launchPadCandidate = launchPad.Union(hittingSetCandidate);

        if ((localVisitees_ == nullptr && isImpliedByMinDep(launchPadCandidate, localVisitees_.get()))
            || isImpliedByMinDep(launchPadCandidate, globalVisitees_.get())) {
            return true;
        }

        if (launchPadIndex_->getAnySubsetEntry(launchPadCandidate).second != nullptr) {
            return true;
        }

        return false;
    };
    {
        std::string pruningSupersetsStr = "[";
        for (auto& pruningSuperset : pruningSupersets) {
            pruningSupersetsStr += pruningSuperset.toString();
        }
        pruningSupersetsStr += "]";
        LOG(TRACE) << boost::format{"Escaping %1% pruned by %2%"} % launchPad.toString() % pruningSupersetsStr;
    }
    auto hittingSet = context_->getSchema()->calculateHittingSet(
            std::move(pruningSupersets),
            boost::make_optional(pruningFunction)
        );
    {
        std::string hittingSetStr = "[";
        for (auto& el : hittingSet) {
            hittingSetStr += el.toString();
        }
        hittingSetStr += "]";
        LOG(TRACE) << boost::format{"* Evaluated hitting set: %1%"} % hittingSetStr;
    }
    for (auto& escaping : hittingSet) {

        auto escapedLaunchPadVertical = launchPad.Union(escaping);

        // assert, который не имплементнуть из-за трансформа
        LOG(TRACE) << "createDependencyCandidate while escaping launch pad";
        DependencyCandidate escapedLaunchPad = strategy_->createDependencyCandidate(escapedLaunchPadVertical);
        LOG(TRACE) << boost::format{"  Escaped: %1%"} % escapedLaunchPadVertical.toString();
        LOG(DEBUG) << boost::format{"  Proposed launch pad arity: %1% should be <= maxLHS: %2%"}
            % escapedLaunchPad.vertical_.getArity() % context_->getConfiguration().maxLHS;
        if (escapedLaunchPad.vertical_.getArity() <= context_->getConfiguration().maxLHS) {
            launchPads_.insert(escapedLaunchPad);
            launchPadIndex_->put(escapedLaunchPad.vertical_, std::make_unique<DependencyCandidate>(escapedLaunchPad));
        }
    }
}

void SearchSpace::addLaunchPad(const DependencyCandidate &launchPad) {
    launchPads_.insert(launchPad);
    launchPadIndex_->put(launchPad.vertical_, std::make_unique<DependencyCandidate>(launchPad));
}

void SearchSpace::returnLaunchPad(DependencyCandidate const &launchPad, bool isDefer) {
    if (isDefer && context_->getConfiguration().isDeferFailedLaunchPads) {
        deferredLaunchPads_.push_back(launchPad);
        LOG(TRACE) << boost::format{"Deferred seed %1%"} % launchPad.vertical_.toString();
    }
    else {
        launchPads_.insert(launchPad);
    }
    launchPadIndex_->put(launchPad.vertical_, std::make_unique<DependencyCandidate>(launchPad));
}

bool SearchSpace::ascend(DependencyCandidate const &launchPad) {

    auto now = std::chrono::system_clock::now();

    LOG(DEBUG) << boost::format{"===== Ascending from %1% ======"} % strategy_->format(launchPad.vertical_);

    if (strategy_->shouldResample(launchPad.vertical_, sampleBoost_)) {
        LOG(TRACE) << "Resampling.";
        context_->createFocusedSample(launchPad.vertical_, sampleBoost_);
    }

    DependencyCandidate traversalCandidate = launchPad;
    boost::optional<double> error;

    while (true) {
        LOG(TRACE) << boost::format{"-> %1%"} % traversalCandidate.vertical_.toString();

        if (context_->getConfiguration().isCheckEstimates) {
            checkEstimate(strategy_.get(), traversalCandidate);
        }

        if (traversalCandidate.isExact()) {
            // TODO: getError()?
            error = traversalCandidate.error_.get();

            bool canBeDependency = *error <= strategy_->maxDependencyError_;
            localVisitees_->put(traversalCandidate.vertical_, std::make_unique<VerticalInfo>(
                    canBeDependency,
                    false,
                    *error
                    ));
            if (canBeDependency) break;
        }
        else {
            if (traversalCandidate.error_.getMin() > strategy_->maxDependencyError_) {
                LOG(TRACE) << boost::format {"  Skipping check form %1% (estimated error: %2%)."}
                    % traversalCandidate.vertical_.toString() % traversalCandidate.error_;
                error.reset();
            }
            else {
                error = context_->getConfiguration().isEstimateOnly
                        ? traversalCandidate.error_.getMean()
                        : strategy_->calculateError(traversalCandidate.vertical_);
                // double errorDiff = *error - traversalCandidate.error_.getMean();

                localVisitees_->put(traversalCandidate.vertical_, std::make_unique<VerticalInfo>(
                        error <= strategy_->maxDependencyError_,
                        false,
                        *error
                        ));

                if (*error <= strategy_->maxDependencyError_) break;

                if (strategy_->shouldResample(traversalCandidate.vertical_, sampleBoost_)) {
                    LOG(TRACE) << "Resampling.";
                    context_->createFocusedSample(traversalCandidate.vertical_, sampleBoost_);
                }
            }
        }

        if (traversalCandidate.vertical_.getArity() >=
            context_->getColumnLayoutRelationData()->getNumColumns() - strategy_->getNumIrrelevantColumns()
            || traversalCandidate.vertical_.getArity() >= context_->getConfiguration().maxLHS) {
            break;
        }

        boost::optional<DependencyCandidate> nextCandidate;
        int numSeenElements = isAscendRandomly_ ? 1 : -1;
        for (auto& extensionColumn : context_->getSchema()->getColumns()) {
            if (traversalCandidate.vertical_.getColumnIndices()[extensionColumn->getIndex()]
                    || strategy_->isIrrelevantColumn(*extensionColumn)) {
                continue;
            }
            auto extendedVertical = traversalCandidate.vertical_.Union(static_cast<Vertical>(*extensionColumn));

            if (scope_ != nullptr && scope_->getSupersetEntries(extendedVertical).empty()) {
                continue;
            }

            bool isSubsetPruned = isImpliedByMinDep(extendedVertical, globalVisitees_.get());
            if (isSubsetPruned) {
                continue;
            }
            LOG(TRACE) << "createDependencyCandidate while ascending";
            DependencyCandidate extendedCandidate = strategy_->createDependencyCandidate(extendedVertical);

            if (!nextCandidate
                    || (numSeenElements == -1 && extendedCandidate.error_.getMean() < nextCandidate->error_.getMean())
                    || (numSeenElements != -1 && context_->nextInt(++numSeenElements) == 0)) {
                nextCandidate = extendedCandidate;
            }
        }

        if (nextCandidate) {
            traversalCandidate = *nextCandidate;
        } else {
            break;
        }
    }
    
    //std::cout << static_cast<std::string>(traversalCandidate) << std::endl;

    if (!error) {
        LOG(TRACE) << boost::format{"  Hit ceiling at %1%."} % traversalCandidate.vertical_.toString();
        error = strategy_->calculateError(traversalCandidate.vertical_);
        [[maybe_unused]] double errorDiff = *error - traversalCandidate.error_.getMean();
        LOG(TRACE) << boost::format{"  Checking candidate... actual error: %1%"} % *error;
    }
    ascending += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - now).count();

    if (*error <= strategy_->maxDependencyError_) {
        LOG(TRACE) << boost::format{"  Key peak in climbing phase e(%1%)=%2% -> Need to minimize."}
            % traversalCandidate.vertical_.toString() % *error;
        trickleDown(traversalCandidate.vertical_, *error);

        if (recursionDepth_ == 0) {
            assert(scope_ == nullptr);
            globalVisitees_->put(traversalCandidate.vertical_,
                                 std::make_unique<VerticalInfo>(VerticalInfo::forMinimalDependency()));
        }

        return true;
    } else {
        if (recursionDepth_ == 0) {
            assert(scope_ == nullptr);
            globalVisitees_->put(traversalCandidate.vertical_,
                                 std::make_unique<VerticalInfo>(VerticalInfo::forMaximalNonDependency()));
            LOG(DEBUG) << boost::format{"[---] %1% is maximum non-dependency (err=%2%)."}
                % traversalCandidate.vertical_.toString() % *error;
        } else {
            localVisitees_->put(traversalCandidate.vertical_,
                               std::make_unique<VerticalInfo>(VerticalInfo::forNonDependency()));
            LOG(DEBUG) << boost::format{"      %1% is local-maximum non-dependency (err=%2%)."}
                % traversalCandidate.vertical_.toString() % *error;
        }
    }
    return false;
}

void SearchSpace::checkEstimate([[maybe_unused]] DependencyStrategy* strategy,
                                [[maybe_unused]] DependencyCandidate const& traversalCandidate) {
    std::cout << "Stepped into method 'checkEstimate' - not implemented yet being a debug method\n";
}

void SearchSpace::trickleDown(Vertical const& mainPeak, double mainPeakError) {
    LOG(DEBUG) << boost::format{"====== Trickling down from %1% ======"} % mainPeak.toString();

    std::unordered_set<Vertical> maximalNonDeps;
    auto allegedMinDeps = std::make_unique<util::VerticalMap<VerticalInfo>>(context_->getSchema());
    auto peaksComparator = [](auto& candidate1, auto& candidate2) -> bool {
        return DependencyCandidate::arityComparator(candidate1, candidate2); };
    std::vector<DependencyCandidate> peaks;
    std::make_heap(peaks.begin(), peaks.end(), peaksComparator);
    peaks.emplace_back(mainPeak, util::ConfidenceInterval(mainPeakError), true);
    std::push_heap(peaks.begin(), peaks.end(), peaksComparator);
    std::unordered_set<Vertical> allegedNonDeps;

    auto now = std::chrono::system_clock::now();

    while (!peaks.empty()) {
        auto peak = peaks.front();

        auto subsetDeps = getSubsetDeps(peak.vertical_, allegedMinDeps.get());
        if (!subsetDeps.empty()) {
            std::pop_heap(peaks.begin(), peaks.end(), peaksComparator);
            peaks.pop_back();

            auto peakHittingSet = context_->getSchema()->calculateHittingSet(
                    std::move(subsetDeps), boost::optional<std::function<bool (Vertical const&)>>());
            std::unordered_set<Vertical> escapedPeakVerticals;

            for (auto& vertical : peakHittingSet) {
                escapedPeakVerticals.insert(peak.vertical_.without(vertical));
            }

            for (auto& escapedPeakVertical : escapedPeakVerticals) {
                if (escapedPeakVertical.getArity() > 0
                    && allegedNonDeps.find(escapedPeakVertical) == allegedNonDeps.end()) {

                    LOG(TRACE) << "createDependencyCandidate as an escaped peak while trickling down";
                    auto escapedPeak = strategy_->createDependencyCandidate(escapedPeakVertical);

                    if (escapedPeak.error_.getMean() > strategy_->minNonDependencyError_) {
                        allegedNonDeps.insert(escapedPeakVertical);
                        continue;
                    }
                    if (isKnownNonDependency(escapedPeakVertical, localVisitees_.get())
                        || isKnownNonDependency(escapedPeakVertical, globalVisitees_.get())) {
                        continue;
                    }
                    peaks.push_back(escapedPeak);
                    std::push_heap(peaks.begin(), peaks.end(), peaksComparator);
                }
            }
            continue;
        }
        auto allegedMinDep = trickleDownFrom(
                std::move(peak), strategy_.get(), allegedMinDeps.get(), allegedNonDeps,
                globalVisitees_.get(), sampleBoost_);
        if (!allegedMinDep.has_value()) {
            std::pop_heap(peaks.begin(), peaks.end(), peaksComparator);
            peaks.pop_back();
        }
    }

    LOG(DEBUG) << boost::format{"* %1% alleged minimum dependencies (%2%)"}
        % allegedMinDeps->getSize() % "UNIMPLEMENTED";

    int numUncertainMinDeps = 0;
    for (auto& [allegedMinDep, info] : allegedMinDeps->entrySet()) {
        if (info->isExtremal_ && !globalVisitees_->containsKey(allegedMinDep)) {
            LOG(DEBUG) << boost::format{"[%1%] Minimum dependency: %2% (error=%3%)"}
                % recursionDepth_ % allegedMinDep.toString() % info->error_;
            // TODO: Костыль -- info в нескольких местах должен храниться. ХЗ, кому он принадлежит, пока копирую
            globalVisitees_->put(allegedMinDep, std::make_unique<VerticalInfo>(*info));
            strategy_->registerDependency(allegedMinDep, info->error_, *context_);
        }
        if (!info->isExtremal_) {
            numUncertainMinDeps++;
        }
    }

    LOG(DEBUG) << boost::format{"* %1%/%2% alleged minimum dependencies might be non-minimal"}
        % numUncertainMinDeps % allegedMinDeps->getSize();

    auto allegedMinDepsSet = allegedMinDeps->keySet();
    // TODO: костыль: calculateHittingSet needs a list, but keySet returns an unordered_set
    // ещё и морока с transform и unordered_set - мб вообще в лист переделать.
    auto allegedMaxNonDepsHS = context_->getSchema()->calculateHittingSet(
            std::vector<Vertical>(allegedMinDepsSet.begin(), allegedMinDepsSet.end()),
            boost::optional<std::function<bool (Vertical const&)>>());
    std::unordered_set<Vertical> allegedMaxNonDeps;


    for (auto& minLeaveOutVertical : allegedMaxNonDepsHS) {
        allegedMaxNonDeps.insert(minLeaveOutVertical.invert(mainPeak));
    }

    LOG(DEBUG) << boost::format{"* %1% alleged maximum non-dependencies (%2%)"}
        % allegedMaxNonDeps.size() % "UNIMPLEMENTED";
    //std::transform(allegedMaxNonDepsHS.begin(), allegedMaxNonDepsHS.end(), allegedMaxNonDeps.begin(),
    //        [mainPeak](auto minLeaveOutVertical) { return minLeaveOutVertical->invert(*mainPeak); });

    // checking the consistency of all data structures
    if (auto allegedMinDepsKeySet = allegedMinDeps->keySet();
            !std::all_of(allegedMinDepsKeySet.begin(), allegedMinDepsKeySet.end(),
                [mainPeak](auto& vertical) -> bool { return mainPeak.contains(vertical); })) {
        throw std::runtime_error("Main peak should contain all alleged min dependencies");
    }
    if (!std::all_of(allegedMaxNonDeps.begin(), allegedMaxNonDeps.end(),
                     [mainPeak](auto& vertical) -> bool { return mainPeak.contains(vertical); })) {
        throw std::runtime_error("Main peak should contain all alleged max non-dependencies");
    }

    for (auto& allegedMaxNonDep : allegedMaxNonDeps) {
        if (allegedMaxNonDep.getArity() == 0) continue;

        if (maximalNonDeps.find(allegedMaxNonDep) != maximalNonDeps.end()
                || isKnownNonDependency(allegedMaxNonDep, localVisitees_.get())
                || isKnownNonDependency(allegedMaxNonDep, globalVisitees_.get())) {
            continue;
        }

        double error = context_->getConfiguration().isEstimateOnly
                ? strategy_->createDependencyCandidate(allegedMaxNonDep).error_.getMean()
                : strategy_->calculateError(allegedMaxNonDep);
        bool isNonDep = error > strategy_->minNonDependencyError_;
        LOG(TRACE) << boost::format {"* Alleged maximal non-dependency %1%: non-dep?: %2%, error: %3%"}
            % allegedMaxNonDep.toString() % isNonDep % error;
        if (isNonDep) {
            maximalNonDeps.insert(allegedMaxNonDep);
            localVisitees_->put(allegedMaxNonDep, std::make_unique<VerticalInfo>(VerticalInfo::forNonDependency()));
        } else {
            peaks.emplace_back(allegedMaxNonDep, util::ConfidenceInterval(error), true);
            std::push_heap(peaks.begin(), peaks.end(), peaksComparator);
        }
    }
    //tricklingDownPart += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - now).count();

    if (peaks.empty()) {
        for (auto& [allegedMinDep, info] : allegedMinDeps->entrySet()) {
            if (!info->isExtremal_ && !globalVisitees_->containsKey(allegedMinDep)) {
                LOG(DEBUG) << boost::format{"[%1%] Minimum dependency: %2% (error=%3%)"}
                    % recursionDepth_ % allegedMinDep.toString() % info->error_;
                // TODO: тут надо сделать non-const - костыльный mutable; опять Info в двух местах хранится
                info->isExtremal_ = true;
                globalVisitees_->put(allegedMinDep, std::make_unique<VerticalInfo>(*info));
                strategy_->registerDependency(allegedMinDep, info->error_, *context_);
            }
        }
        tricklingDown += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - now).count();
    } else {
        LOG(DEBUG) << boost::format{"* %1% new peaks (%2%)"} % peaks.size() % "UNIMPLEMENTED";
        auto newScope = std::make_unique<util::VerticalMap<Vertical>>(context_->getSchema());
        std::sort_heap(peaks.begin(), peaks.end(), peaksComparator);
        for (auto& peak : peaks) {
            newScope->put(peak.vertical_, std::make_unique<Vertical>(peak.vertical_));
        }

        double newSampleBoost = sampleBoost_ * sampleBoost_;
        LOG(DEBUG) << boost::format{"* Increasing sampling boost factor to %1%"} % newSampleBoost;

        auto scopeVerticals  = newScope->keySet();
        // TODO: что делать с strategy, globalVisitees?
        auto nestedSearchSpace = std::make_unique<SearchSpace> (
                -1, strategy_->createClone(), std::move(newScope), std::move(globalVisitees_),
                context_->getSchema(), launchPads_.key_comp(), recursionDepth_ + 1,
                sampleBoost_ * context_->getConfiguration().sampleBooster
                );
        nestedSearchSpace->setContext(context_);

        std::unordered_set<Column> scopeColumns;
        for (auto& vertical : scopeVerticals) {
            for (auto column : vertical.getColumns()) {
                scopeColumns.insert(*column);
            }
        }
        for (auto& scopeColumn : scopeColumns) {
            LOG(TRACE) << "createDependencyCandidate while building a nested search space";
            // TODO: again problems with conversion: Column* -> Vertical*. If this is too inefficient, consider refactoring
            nestedSearchSpace->addLaunchPad(strategy_->createDependencyCandidate(
                    static_cast<Vertical>(scopeColumn)
                    ));
        }
        auto prev = std::chrono::system_clock::now();
        numNested++;
        //std::cout << static_cast<std::string>(*strategy_) << ' ';
        //std::cout << numNested << std::endl;
        tricklingDown += std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::system_clock::now() - now).count();
        nestedSearchSpace->moveInLocalVisitees(std::move(localVisitees_));
        nestedSearchSpace->discover();
        tricklingDownPart += std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::system_clock::now() - prev).count();
        globalVisitees_ = nestedSearchSpace->moveOutGlobalVisitees();
        localVisitees_ = nestedSearchSpace->moveOutLocalVisitees();

        for (auto& [allegedMinDep, info] : allegedMinDeps->entrySet()) {
            if (!isImpliedByMinDep(allegedMinDep, globalVisitees_.get())) {
                LOG(DEBUG) << boost::format{"[%1%] Minimum dependency: %2% (error=%3%) (was right after all)"}
                    % recursionDepth_ % allegedMinDep.toString() % info->error_;
                // TODO: тут надо сделать non-const - костыльный mutable; опять Info в двух местах хранится
                info->isExtremal_ = true;
                globalVisitees_->put(allegedMinDep, std::make_unique<VerticalInfo>(*info));
                strategy_->registerDependency(allegedMinDep, info->error_, *context_);
            }
        }
    }
}

std::optional<Vertical> SearchSpace::trickleDownFrom(
        DependencyCandidate minDepCandidate, DependencyStrategy* strategy,
        util::VerticalMap<VerticalInfo>* allegedMinDeps,
        std::unordered_set<Vertical> & allegedNonDeps,
        util::VerticalMap<VerticalInfo>* globalVisitees, double boostFactor) {
    auto now = std::chrono::system_clock::now();
    if (minDepCandidate.error_.getMin() > strategy->maxDependencyError_) {
        throw std::runtime_error("Error in trickleDownFrom: minDepCandidate's error should be <= maxError");
    }

    bool areAllParentsKnownNonDeps = true;
    if (minDepCandidate.vertical_.getArity() > 1) {
        std::priority_queue<DependencyCandidate, std::vector<DependencyCandidate>,
            std::function<bool (DependencyCandidate&, DependencyCandidate&)>>
            parentCandidates([](auto& candidate1, auto& candidate2)
                {return DependencyCandidate::minErrorComparator(candidate1, candidate2); });
        for (auto& parentVertical : minDepCandidate.vertical_.getParents()) {
            if (isKnownNonDependency(parentVertical, localVisitees_.get())
                    || isKnownNonDependency(parentVertical, globalVisitees))
                continue;
            if (allegedNonDeps.count(parentVertical) != 0) {
                areAllParentsKnownNonDeps = false;
                continue;
            }
            // TODO: construction methods should return unique_ptr<...>
            LOG(TRACE) << "createDependencyCandidate while trickling down from";
            parentCandidates.push(strategy->createDependencyCandidate(parentVertical));
        }

        while (!parentCandidates.empty()) {
            auto parentCandidate = parentCandidates.top();
            parentCandidates.pop();

            if (parentCandidate.error_.getMin() > strategy->minNonDependencyError_) {
                do {
                    if (parentCandidate.isExact()) {
                        localVisitees_->put(parentCandidate.vertical_,
                                std::make_unique<VerticalInfo>(VerticalInfo::forNonDependency()));
                    } else {
                        allegedNonDeps.insert(parentCandidate.vertical_);
                        areAllParentsKnownNonDeps = false;
                    }
                    if (!parentCandidates.empty()) {
                        parentCandidate = parentCandidates.top();
                        parentCandidates.pop();
                    }
                } while (!parentCandidates.empty());
                break;
            }

            tricklingDownFrom += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::system_clock::now() - now).count();

            auto allegedMinDep = trickleDownFrom(
                    std::move(parentCandidate),
                    strategy,
                    allegedMinDeps,
                    allegedNonDeps,
                    globalVisitees,
                    boostFactor
                    );

            now = std::chrono::system_clock::now();

            if (allegedMinDep.has_value()){
                return allegedMinDep;
            }

            if (!minDepCandidate.isExact()) {
                double error = strategy->calculateError(minDepCandidate.vertical_);
                // TODO: careful with reference shenanigans - looks like it works this way in the original
                minDepCandidate = DependencyCandidate(minDepCandidate.vertical_,
                                                      util::ConfidenceInterval(error), true);
                if (error > strategy->minNonDependencyError_) break;
            }
        }
    }

    double candidateError = minDepCandidate.isExact()
        ? minDepCandidate.error_.get()
        : strategy->calculateError(minDepCandidate.vertical_);
    [[maybe_unused]] double errorDiff = candidateError - minDepCandidate.error_.getMean();
    if (candidateError <= strategy->maxDependencyError_) {
        LOG(TRACE) << boost::format{"* Found %1%-ary minimum dependency candidate: %2%"}
            % minDepCandidate.vertical_.getArity() % minDepCandidate;
        allegedMinDeps->removeSupersetEntries(minDepCandidate.vertical_);
        allegedMinDeps->put(minDepCandidate.vertical_,
                std::make_unique<VerticalInfo>(true, areAllParentsKnownNonDeps, candidateError));
        if (areAllParentsKnownNonDeps && context_->getConfiguration().isCheckEstimates) {
            requireMinimalDependency(strategy, minDepCandidate.vertical_);
        }
        tricklingDownFrom += std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::system_clock::now() - now).count();
        return minDepCandidate.vertical_;
    } else {
        LOG(TRACE) << boost::format{"* Guessed incorrect %1%-ary minimum dependency candidate."}
            % minDepCandidate.vertical_.getArity();
        localVisitees_->put(minDepCandidate.vertical_, std::make_unique<VerticalInfo>(VerticalInfo::forNonDependency()));

        if (strategy->shouldResample(minDepCandidate.vertical_, boostFactor)) {
            context_->createFocusedSample(minDepCandidate.vertical_, boostFactor);
        }
        tricklingDownFrom += std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::system_clock::now() - now).count();
        return std::optional<Vertical>();
    }
}

void SearchSpace::requireMinimalDependency(DependencyStrategy* strategy, Vertical const& minDependency) {
    double error = strategy->calculateError(minDependency);
    if (error > strategy->maxDependencyError_) {
        throw std::runtime_error("Wrong minimal dependency estimate");
    }
    if (minDependency.getArity() > 1) {
        for (auto& parent : minDependency.getParents()) {
            double parentError = strategy->calculateError(parent);
            if (parentError <= strategy->minNonDependencyError_) {
                throw std::runtime_error("Wrong minimal dependency estimate");
            }
        }
    }
}

std::vector<Vertical> SearchSpace::getSubsetDeps(
        Vertical const& vertical, util::VerticalMap<VerticalInfo>* verticalInfos) {

    auto subsetEntries = verticalInfos->getSubsetEntries(vertical);
    auto subsetEntriesEnd = std::remove_if(subsetEntries.begin(), subsetEntries.end(),
            [](auto& entry){ return !entry.second->isDependency_;});
    std::vector<Vertical> subsetDeps;

    std::transform(subsetEntries.begin(), subsetEntriesEnd, std::inserter(subsetDeps, subsetDeps.begin()),
            [](auto& entry){
        return entry.first;
    });

    return subsetDeps;
}

bool SearchSpace::isImpliedByMinDep(Vertical const& vertical,
                                    util::VerticalMap<VerticalInfo>* verticalInfos) {
    // TODO: function<bool(Vertical, ...)> --> function<bool(Vertical&, ...)>
    return verticalInfos->getAnySubsetEntry(
            vertical, []([[maybe_unused]] auto vertical, auto info) {
                return info->isDependency_ && info->isExtremal_;
            }).second != nullptr;
}

bool SearchSpace::isKnownNonDependency(Vertical const& vertical,
                                       util::VerticalMap<VerticalInfo>* verticalInfos) {
    return verticalInfos->getAnySupersetEntry(
            vertical, []([[maybe_unused]] auto vertical, auto info) {
                return !info->isDependency_;
            }).second != nullptr;
}

void SearchSpace::printStats() const {
    using std::cout, std::endl;
    cout << "Trickling down from: " << tricklingDownFrom / 1000000 << endl;
    cout << "Trickling down: " << tricklingDown / 1000000 - tricklingDownFrom / 1000000 << endl;
    cout << "Trickling down nested:" << tricklingDownPart / 1000000 << endl;
    cout << "Num nested: " << numNested / 1000000 << endl;
    cout << "Ascending: " << ascending / 1000000<< endl;
    cout << "Polling: " << pollingLaunchPads / 1000000 << endl;
    cout << "Returning launch pad: " << returningLaunchPad / 1000000 << endl;
}

void SearchSpace::ensureInitialized() {
    strategy_->ensureInitialized(this);
    std::string initializedLaunchPads;
    for (auto const& pad : launchPads_) {
        initializedLaunchPads += std::string(pad) + " ";
    }
    LOG(TRACE) << "Initialized with launch pads: " + initializedLaunchPads;
}

