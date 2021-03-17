#include "logging/easylogging++.h"
#include "SearchSpace.h"
#include <queue>
/*SearchSpace::SearchSpace(int id, std::shared_ptr<DependencyStrategy> strategy,
                         std::unique_ptr<VerticalMap<Vertical>> scope, VerticalMap<VerticalInfo> globalVisitees,
                         std::shared_ptr<RelationalSchema> schema, std::function<bool(DependencyCandidate const &,
                                                                                      DependencyCandidate const &)> const &dependencyCandidateComparator,
                         int recursionDepth, double sampleBoost) : id_(id), strategy_(strategy), scope_(std::move(scope)), globalVisitees_(globalVisitees), recursionDepth_(recursionDepth),
                                                                   sampleBoost_(sampleBoost), launchPadIndex_(schema), launchPads_(dependencyCandidateComparator) {}*/

// TODO: extra careful with const& -> shared_ptr conversions via make_shared-smart pointer may delete the object - pass empty deleter [](*) {}

// TODO: consider storing only containers of shared_ptrs

void SearchSpace::discover(std::unique_ptr<VerticalMap<VerticalInfo>> localVisitees) {
    std::cout << "Discovering in: " << static_cast<std::string>(*strategy_) << std::endl;
    while (true) {  // на второй итерации дропается
        auto now = std::chrono::system_clock::now();
        std::optional<DependencyCandidate> launchPad = pollLaunchPad(localVisitees.get());
        if (!launchPad.has_value()) break;

        if (localVisitees == nullptr)
            localVisitees = std::make_unique<VerticalMap<VerticalInfo>>(context_->getSchema());

        bool isDependencyFound = ascend(*launchPad, localVisitees.get());
        pollingLaunchPads += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - now).count();
        returnLaunchPad(*launchPad, !isDependencyFound);
    }
}

std::optional<DependencyCandidate> SearchSpace::pollLaunchPad(
        VerticalMap<VerticalInfo>* localVisitees) {
    while (true) {
        if (launchPads_.empty()) {
            if (deferredLaunchPads_.empty()) return std::optional<DependencyCandidate>();

            launchPads_.insert(deferredLaunchPads_.begin(), deferredLaunchPads_.end());
            deferredLaunchPads_.clear();
        }

        auto launchPad = launchPads_.extract(launchPads_.begin()).value();

        launchPads_.erase(launchPads_.begin());
        launchPadIndex_->remove(launchPad.vertical_);

        if (isImpliedByMinDep(launchPad.vertical_, globalVisitees_.get())
            || (localVisitees != nullptr && isImpliedByMinDep(launchPad.vertical_, localVisitees))) {
            launchPadIndex_->remove(launchPad.vertical_);
            LOG(TRACE) << "* Removing subset-pruned launch pad {" << launchPad.vertical_.toString() << '}';
            continue;
        }

        auto supersetEntries = globalVisitees_->getSupersetEntries(launchPad.vertical_);
        if (localVisitees != nullptr) {
            auto localSupersetEntries = localVisitees->getSupersetEntries(launchPad.vertical_);
            auto endIterator = std::remove_if(localSupersetEntries.begin(), localSupersetEntries.end(), [](auto entry) { return !entry.second->isPruningSubsets(); });

            std::for_each(localSupersetEntries.begin(), endIterator, [&supersetEntries](auto entry) { supersetEntries.push_back(entry); });
        }
        if (supersetEntries.empty()) return launchPad;
        LOG(TRACE) << boost::format{"* Escaping launchPad %1% from: %2%"} % launchPad.vertical_.toString() % "[UNIMPLEMENTED]";
        std::vector<Vertical> supersetVerticals;

        for (auto& entry : supersetEntries) {
            supersetVerticals.push_back(entry.first);
        }

        escapeLaunchPad(launchPad.vertical_, std::move(supersetVerticals), localVisitees);
    }
}

// this move looks legit IMO
void SearchSpace::escapeLaunchPad(Vertical const& launchPad,
                                  std::vector<Vertical> pruningSupersets,
                                  VerticalMap<VerticalInfo>* localVisitees) {
    // TODO: list<не shared_ptr>, чтобы не создавать новые поинтеры? Мб создать новый лист list<Vertical>

    std::transform(pruningSupersets.begin(), pruningSupersets.end(), pruningSupersets.begin(),
            [this](auto superset) { return superset->invert()->without(strategy_->getIrrelevantColumns()); } );

    std::function<bool (Vertical const&)> pruningFunction = [this, &launchPad, &localVisitees] (auto const& hittingSetCandidate) -> bool {
        if (scope_ != nullptr && scope_->getAnySupersetEntry(hittingSetCandidate).second == nullptr) {
            return true;
        }

        auto launchPadCandidate_ptr = launchPad.Union(hittingSetCandidate);

        if ((localVisitees == nullptr && isImpliedByMinDep(launchPadCandidate_ptr, localVisitees))
            || isImpliedByMinDep(launchPadCandidate_ptr, globalVisitees_)) {
            return true;
        }

        if (launchPadIndex_->getAnySubsetEntry(*launchPadCandidate_ptr).second != nullptr) {
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
        LOG(TRACE) << boost::format{"Escaping %1% pruned by %2%"} % launchPad->toString() % pruningSupersetsStr;
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
            % escapedLaunchPad.vertical_.getArity() % context_->configuration_.maxLHS;
        if (escapedLaunchPad.vertical_.getArity() <= context_->configuration_.maxLHS) {
            launchPads_.insert(escapedLaunchPad);
            launchPadIndex_->put(escapedLaunchPad.vertical_, std::make_unique<DependencyCandidate>(escapedLaunchPad));
        }
    }
}

void SearchSpace::addLaunchPad(const DependencyCandidate &launchPad) {
    launchPads_.insert(launchPad);
    launchPadIndex_.put(*launchPad.vertical_, std::make_shared<DependencyCandidate>(launchPad));
}

void SearchSpace::returnLaunchPad(DependencyCandidate const &launchPad, bool isDefer) {
    if (isDefer && context_->configuration_.isDeferFailedLaunchPads) {
        deferredLaunchPads_.push_back(launchPad);
        LOG(TRACE) << boost::format{"Deferred seed %1%"} % launchPad.vertical_->toString();
    }
    else {
        launchPads_.insert(launchPad);
    }
    launchPadIndex_.put(*launchPad.vertical_, std::make_shared<DependencyCandidate>(launchPad));
}

bool SearchSpace::ascend(DependencyCandidate const &launchPad,
                         std::shared_ptr<VerticalMap<std::shared_ptr<VerticalInfo>>> localVisitees) {

    auto now = std::chrono::system_clock::now();

    LOG(DEBUG) << boost::format{"===== Ascending from %1% ======"} % strategy_->format(launchPad.vertical_);

    if (strategy_->shouldResample(launchPad.vertical_, sampleBoost_)) {
        LOG(TRACE) << "Resampling.";
        context_->createFocusedSample(launchPad.vertical_, sampleBoost_);
    }

    DependencyCandidate traversalCandidate = launchPad;
    boost::optional<double> error;

    while (true) {
        LOG(TRACE) << boost::format{"-> %1%"} % traversalCandidate.vertical_->toString();

        if (context_->configuration_.isCheckEstimates) {
            checkEstimate(strategy_, traversalCandidate);
        }

        if (traversalCandidate.isExact()) {
            // TODO: getError()?
            error = traversalCandidate.error_.get();

            bool canBeDependency = *error <= strategy_->maxDependencyError_;
            localVisitees->put(*traversalCandidate.vertical_, std::make_shared<VerticalInfo>(
                    canBeDependency,
                    false,
                    *error
                    ));
            if (canBeDependency) break;
        }
        else {
            if (traversalCandidate.error_.getMin() > strategy_->maxDependencyError_) {
                LOG(TRACE) << boost::format {"  Skipping check form %1% (estimated error: %2%)."}
                    % traversalCandidate.vertical_->toString() % traversalCandidate.error_;
                error.reset();
            }
            else {
                error = context_->configuration_.isEstimateOnly
                        ? traversalCandidate.error_.getMean()
                        : strategy_->calculateError(traversalCandidate.vertical_);
                // double errorDiff = *error - traversalCandidate.error_.getMean();

                localVisitees->put(*traversalCandidate.vertical_, std::make_shared<VerticalInfo>(
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

        if (traversalCandidate.vertical_->getArity() >=
            context_->relationData_->getNumColumns() - strategy_->getNumIrrelevantColumns()
            || traversalCandidate.vertical_->getArity() >= context_->configuration_.maxLHS) {
            break;
        }

        boost::optional<DependencyCandidate> nextCandidate;
        int numSeenElements = isAscendRandomly_ ? 1 : -1;
        for (auto& extensionColumn : context_->getSchema()->getColumns()) {
            if (traversalCandidate.vertical_->getColumnIndices()[extensionColumn->getIndex()]
                    || strategy_->isIrrelevantColumn(extensionColumn)) {
                continue;
            }
            auto extendedVertical = traversalCandidate.vertical_->Union(static_cast<Vertical>(*extensionColumn));

            if (scope_ != nullptr && scope_->getSupersetEntries(*extendedVertical).empty()) {
                continue;
            }

            bool isSubsetPruned = isImpliedByMinDep(extendedVertical, globalVisitees_);
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
        LOG(TRACE) << boost::format{"  Hit ceiling at %1%."} % traversalCandidate.vertical_->toString();
        error = strategy_->calculateError(traversalCandidate.vertical_);
        double errorDiff = *error - traversalCandidate.error_.getMean();
        LOG(TRACE) << boost::format{"  Checking candidate... actual error: %1%"} % *error;
    }
    ascending += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - now).count();

    if (*error <= strategy_->maxDependencyError_) {
        LOG(TRACE) << boost::format{"  Key peak in climbing phase e(%1%)=%2% -> Need to minimize."}
            % traversalCandidate.vertical_->toString() % *error;
        trickleDown(traversalCandidate.vertical_, *error, localVisitees);

        if (recursionDepth_ == 0) {
            assert(scope_ == nullptr);
            globalVisitees_->put(*traversalCandidate.vertical_, std::make_shared<VerticalInfo>(VerticalInfo::forMinimalDependency()));
        }

        return true;
    } else {
        if (recursionDepth_ == 0) {
            assert(scope_ == nullptr);
            globalVisitees_->put(*traversalCandidate.vertical_, std::make_shared<VerticalInfo>(VerticalInfo::forMaximalNonDependency()));
            LOG(DEBUG) << boost::format{"[---] %1% is maximum non-dependency (err=%2%)."}
                % traversalCandidate.vertical_->toString() % *error;
        } else {
            localVisitees->put(*traversalCandidate.vertical_, std::make_shared<VerticalInfo>(VerticalInfo::forNonDependency()));
            LOG(DEBUG) << boost::format{"      %1% is local-maximum non-dependency (err=%2%)."}
                % traversalCandidate.vertical_->toString() % *error;
        }
    }
    return false;
}

void SearchSpace::checkEstimate(std::shared_ptr<DependencyStrategy> strategy,
                                DependencyCandidate const &traversalCandidate) {
    std::cout << "Stepped into method 'checkEstimate' - not implemented yet being a debug method\n";
}

void SearchSpace::trickleDown(std::shared_ptr<Vertical> mainPeak, double mainPeakError,
                              std::shared_ptr<VerticalMap<std::shared_ptr<VerticalInfo>>> localVisitees) {
    LOG(DEBUG) << boost::format{"====== Trickling down from %1% ======"} % mainPeak->toString();

    std::unordered_set<Vertical> maximalNonDeps;
    auto allegedMinDeps = std::make_shared<VerticalMap<std::shared_ptr<VerticalInfo>>>(context_->getSchema());
    // arityComparator returns true if candidate1 < candidate2. We need a candidate with the smallest arity to be the highest value in the heap,
    // i.e. such a candidate that arityComparator returns false for every other candidate => peaksComparator should return !(candidate1 < candidate2)
    // In one word, Java priority queue highlights the least element, C++ - the largest
    std::function<bool(std::shared_ptr<DependencyCandidate>, std::shared_ptr<DependencyCandidate>)> peaksComparator =
            [](auto candidate1, auto candidate2) -> bool { return DependencyCandidate::arityComparator(*candidate1, *candidate2); };
    std::vector<std::shared_ptr<DependencyCandidate>> peaks;
    std::make_heap(peaks.begin(), peaks.end(), peaksComparator);
    peaks.push_back(std::make_shared<DependencyCandidate>(mainPeak, ConfidenceInterval(mainPeakError), true));
    std::push_heap(peaks.begin(), peaks.end(), peaksComparator);
    std::unordered_set<Vertical> allegedNonDeps;

    auto now = std::chrono::system_clock::now();

    while (!peaks.empty()) {
        // пока ситуация такая: при спуске из АС правильно считает C за не-ФЗ, АС за ФЗ, но не выходит из этого цикла,
        // т.к. не убирает АС: subsetDeps - пусто (должно ли оно быть пустым?)
        auto peak = peaks.front();

        auto subsetDeps = getSubsetDeps(peak->vertical_, allegedMinDeps);
        if (!subsetDeps.empty()) {
            std::pop_heap(peaks.begin(), peaks.end(), peaksComparator);
            peaks.pop_back();

            auto peakHittingSet = context_->getSchema()->calculateHittingSet(
                    std::move(subsetDeps), boost::optional<std::function<bool (Vertical const&)>>());
            std::unordered_set<std::shared_ptr<Vertical>> escapedPeakVerticals;

            for (auto vertical : peakHittingSet) {
                escapedPeakVerticals.insert(peak->vertical_->without(*vertical));
            }
//            std::transform(peakHittingSet.begin(), peakHittingSet.end(), escapedPeakVerticals.begin(),
//                [&peak](auto vertical) { return std::make_shared<Vertical>(peak->vertical_->without(*vertical)); });

            for (auto& escapedPeakVertical : escapedPeakVerticals) {
                if (escapedPeakVertical->getArity() > 0
                    && allegedNonDeps.find(*escapedPeakVertical) == allegedNonDeps.end()) {
                    // TODO: escapedPeakVertical, [](Vertical* v) {} ?
                    //auto escapedPeakVertical_ptr = std::make_shared<Vertical>(escapedPeakVertical);
                    LOG(TRACE) << "createDependencyCandidate as an escaped peak while trickling down";
                    auto escapedPeak = strategy_->createDependencyCandidate(escapedPeakVertical);

                    if (escapedPeak.error_.getMean() > strategy_->minNonDependencyError_) {
                        allegedNonDeps.insert(*escapedPeakVertical);
                        continue;
                    }
                    if (isKnownNonDependency(escapedPeakVertical, localVisitees)
                        || isKnownNonDependency(escapedPeakVertical, globalVisitees_)) {
                        continue;
                    }
                    peaks.push_back(std::make_shared<DependencyCandidate>(escapedPeak));
                    std::push_heap(peaks.begin(), peaks.end(), peaksComparator);
                }
            }
            continue;
        }
        auto allegedMinDep = trickleDownFrom(*peak, strategy_,
                                             allegedMinDeps, allegedNonDeps, localVisitees, globalVisitees_, sampleBoost_);
        if (allegedMinDep == nullptr) {
            std::pop_heap(peaks.begin(), peaks.end(), peaksComparator);
            peaks.pop_back();
        }
    }

    LOG(DEBUG) << boost::format{"* %1% alleged minimum dependencies (%2%)"} % allegedMinDeps->getSize() % "UNIMPLEMENTED";

    int numUncertainMinDeps = 0;
    for (auto& [allegedMinDep, info] : allegedMinDeps->entrySet()) {
        if (info->isExtremal_ && !globalVisitees_->containsKey(*allegedMinDep)) {
            LOG(DEBUG) << boost::format{"[%1%] Minimum dependency: %2% (error=%3%)"}
                % recursionDepth_ % allegedMinDep->toString() % info->error_;
            globalVisitees_->put(*allegedMinDep, info);
            strategy_->registerDependency(allegedMinDep, info->error_, *context_);  // TODO: is it expensive to dereference ptr like this?
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
            std::list<std::shared_ptr<Vertical>>(allegedMinDepsSet.begin(), allegedMinDepsSet.end()),
            boost::optional<std::function<bool (Vertical const&)>>());
    std::unordered_set<std::shared_ptr<Vertical>> allegedMaxNonDeps;


    for (auto minLeaveOutVertical : allegedMaxNonDepsHS) {
        allegedMaxNonDeps.insert(minLeaveOutVertical->invert(*mainPeak));
    }

    LOG(DEBUG) << boost::format{"* %1% alleged maximum non-dependencies (%2%)"}
        % allegedMaxNonDeps.size() % "UNIMPLEMENTED";
    //std::transform(allegedMaxNonDepsHS.begin(), allegedMaxNonDepsHS.end(), allegedMaxNonDeps.begin(),
    //        [mainPeak](auto minLeaveOutVertical) { return minLeaveOutVertical->invert(*mainPeak); });

    // checking the consistency of all data structures
    if (auto allegedMinDepsKeySet = allegedMinDeps->keySet();
            !std::all_of(allegedMinDepsKeySet.begin(), allegedMinDepsKeySet.end(),
                [mainPeak](auto vertical) -> bool { return mainPeak->contains(*vertical); })) {
        throw std::runtime_error("Main peak should contain all alleged min dependencies");
    }
    if (!std::all_of(allegedMaxNonDeps.begin(), allegedMaxNonDeps.end(),
                     [mainPeak](auto vertical) -> bool { return mainPeak->contains(*vertical); })) {
        throw std::runtime_error("Main peak should contain all alleged max non-dependencies");
    }

    for (auto allegedMaxNonDep : allegedMaxNonDeps) {
        if (allegedMaxNonDep->getArity() == 0) continue;

        if (maximalNonDeps.find(*allegedMaxNonDep) != maximalNonDeps.end()
                || isKnownNonDependency(allegedMaxNonDep, localVisitees)
                || isKnownNonDependency(allegedMaxNonDep, globalVisitees_)) {
            continue;
        }

        double error = context_->configuration_.isEstimateOnly
                ? strategy_->createDependencyCandidate(allegedMaxNonDep).error_.getMean()
                : strategy_->calculateError(allegedMaxNonDep);
        bool isNonDep = error > strategy_->minNonDependencyError_;
        LOG(TRACE) << boost::format {"* Alleged maximal non-dependency %1%: non-dep?: %2%, error: %3%"}
            % allegedMaxNonDep->toString() % isNonDep % error;
        if (isNonDep) {
            maximalNonDeps.insert(*allegedMaxNonDep);
            localVisitees->put(*allegedMaxNonDep, std::make_shared<VerticalInfo>(VerticalInfo::forNonDependency()));
        } else {
            peaks.push_back(std::make_shared<DependencyCandidate>(allegedMaxNonDep, ConfidenceInterval(error), true));
            std::push_heap(peaks.begin(), peaks.end(), peaksComparator);
        }
    }
    //tricklingDownPart += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - now).count();

    if (peaks.empty()) {
        for (auto& [allegedMinDep, info] : allegedMinDeps->entrySet()) {
            if (!info->isExtremal_ && !globalVisitees_->containsKey(*allegedMinDep)) {
                LOG(DEBUG) << boost::format{"[%1%] Minimum dependency: %2% (error=%3%)"}
                    % recursionDepth_ % allegedMinDep->toString() % info->error_;
                info->isExtremal_ = true;
                globalVisitees_->put(*allegedMinDep, info);
                strategy_->registerDependency(allegedMinDep, info->error_, *context_);
            }
        }
        tricklingDown += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - now).count();
    } else {
        LOG(DEBUG) << boost::format{"* %1% new peaks (%2%)"} % peaks.size() % "UNIMPLEMENTED";
        auto newScope = std::make_unique<VerticalMap<std::shared_ptr<Vertical>>>(context_->getSchema());
        std::sort_heap(peaks.begin(), peaks.end(), peaksComparator);
        for (auto& peak : peaks) {
            newScope->put(*peak->vertical_, peak->vertical_);
        }

        double newSampleBoost = sampleBoost_ * sampleBoost_;
        LOG(DEBUG) << boost::format{"* Increasing sampling boost factor to %1%"} % newSampleBoost;

        // Nested слишком часто создаётся - что-то не так
        auto scopeVerticals  = newScope->keySet();
        auto nestedSearchSpace = std::make_shared<SearchSpace> (
                -1, strategy_, std::move(newScope), globalVisitees_, context_->getSchema(),
                launchPads_.key_comp(), recursionDepth_ + 1,
                sampleBoost_ * context_->configuration_.sampleBooster
                );
        nestedSearchSpace->setContext(context_);

        std::unordered_set<std::shared_ptr<Column>> scopeColumns;
        for (auto vertical : scopeVerticals) {
            for (auto column : vertical->getColumns()) {
                scopeColumns.insert(column);
            }
        }
        for (auto scopeColumn : scopeColumns) {
            LOG(TRACE) << "createDependencyCandidate while building a nested search space";
            // TODO: again problems with conversion: Column* -> Vertical*. If this is too inefficient, consider refactoring
            nestedSearchSpace->addLaunchPad(strategy_->createDependencyCandidate(
                    std::make_shared<Vertical>(static_cast<Vertical>(*scopeColumn))
                    ));
        }
        auto prev = std::chrono::system_clock::now();
        numNested++;
        //std::cout << static_cast<std::string>(*strategy_) << ' ';
        //std::cout << numNested << std::endl;
        tricklingDown += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - now).count();
        nestedSearchSpace->discover(localVisitees);
        tricklingDownPart += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - prev).count();

        for (auto& [allegedMinDep, info] : allegedMinDeps->entrySet()) {
            if (!isImpliedByMinDep(allegedMinDep, globalVisitees_)) {
                LOG(DEBUG) << boost::format{"[%1%] Minimum dependency: %2% (error=%3%) (was right after all)"}
                    % recursionDepth_ % allegedMinDep->toString() % info->error_;
                info->isExtremal_ = true;
                globalVisitees_->put(*allegedMinDep, info);
                strategy_->registerDependency(allegedMinDep, info->error_, *context_);
            }
        }
    }
}

std::shared_ptr<Vertical>
SearchSpace::trickleDownFrom(DependencyCandidate &minDepCandidate, std::shared_ptr<DependencyStrategy> strategy,
                             std::shared_ptr<VerticalMap<std::shared_ptr<VerticalInfo>>> allegedMinDeps,
                             std::unordered_set<Vertical> &allegedNonDeps,
                             std::shared_ptr<VerticalMap<std::shared_ptr<VerticalInfo>>> localVisitees,
                             std::shared_ptr<VerticalMap<std::shared_ptr<VerticalInfo>>> globalVisitees,
                             double boostFactor) {
    auto now = std::chrono::system_clock::now();
    if (minDepCandidate.error_.getMin() > strategy->maxDependencyError_)
        throw std::runtime_error("Error in trickleDownFrom: minDepCandidate's error should be <= maxError");

    bool areAllParentsKnownNonDeps = true;
    if (minDepCandidate.vertical_->getArity() > 1) {
        std::priority_queue<std::shared_ptr<DependencyCandidate>, std::vector<std::shared_ptr<DependencyCandidate>>,
            std::function<bool (std::shared_ptr<DependencyCandidate>, std::shared_ptr<DependencyCandidate>)>>
            parentCandidates([](auto candidate1, auto candidate2)
                {return DependencyCandidate::minErrorComparator(*candidate1, *candidate2); });
        for (auto parentVertical : minDepCandidate.vertical_->getParents()) {
            if (isKnownNonDependency(parentVertical, localVisitees)
                    || isKnownNonDependency(parentVertical, globalVisitees))
                continue;
            if (allegedNonDeps.count(*parentVertical) != 0) {
                areAllParentsKnownNonDeps = false;
                continue;
            }
            // TODO: construction methods should return unique_ptr<...>
            LOG(TRACE) << "createDependencyCandidate while trickling down from";
            parentCandidates.push(std::make_shared<DependencyCandidate>(strategy->createDependencyCandidate(parentVertical)));
        }

        while (!parentCandidates.empty()) {
            auto parentCandidate = parentCandidates.top();
            parentCandidates.pop();

            if (parentCandidate->error_.getMin() > strategy->minNonDependencyError_) {
                do {
                    if (parentCandidate->isExact()) {
                        localVisitees->put(*parentCandidate->vertical_,
                                std::make_shared<VerticalInfo>(VerticalInfo::forNonDependency()));
                    } else {
                        allegedNonDeps.insert(*parentCandidate->vertical_);
                        areAllParentsKnownNonDeps = false;
                    }
                    if (!parentCandidates.empty()) {
                        parentCandidate = parentCandidates.top();
                        parentCandidates.pop();
                    }
                } while (!parentCandidates.empty());
                break;
            }

            tricklingDownFrom += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - now).count();

            auto allegedMinDep = trickleDownFrom(
                    *parentCandidate,
                    strategy,
                    allegedMinDeps,
                    allegedNonDeps,
                    localVisitees,
                    globalVisitees,
                    boostFactor
                    );

            now = std::chrono::system_clock::now();

            if (allegedMinDep != nullptr){
                return allegedMinDep;
            }

            if (!minDepCandidate.isExact()) {
                double error = strategy->calculateError(minDepCandidate.vertical_);
                // TODO: careful with reference shenanigans - looks like it works this way in the original
                minDepCandidate = DependencyCandidate(minDepCandidate.vertical_, ConfidenceInterval(error), true);
                if (error > strategy->minNonDependencyError_) break;
            }
        }
    }

    double candidateError = minDepCandidate.isExact()
        ? minDepCandidate.error_.get()
        : strategy->calculateError(minDepCandidate.vertical_);
    double errorDiff = candidateError - minDepCandidate.error_.getMean();
    if (candidateError <= strategy->maxDependencyError_) {
        LOG(TRACE) << boost::format{"* Found %1%-ary minimum dependency candidate: %2%"}
            % minDepCandidate.vertical_->getArity() % minDepCandidate;
        allegedMinDeps->removeSupersetEntries(*minDepCandidate.vertical_);
        allegedMinDeps->put(*minDepCandidate.vertical_,
                std::make_shared<VerticalInfo>(true, areAllParentsKnownNonDeps, candidateError));
        if (areAllParentsKnownNonDeps && context_->configuration_.isCheckEstimates) {
            requireMinimalDependency(strategy, minDepCandidate.vertical_);
        }
        tricklingDownFrom += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - now).count();
        return minDepCandidate.vertical_;
    } else {
        LOG(TRACE) << boost::format{"* Guessed incorrect %1%-ary minimum dependency candidate."}
            % minDepCandidate.vertical_->getArity();
        localVisitees->put(*minDepCandidate.vertical_, std::make_shared<VerticalInfo>(VerticalInfo::forNonDependency()));

        if (strategy->shouldResample(minDepCandidate.vertical_, boostFactor)) {
            context_->createFocusedSample(minDepCandidate.vertical_, boostFactor);
        }
        tricklingDownFrom += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - now).count();
        return nullptr;
    }
}

void SearchSpace::requireMinimalDependency(std::shared_ptr<DependencyStrategy> strategy,
                                           std::shared_ptr<Vertical> minDependency) {
    double error = strategy->calculateError(minDependency);
    if (error > strategy->maxDependencyError_) {
        throw std::runtime_error("Wrong minimal dependency estimate");
    }
    if (minDependency->getArity() > 1) {
        for (auto parent : minDependency->getParents()) {
            double parentError = strategy->calculateError(parent);
            if (parentError <= strategy->minNonDependencyError_) {
                throw std::runtime_error("Wrong minimal dependency estimate");
            }
        }
    }
}

std::list<std::shared_ptr<Vertical>> SearchSpace::getSubsetDeps(
        Vertical const& vertical, VerticalMap<VerticalInfo>* verticalInfos) {

    auto subsetEntries = verticalInfos->getSubsetEntries(vertical);
    auto subsetEntriesEnd = std::remove_if(subsetEntries.begin(), subsetEntries.end(),
            [](auto& entry){ return !entry.second->isDependency_;});
    std::list<std::shared_ptr<Vertical>> subsetDeps;

    std::transform(subsetEntries.begin(), subsetEntriesEnd, std::inserter(subsetDeps, subsetDeps.begin()),
            [](auto& entry){
        return entry.first;
    });

    return subsetDeps;
}

bool SearchSpace::isImpliedByMinDep(Vertical const& vertical, VerticalMap<VerticalInfo>* verticalInfos) {
    // TODO: function<bool(Vertical, ...)> --> function<bool(Vertical&, ...)>
    return verticalInfos->getAnySubsetEntry(
            vertical,[](auto vertical, auto info) -> bool { return info->isDependency_ && info->isExtremal_; }
            ).second != nullptr;
}

bool SearchSpace::isKnownNonDependency(Vertical const& vertical, VerticalMap<VerticalInfo>* verticalInfos) {
    return verticalInfos->getAnySupersetEntry(
            vertical, [](auto vertical, auto info) -> bool { return !info->isDependency_; }).second != nullptr;
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
    strategy_->ensureInitialized(shared_from_this());
    std::string initializedLaunchPads;
    for (auto const& pad : launchPads_) {
        initializedLaunchPads += std::string(pad) + " ";
    }
    LOG(TRACE) << "Initialized with launch pads: " + initializedLaunchPads;
}








