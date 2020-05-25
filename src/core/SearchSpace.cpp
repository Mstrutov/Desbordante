#include "SearchSpace.h"

/*SearchSpace::SearchSpace(int id, std::shared_ptr<DependencyStrategy> strategy,
                         std::unique_ptr<VerticalMap<Vertical>> scope, VerticalMap<VerticalInfo> globalVisitees,
                         std::shared_ptr<RelationalSchema> schema, std::function<bool(DependencyCandidate const &,
                                                                                      DependencyCandidate const &)> const &dependencyCandidateComparator,
                         int recursionDepth, double sampleBoost) : id_(id), strategy_(strategy), scope_(std::move(scope)), globalVisitees_(globalVisitees), recursionDepth_(recursionDepth),
                                                                   sampleBoost_(sampleBoost), launchPadIndex_(schema), launchPads_(dependencyCandidateComparator) {}*/

void SearchSpace::discover(std::shared_ptr<VerticalMap<std::shared_ptr<VerticalInfo>>> localVisitees) {
    while (true) {
        std::shared_ptr<DependencyCandidate> launchPad = pollLaunchPad(localVisitees);
        if (launchPad == nullptr) break;

        localVisitees = localVisitees != nullptr ? localVisitees : std::make_shared<VerticalMap<std::shared_ptr<VerticalInfo>>>(context_->getSchema());
        bool isDependencyFound = ascend(*launchPad, localVisitees);

        returnLaunchPad(*launchPad, !isDependencyFound);
    }
}

// этот метод больше для распараллеленного варианта
std::shared_ptr<DependencyCandidate> SearchSpace::pollLaunchPad(
        std::shared_ptr<VerticalMap<std::shared_ptr<VerticalInfo>>> localVisitees) {
    while (true) {
        // TODO: boost::optional
        std::shared_ptr<DependencyCandidate> launchPad;
        if (launchPads_.empty()) {
            if (deferredLaunchPads_.empty()) return nullptr;

            launchPads_.insert(deferredLaunchPads_.begin(), deferredLaunchPads_.end());
            deferredLaunchPads_.clear();
        }

        if (isImpliedByMinDep(launchPad->vertical_, globalVisitees_)
            || (localVisitees != nullptr && isImpliedByMinDep(launchPad->vertical_, localVisitees))) {
            launchPadIndex_.remove(*launchPad->vertical_);
            continue;
        }

        auto supersetEntries = globalVisitees_->getSupersetEntries(*launchPad->vertical_);
        if (localVisitees != nullptr) {
            auto localSupersetEntries = localVisitees->getSupersetEntries(*launchPad->vertical_);
            auto endIterator = std::remove_if(localSupersetEntries.begin(), localSupersetEntries.end(), [](auto entry) { return entry.second->isPruningSubsets(); });
            std::for_each(localSupersetEntries.begin(), localSupersetEntries.end(), [&supersetEntries](auto entry) { supersetEntries.push_back(entry); });
        }
        if (supersetEntries.empty()) return launchPad;

        std::list<std::shared_ptr<Vertical>> supersetVerticals;
        std::transform(supersetEntries.begin(), supersetEntries.end(), supersetVerticals,
                [](auto entry) {return std::make_shared<Vertical>(entry->second); });
        escapeLaunchPad(launchPad->vertical_, std::move(supersetVerticals), localVisitees);
        //continue;

        //launchPads_.insert(*launchPad);
    }
}

// this move looks legit IMO
void SearchSpace::escapeLaunchPad(std::shared_ptr<Vertical> launchPad,
                                  std::list<std::shared_ptr<Vertical>>&& pruningSupersets,
                                  std::shared_ptr<VerticalMap<std::shared_ptr<VerticalInfo>>> localVisitees) {
    // TODO: list<не shared_ptr>, чтобы не создавать новые поинтеры? Мб создать новый лист list<Vertical>
    std::transform(pruningSupersets.begin(), pruningSupersets.end(), pruningSupersets.begin(),
            [this](auto superset) { return std::make_shared<Vertical>(superset->invert().without(strategy_->getIrrelevantColumns())); } );

    std::function<bool (Vertical const&)> pruningFunction = [this, &launchPad, &localVisitees] (auto const& hittingSetCandidate) -> bool {
        if (scope_ != nullptr && scope_->getAnySupersetEntry(hittingSetCandidate)) {
            return true;
        }

        Vertical launchPadCandidate = launchPad->Union(hittingSetCandidate);

        // TODO: после создания такого поинтера может произойти непреднамеренное удаление объекта при удалении поинтера?
        auto launchPadCandidate_ptr = std::make_shared<Vertical>(launchPadCandidate);
        if ((localVisitees == nullptr && isImpliedByMinDep(launchPadCandidate_ptr, localVisitees))
            || isImpliedByMinDep(launchPadCandidate_ptr, globalVisitees_)) {
            return true;
        }

        if (launchPadIndex_.getAnySubsetEntry(launchPadCandidate).second != nullptr) {
            return true;
        }

        return false;
    };
    auto hittingSet = context_->getSchema()->calculateHittingSet(
            std::move(pruningSupersets),
            boost::make_optional(pruningFunction)
        );

    for (auto& escaping : hittingSet) {
        auto escapedLaunchPadVertical = std::make_shared<Vertical>(launchPad->Union(escaping));
        // assert, который не имплементнуть из-за трансформа
        DependencyCandidate escapedLaunchPad = strategy_->createDependencyCandidate(escapedLaunchPadVertical);
        launchPads_.insert(escapedLaunchPad);
        launchPadIndex_.put(*escapedLaunchPad.vertical_, std::make_shared<DependencyCandidate>(escapedLaunchPad));
    }
}

void SearchSpace::addLaunchPad(const DependencyCandidate &launchPad) {
    launchPads_.insert(launchPad);
    launchPadIndex_.put(*launchPad.vertical_, std::make_shared<DependencyCandidate>(launchPad));
}

void SearchSpace::returnLaunchPad(DependencyCandidate const &launchPad, bool isDefer) {
    if (isDefer && context_->configuration_.isDeferFailedLaunchPads) {
        deferredLaunchPads_.push_back(launchPad);
    }
    else {
        launchPads_.insert(launchPad);
    }
    launchPadIndex_.put(*launchPad.vertical_, std::make_shared<DependencyCandidate>(launchPad));
}

bool SearchSpace::ascend(DependencyCandidate const &launchPad,
                         std::shared_ptr<VerticalMap<std::shared_ptr<VerticalInfo>>> localVisitees) {
    if (strategy_->shouldResample(launchPad.vertical_, sampleBoost_)) {
        context_->createFocusedSample(launchPad.vertical_, sampleBoost_);
    }

    DependencyCandidate traversalCandidate = launchPad;
    boost::optional<double> error;
    while (true) {
        if (context_->configuration_.isCheckEstimates) {
            checkEstimate(strategy_, traversalCandidate);
        }

        if (traversalCandidate.isExact()) {
            // TODO: getError()?
            error = traversalCandidate.error_.get();

            bool canBeDependency = error <= strategy_->maxDependencyError_;
            localVisitees->put(*traversalCandidate.vertical_, std::make_shared<VerticalInfo>(
                    canBeDependency,
                    false,
                    *error
                    ));
            if (canBeDependency) break;
        }
        else {
            if (traversalCandidate.error_.getMin() > strategy_->maxDependencyError_) {
                error.reset();
            }
            else {
                error = context_->configuration_.isEstimateOnly
                        ? traversalCandidate.error_.getMean()
                        : strategy_->calculateError(traversalCandidate.vertical_);
                double errorDiff = *error - traversalCandidate.error_.getMean();
                // TODO:: context_->profilingData is only for logging, right?

                localVisitees->put(*traversalCandidate.vertical_, std::make_shared<VerticalInfo>(
                        error <= strategy_->maxDependencyError_,
                        false,
                        *error
                        ));

                if (error <= strategy_->maxDependencyError_) break;

                if (strategy_->shouldResample(traversalCandidate.vertical_, sampleBoost_)) {
                    context_->createFocusedSample(traversalCandidate.vertical_, sampleBoost_);
                }
            }
        }

        if (traversalCandidate.vertical_->getArity() >= context_->relationData_->getNumColumns() - strategy_->getNumIrrelevantColumns()) {
            break;
        }

        boost::optional<DependencyCandidate> nextCandidate;
        int numSeenElements = isAscendRandomly_ ? 1 : -1;
        for (auto& extensionColumn : context_->getSchema()->getColumns()) {
            if (traversalCandidate.vertical_->getColumnIndices()[extensionColumn->getIndex()]
                    || strategy_->isIrrelevantColumn(extensionColumn)) {
                continue;
            }
            auto extendedVertical = std::make_shared<Vertical>(
                    traversalCandidate.vertical_->Union(static_cast<Vertical>(*extensionColumn)));

            if (scope_ != nullptr && scope_->getSupersetEntries(*extendedVertical).empty()) {
                continue;
            }

            bool isSubsetPruned = isImpliedByMinDep(extendedVertical, globalVisitees_);
            if (isSubsetPruned) {
                continue;
            }
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

    if (!error) {
        error = strategy_->calculateError(traversalCandidate.vertical_);
        double errorDiff = *error - traversalCandidate.error_.getMean();
    }

    if (error <= strategy_->maxDependencyError_) {
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
        } else {
            localVisitees->put(*traversalCandidate.vertical_, std::make_shared<VerticalInfo>(VerticalInfo::forNonDependency()));
        }
    }
    return false;
}

