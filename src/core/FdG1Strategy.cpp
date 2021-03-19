#include <unordered_map>

#include "FdG1Strategy.h"
#include "SearchSpace.h"
#include "PLICache.h"

#include "logging/easylogging++.h"

unsigned long long FdG1Strategy::nanos_ = 0;

double FdG1Strategy::calculateG1(PositionListIndex* lhsPLI) const {
    unsigned long long numViolations = 0;
    std::unordered_map<int, int> valueCounts;
    std::vector<int> const& probingTable =
            context_->getColumnLayoutRelationData()->getColumnData(rhs_->getIndex()).getProbingTable();

    LOG(DEBUG) << boost::format{"Probing table size for %1%: %2%"} % rhs_->toString()
        % std::to_string(probingTable.size());

    // Perform probing
    int probingTableValueId;
    for (auto const& cluster : lhsPLI->getIndex()) {
        valueCounts.clear();
        for (auto const& position : cluster) {
            probingTableValueId = probingTable[position];
            //    auto now = std::chrono::system_clock::now();
            if (probingTableValueId != PositionListIndex::singletonValueId) {
                auto location = valueCounts.find(probingTableValueId);
                if (location == valueCounts.end()) {
                    valueCounts.emplace_hint(location, probingTableValueId, 1);
                } else {
                    location->second += 1;
                }
            }
            //    nanos_ += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - now).count();
        }

        unsigned long long numViolationsInCluster = cluster.size() * (cluster.size() - 1) / 2;
        for (auto const& [key, refinedClusterSize] : valueCounts) {
            numViolationsInCluster -= static_cast<unsigned long long>(refinedClusterSize) * (refinedClusterSize - 1) / 2;
        }
        numViolations += numViolationsInCluster;
    }

    return calculateG1(numViolations);
}

double FdG1Strategy::calculateG1(double numViolatingTuplePairs) const {
    unsigned long long numTuplePairs = context_->getColumnLayoutRelationData()->getNumTuplePairs();
    if (numTuplePairs == 0) return 0;
    double g1 = numViolatingTuplePairs / numTuplePairs;
    return round(g1);
}

void FdG1Strategy::ensureInitialized(SearchSpace* searchSpace) const {
    if (searchSpace->isInitialized_) return;

    double zeroFdError = calculateError(*context_->getColumnLayoutRelationData()->getSchema()->emptyVertical);
    searchSpace->addLaunchPad(
            DependencyCandidate(
                    *context_->getColumnLayoutRelationData()->getSchema()->emptyVertical,
                    ConfidenceInterval(zeroFdError),
                    true
                    )
            );

    searchSpace->isInitialized_ = true;
}

double FdG1Strategy::calculateError(Vertical const& lhs) const {
    double error = 0;
    if (lhs.getArity() == 0) {
        auto rhsPli = context_->getPLICache()->get(static_cast<Vertical>(*rhs_));
        if (rhsPli == nullptr) throw std::runtime_error("Couldn't get rhsPLI from PLICache while calculating FD error");
        error = calculateG1(rhsPli->getNip());
    } else {
        auto lhsPli = context_->getPLICache()->getOrCreateFor(lhs, context_);
        auto jointPli = context_->getPLICache()->get(lhs.Union(static_cast<Vertical>(*rhs_)));
        error = jointPli == nullptr
                ? calculateG1(lhsPli)
                : calculateG1(lhsPli->getNepAsLong() - jointPli->getNepAsLong());
    }
    calcCount_++;
    return error;
}

ConfidenceInterval FdG1Strategy::calculateG1(ConfidenceInterval const &numViolations) const {
    return ConfidenceInterval(
            calculateG1(numViolations.getMin()),
            calculateG1(numViolations.getMean()),
            calculateG1(numViolations.getMax()));
}

DependencyCandidate FdG1Strategy::createDependencyCandidate(Vertical const& vertical) const {
    if (context_->isAgreeSetSamplesEmpty()) {
        return DependencyCandidate(vertical, ConfidenceInterval(0, .5, 1), false);
    }

    AgreeSetSample const* agreeSetSample = context_->getAgreeSetSample(vertical);
    ConfidenceInterval numViolatingTuplePairs = agreeSetSample
            ->estimateMixed(vertical, static_cast<Vertical>(*rhs_), context_->getConfiguration().estimateConfidence)
            .multiply(context_->getColumnLayoutRelationData()->getNumTuplePairs());
    //LOG(DEBUG) << boost::format{"Creating dependency candidate %1% with %2% violating pairs"}
    //    % vertical->toString() % numViolatingTuplePairs;
    ConfidenceInterval g1 = calculateG1(numViolatingTuplePairs);
    //LOG(DEBUG) << boost::format {"Creating dependency candidate %1% with error ~ %2%"}
    //    % vertical->toString() % g1;
    return DependencyCandidate(vertical, g1, false);
}

void FdG1Strategy::registerDependency(
        Vertical const& vertical, double error, DependencyConsumer const& discoveryUnit) const {
    discoveryUnit.registerFd(vertical, *rhs_, error, 0); // TODO: calculate score?
}

std::unique_ptr<DependencyStrategy> FdG1Strategy::createClone() {
    return std::make_unique<FdG1Strategy>(
        rhs_, (maxDependencyError_ + minNonDependencyError_) / 2,
        (maxDependencyError_ - minNonDependencyError_) / 2);
}