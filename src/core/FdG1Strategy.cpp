#include <unordered_map>
#include "FdG1Strategy.h"
#include "SearchSpace.h"
#include "PLICache.h"

#include "logging/easylogging++.h"

double FdG1Strategy::calculateG1(std::shared_ptr<PositionListIndex> lhsPLI) {
    unsigned long long numViolations = 0;
    std::map<int, int> valueCounts;
    std::vector<int> probingTable = context_->relationData_->getColumnData(rhs_->getIndex())->getProbingTable();

    // Perform probing
    for (auto const& cluster : lhsPLI->getIndex()) {
        valueCounts.clear();
        for (auto const& position : cluster) {
            int probingTableValueId = probingTable[position];

            if (probingTableValueId != PositionListIndex::singletonValueId) {
                if (valueCounts.find(probingTableValueId) == valueCounts.end()) {
                    valueCounts.emplace(probingTableValueId, 0);
                }
                valueCounts[probingTableValueId] += 1;
            }
        }

        unsigned long long numViolationsInCluster = cluster.size() * (cluster.size() - 1) / 2;
        for (auto const& [key, refinedClusterSize] : valueCounts) {
            numViolationsInCluster -= static_cast<unsigned long long>(refinedClusterSize) * (refinedClusterSize - 1) / 2;
        }
        numViolations += numViolationsInCluster;
    }

    return calculateG1(numViolations);
}

double FdG1Strategy::calculateG1(double numViolatingTuplePairs) {
    unsigned long long numTuplePairs = context_->relationData_->getNumTuplePairs();
    if (numTuplePairs == 0) return 0;
    double g1 = numViolatingTuplePairs / numTuplePairs;
    return round(g1);
}

void FdG1Strategy::ensureInitialized(std::shared_ptr<SearchSpace> searchSpace) {
    if (searchSpace->isInitialized_) return;

    double zeroFdError = calculateError(context_->relationData_->getSchema()->emptyVertical);
    searchSpace->addLaunchPad(
            DependencyCandidate(
                    context_->relationData_->getSchema()->emptyVertical,
                    ConfidenceInterval(zeroFdError),
                    true
                    )
            );

    searchSpace->isInitialized_ = true;
}

double FdG1Strategy::calculateError(std::shared_ptr<Vertical> lhs) {
    double error = 0;
    if (lhs->getArity() == 0) {
        auto rhsPli = context_->pliCache_->get(static_cast<Vertical>(*rhs_));
        if (rhsPli == nullptr) throw std::runtime_error("Couldn't get rhsPLI from PLICache while calculating FD error");
        error = calculateG1(rhsPli->getNip());
    } else {
        auto lhsPli = context_->pliCache_->getOrCreateFor(*lhs, *context_);
        auto jointPli = context_->pliCache_->get(*lhs->Union(static_cast<Vertical>(*rhs_)));
        error = jointPli == nullptr ? calculateG1(lhsPli) : calculateG1(lhsPli->getNepAsLong() - jointPli->getNepAsLong());
    }
    calcCount_++;
    return error;
}

ConfidenceInterval FdG1Strategy::calculateG1(ConfidenceInterval const &numViolations) {
    return ConfidenceInterval(calculateG1(numViolations.getMin()), calculateG1(numViolations.getMean()), calculateG1(numViolations.getMax()));
}

DependencyCandidate FdG1Strategy::createDependencyCandidate(std::shared_ptr<Vertical> vertical) {
    if (context_->agreeSetSamples_ == nullptr) {
        return DependencyCandidate(vertical, ConfidenceInterval(0, .5, 1), false);
    }

    std::shared_ptr<AgreeSetSample> agreeSetSample = context_->getAgreeSetSample(vertical);
    ConfidenceInterval numViolatingTuplePairs = agreeSetSample
            ->estimateMixed(vertical, std::make_shared<Vertical>(static_cast<Vertical>(*rhs_)),context_->configuration_.estimateConfidence)
            .multiply(context_->relationData_->getNumTuplePairs());
    //LOG(DEBUG) << boost::format{"Creating dependency candidate %1% with %2% violating pairs"}
    //    % vertical->toString() % numViolatingTuplePairs;
    ConfidenceInterval g1 = calculateG1(numViolatingTuplePairs);
    //LOG(DEBUG) << boost::format {"Creating dependency candidate %1% with error ~ %2%"}
    //    % vertical->toString() % g1;
    return DependencyCandidate(vertical, g1, false);
}

void FdG1Strategy::registerDependency(std::shared_ptr<Vertical> vertical, double error,
                                      DependencyConsumer const& discoveryUnit) {
    discoveryUnit.registerFd(vertical, rhs_, error, 0); // TODO: calculate score?
}
