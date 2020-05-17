#pragma once
#include <boost/format.hpp>
#include "DependencyStrategy.h"

class KeyG1Strategy : DependencyStrategy {
private:
    double calculateKeyError(std::shared_ptr<PositionListIndex> pli);
    double calculateKeyError(double numViolatingTuplePairs);
    ConfidenceInterval calculateKeyError(ConfidenceInterval const& numViolations);
public:
    KeyG1Strategy(double maxError, double deviation) : DependencyStrategy(maxError, deviation) {}

    void ensureInitialized(std::shared_ptr<SearchSpace> searchSpace) override;
    double calculateError(std::shared_ptr<Vertical> keyCandidate) override;
    DependencyCandidate createDependencyCandidate(std::shared_ptr<Vertical> vertical) override;
    std::string format(std::shared_ptr<Vertical> vertical) override { return (boost::format("key(%s)") % std::string(*vertical)).str(); }
    explicit operator std::string() const override
    { return (boost::format("key[g1\u2264(%.3f..%.3f)]") % minNonDependencyError_ % maxDependencyError_).str(); }
    void registerDependency(std::shared_ptr<Vertical> vertical, double error, const DependencyConsumer &discoveryUnit) override;
    bool isIrrelevantColumn(int columnIndex) override { return false; }
    int getNumIrrelevantColumns() override { return 0; }
    Vertical getIrrelevantColumns() override { return *context_->relationData_->getSchema()->emptyVertical; }
};
