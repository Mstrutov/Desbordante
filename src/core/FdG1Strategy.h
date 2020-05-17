#pragma once

#include <boost/format.hpp>
#include "DependencyStrategy.h"

class FdG1Strategy : DependencyStrategy {
private:
    std::shared_ptr<Column> rhs_;

    double calculateG1(std::shared_ptr<PositionListIndex> lhsPLI);
    double calculateG1(double numViolatingTuplePairs);
    ConfidenceInterval calculateG1(ConfidenceInterval const& numViolations);
public:
    FdG1Strategy(std::shared_ptr<Column> rhs, double maxError, double deviation) : DependencyStrategy(maxError, deviation), rhs_(rhs) {}

    void ensureInitialized(std::shared_ptr<SearchSpace> searchSpace) override;
    double calculateError(std::shared_ptr<Vertical> lhs) override;
    DependencyCandidate createDependencyCandidate(std::shared_ptr<Vertical> vertical) override;
    std::string format(std::shared_ptr<Vertical> vertical) override { return (boost::format("%s\u2192%s") % std::string(*vertical) % std::string(*rhs_)).str(); }
    explicit operator std::string() const override
        { return (boost::format("FD[RHS=%s, g1\u2264(%.3f..%.3f)]") % rhs_->getName() % minNonDependencyError_ % maxDependencyError_).str(); }
    void registerDependency(std::shared_ptr<Vertical> vertical, double error, const DependencyConsumer &discoveryUnit) override;
    bool isIrrelevantColumn(int columnIndex) override { return rhs_->getIndex() == columnIndex; }
    int getNumIrrelevantColumns() override { return 1; }
    Vertical getIrrelevantColumns() override { return static_cast<Vertical>(*rhs_); }
};
