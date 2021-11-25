#pragma once

#include <boost/format.hpp>
#include "DependencyStrategy.h"

class FdG1Strategy : public DependencyStrategy {
private:
    Column const* rhs_;

    double calculateG1(util::PositionListIndex* lhsPLI) const;
    double calculateG1(double numViolatingTuplePairs) const;
    util::ConfidenceInterval calculateG1(util::ConfidenceInterval const& numViolations) const;
public:
    static unsigned long long nanos_;

    FdG1Strategy(Column const* rhs, double maxError, double deviation)
        : DependencyStrategy(maxError, deviation), rhs_(rhs) {}

    void ensureInitialized(SearchSpace* searchSpace) const override;
    double calculateError(Vertical const& lhs) const override;
    DependencyCandidate createDependencyCandidate(Vertical const& vertical) const override;
    std::string format(Vertical const& vertical) const override {
        return (boost::format("%s\u2192%s") % std::string(vertical) % std::string(*rhs_)).str();
    }
    explicit operator std::string() const override {
        return (boost::format("FD[RHS=%s, g1\u2264(%.3f..%.3f)]") %
                rhs_->getName() % minNonDependencyError_ % maxDependencyError_).str();
    }

    // TODO: can it be const though? Dependency registers --> some state somewhere changes.
    // Non-const discoveryUnit?
    void registerDependency(Vertical const& vertical, double error,
                            DependencyConsumer const& discoveryUnit) const override;
    bool isIrrelevantColumn(unsigned int columnIndex) const override {
        return rhs_->getIndex() == columnIndex;
    }
    unsigned int getNumIrrelevantColumns() const override { return 1; }

    std::unique_ptr<DependencyStrategy> createClone() override;

    Vertical getIrrelevantColumns() const override { return static_cast<Vertical>(*rhs_); }
};
