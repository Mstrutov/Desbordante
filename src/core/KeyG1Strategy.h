#pragma once
#include <boost/format.hpp>
#include "DependencyStrategy.h"

class KeyG1Strategy : public DependencyStrategy {
private:
    double calculateKeyError(util::PositionListIndex* pli) const;
    double calculateKeyError(double numViolatingTuplePairs) const;
    util::ConfidenceInterval calculateKeyError(util::ConfidenceInterval const& numViolations) const;
public:
    KeyG1Strategy(double maxError, double deviation) : DependencyStrategy(maxError, deviation) {}

    void ensureInitialized(SearchSpace* searchSpace) const override;
    double calculateError(Vertical const& keyCandidate) const override;
    DependencyCandidate createDependencyCandidate(Vertical const& vertical) const override;
    std::string format(Vertical const& vertical) const override {
        return (boost::format("key(%s)") % std::string(vertical)).str(); }
    explicit operator std::string() const override {
        return (boost::format("key[g1\u2264(%.3f..%.3f)]") %
                minNonDependencyError_ % maxDependencyError_).str();
    }
    void registerDependency(Vertical const& vertical, double error,
                            DependencyConsumer const& discoveryUnit) const override;
    bool isIrrelevantColumn([[maybe_unused]] unsigned int columnIndex) const override {
        return false;
    }
    unsigned int getNumIrrelevantColumns() const override { return 1; }
    Vertical getIrrelevantColumns() const override {
        return *context_->getColumnLayoutRelationData()->getSchema()->emptyVertical;
    }

    std::unique_ptr<DependencyStrategy> createClone() override;
};
