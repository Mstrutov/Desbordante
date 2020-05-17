#pragma once
#include "core/ProfilingContext.h"
#include "core/DependencyCandidate.h"
#include "core/DependencyConsumer.h"
#include "model/Vertical.h"
class SearchSpace;

class DependencyStrategy {
protected:
    double maxDependencyError_;
    double minNonDependencyError_;
    std::shared_ptr<ProfilingContext> context_;

    DependencyStrategy(double maxError, double deviation) : maxDependencyError_(maxError + deviation), minNonDependencyError_(maxError - deviation) {}
public:
    /*
     * Create the initial candidate for the given SearchSpace
     * */
    virtual void ensureInitialized(std::shared_ptr<SearchSpace> searchSpace) = 0;
    virtual DependencyCandidate createDependencyCandidate(std::shared_ptr<Vertical> candidate) = 0;
    virtual double calculateError(std::shared_ptr<Vertical> candidate) = 0;
    virtual std::string format(std::shared_ptr<Vertical> vertical) = 0;
    explicit virtual operator std::string() const = 0;
    virtual void registerDependency(std::shared_ptr<Vertical> vertical, double error, DependencyConsumer const& discoveryUnit) = 0;
    virtual bool isIrrelevantColumn(int columnIndex) = 0;
    virtual int getNumIrrelevantColumns() = 0;
    virtual Vertical getIrrelevantColumns() = 0;

    virtual ~DependencyStrategy() = default;

    bool shouldResample(std::shared_ptr<Vertical> vertical, double boostFactor);
    bool isIrrelevantColumn(std::shared_ptr<Column> column) { return this->isIrrelevantColumn(column->getIndex()); }
};