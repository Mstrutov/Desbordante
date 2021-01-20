#pragma once
#include "ProfilingContext.h"
#include "DependencyCandidate.h"
#include "DependencyConsumer.h"
#include "Vertical.h"
class SearchSpace;

class DependencyStrategy {
protected:

    DependencyStrategy(double maxError, double deviation) :
        minNonDependencyError_(maxError - deviation), maxDependencyError_(maxError + deviation) {}
public:
    double minNonDependencyError_;
    double maxDependencyError_;
    std::shared_ptr<ProfilingContext> context_;
    unsigned int calcCount_ = 0;
    /*
     * Create the initial candidate for the given SearchSpace
     * */
    virtual void ensureInitialized(std::shared_ptr<SearchSpace> searchSpace) = 0;
    virtual DependencyCandidate createDependencyCandidate(std::shared_ptr<Vertical> candidate) = 0;
    virtual double calculateError(std::shared_ptr<Vertical> candidate) = 0;
    virtual std::string format(std::shared_ptr<Vertical> vertical) = 0;
    explicit virtual operator std::string() const = 0;
    virtual void registerDependency(std::shared_ptr<Vertical> vertical, double error, DependencyConsumer const& discoveryUnit) = 0;
    virtual bool isIrrelevantColumn(unsigned int columnIndex) = 0;
    virtual unsigned int getNumIrrelevantColumns() = 0;
    virtual Vertical getIrrelevantColumns() = 0;

    virtual ~DependencyStrategy() = default;

    bool shouldResample(std::shared_ptr<Vertical> vertical, double boostFactor);
    bool isIrrelevantColumn(std::shared_ptr<Column> column) { return this->isIrrelevantColumn(column->getIndex()); }

    static double round(double error) { return std::ceil(error * 32768) / 32768; }
};