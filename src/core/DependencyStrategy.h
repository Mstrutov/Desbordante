#pragma once
#include "ProfilingContext.h"
#include "DependencyCandidate.h"
#include "DependencyConsumer.h"
#include "Vertical.h"

class SearchSpace;

class DependencyStrategy {
protected:
    DependencyStrategy(double maxError, double deviation)
        : minNonDependencyError_(maxError - deviation),
          maxDependencyError_(maxError + deviation) {}
public:
    // TODO: public --> protected fields
    double minNonDependencyError_;
    double maxDependencyError_;
    ProfilingContext* context_;
    mutable unsigned int calcCount_ = 0;
    /*
     * Create the initial candidate for the given SearchSpace
     * */
    virtual void ensureInitialized(SearchSpace* searchSpace) const = 0;
    virtual DependencyCandidate createDependencyCandidate(Vertical const& candidate) const = 0;
    virtual double calculateError(Vertical const& candidate) const = 0;
    virtual std::string format(Vertical const& vertical) const = 0;
    explicit virtual operator std::string() const = 0;
    virtual void registerDependency(Vertical const& vertical, double error,
                                    DependencyConsumer const& discoveryUnit) const = 0;
    virtual bool isIrrelevantColumn(unsigned int columnIndex) const = 0;
    virtual unsigned int getNumIrrelevantColumns() const = 0;
    virtual Vertical getIrrelevantColumns() const = 0;

    virtual ~DependencyStrategy() = default;

    virtual std::unique_ptr<DependencyStrategy> createClone() = 0;

    bool shouldResample(Vertical const& vertical, double boostFactor) const;
    bool isIrrelevantColumn(Column const& column) const {
        return this->isIrrelevantColumn(column.getIndex());
    }

    static double round(double error) { return std::ceil(error * 32768) / 32768; }
};
