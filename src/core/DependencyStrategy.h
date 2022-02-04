#pragma once
#include "ProfilingContext.h"
#include "DependencyCandidate.h"
#include "DependencyConsumer.h"
#include "Vertical.h"

class SearchSpace;

class DependencyStrategy {
protected:
    DependencyStrategy(double max_error, double deviation)
        : min_non_dependency_error_(max_error - deviation),
          max_dependency_error_(max_error + deviation) {}
public:
    // TODO: public --> protected fields
    double min_non_dependency_error_;
    double max_dependency_error_;
    ProfilingContext* context_;
    mutable unsigned int calc_count_ = 0;
    /*
     * Create the initial candidate for the given SearchSpace
     * */
    virtual void EnsureInitialized(SearchSpace* search_space) const = 0;
    virtual DependencyCandidate CreateDependencyCandidate(Vertical const& candidate) const = 0;
    virtual double CalculateError(Vertical const& candidate) const = 0;
    virtual std::string Format(Vertical const& vertical) const = 0;
    explicit virtual operator std::string() const = 0;
    virtual void RegisterDependency(Vertical const& vertical, double error,
                                    DependencyConsumer const& discovery_unit) const = 0;
    virtual bool IsIrrelevantColumn(unsigned int column_index) const = 0;
    virtual unsigned int GetNumIrrelevantColumns() const = 0;
    virtual Vertical GetIrrelevantColumns() const = 0;

    virtual ~DependencyStrategy() = default;

    virtual std::unique_ptr<DependencyStrategy> CreateClone() = 0;

    bool ShouldResample(Vertical const& vertical, double boost_factor) const;
    bool IsIrrelevantColumn(Column const& column) const {
        return this->IsIrrelevantColumn(column.GetIndex());
    }

    static double Round(double error) { return std::ceil(error * 32768) / 32768; }
};
