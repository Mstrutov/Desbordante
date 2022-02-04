#pragma once

#include <boost/format.hpp>
#include "DependencyStrategy.h"

class FdG1Strategy : public DependencyStrategy {
private:
    Column const* rhs_;

    double CalculateG1(util::PositionListIndex* lhs_pli) const;
    double CalculateG1(double num_violating_tuple_pairs) const;
    util::ConfidenceInterval CalculateG1(util::ConfidenceInterval const& num_violations) const;
public:
    static unsigned long long nanos_;

    FdG1Strategy(Column const* rhs, double max_error, double deviation)
        : DependencyStrategy(max_error, deviation), rhs_(rhs) {}

    void EnsureInitialized(SearchSpace* search_space) const override;
    double CalculateError(Vertical const& lhs) const override;
    DependencyCandidate CreateDependencyCandidate(Vertical const& vertical) const override;
    std::string Format(Vertical const& vertical) const override {
        return (boost::format("%s\u2192%s") % std::string(vertical) % std::string(*rhs_)).str();
    }
    explicit operator std::string() const override {
        return (boost::format("FD[RHS=%s, g1\u2264(%.3f..%.3f)]") %
            rhs_->GetName() % min_non_dependency_error_ % max_dependency_error_).str();
    }

    // TODO: can it be const though? Dependency registers --> some state somewhere changes.
    // Non-const discovery_unit?
    void RegisterDependency(Vertical const& vertical, double error,
                            DependencyConsumer const& discovery_unit) const override;
    bool IsIrrelevantColumn(unsigned int column_index) const override {
        return rhs_->GetIndex() == column_index;
    }
    unsigned int GetNumIrrelevantColumns() const override { return 1; }

    std::unique_ptr<DependencyStrategy> CreateClone() override;

    Vertical GetIrrelevantColumns() const override { return static_cast<Vertical>(*rhs_); }
};
