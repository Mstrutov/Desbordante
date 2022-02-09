#pragma once
#include <boost/format.hpp>
#include "DependencyStrategy.h"

class KeyG1Strategy : public DependencyStrategy {
private:
    double CalculateKeyError(util::PositionListIndex* pli) const;
    double CalculateKeyError(double num_violating_tuple_pairs) const;
    util::ConfidenceInterval CalculateKeyError(util::ConfidenceInterval const& num_violations) const;
public:
    KeyG1Strategy(double max_error, double deviation) : DependencyStrategy(max_error, deviation) {}

    void EnsureInitialized(SearchSpace* search_space) const override;
    double CalculateError(Vertical const& key_candidate) const override;
    DependencyCandidate CreateDependencyCandidate(Vertical const& vertical) const override;
    std::string Format(Vertical const& vertical) const override {
        return (boost::format("key(%s)") % std::string(vertical)).str();
    }
    explicit operator std::string() const override {
        return (boost::format("key[g1\u2264(%.3f..%.3f)]") % min_non_dependency_error_ %
                max_dependency_error_)
            .str();
    }
    void RegisterDependency(Vertical const& vertical, double error,
                            DependencyConsumer const& discovery_unit) const override;
    bool IsIrrelevantColumn([[maybe_unused]] unsigned int column_index) const override {
        return false;
    }
    unsigned int GetNumIrrelevantColumns() const override { return 1; }
    Vertical GetIrrelevantColumns() const override {
        return *context_->GetColumnLayoutRelationData()->GetSchema()->empty_vertical_;
    }

    std::unique_ptr<DependencyStrategy> CreateClone() override;
};
