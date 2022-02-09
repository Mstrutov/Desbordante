//
// Created by Ilya Vologin
// https://github.com/cupertank
//

#pragma once

#include <boost/dynamic_bitset.hpp>

#include "ColumnLayoutRelationData.h"
#include "ConfidenceInterval.h"
#include "custom/CustomRandom.h"
#include "Vertical.h"

namespace util {

//abstract base class for Agree Set Sample implementations (trie <- not used, list)
class AgreeSetSample {
public:

    virtual unsigned long long GetNumAgreeSupersets(Vertical const& agreement) const = 0;
    virtual unsigned long long GetNumAgreeSupersets(Vertical const& agreement,
                                                    Vertical const& disagreement) const = 0;
    virtual std::unique_ptr<std::vector<unsigned long long>> GetNumAgreeSupersetsExt(
        Vertical const& agreement, Vertical const& disagreement) const;

    double EstimateAgreements(Vertical const& agreement) const;
    ConfidenceInterval EstimateAgreements(Vertical const& agreement, double confidence) const;
    ConfidenceInterval EstimateMixed(Vertical const& agreement, Vertical const& disagreement,
                                     double confidence) const;

    double GetSamplingRatio() const { return sample_size_ / static_cast<double>(population_size_); }
    bool IsExact() const { return population_size_ == sample_size_; }

    virtual ~AgreeSetSample() = default;

protected:
    ::ColumnLayoutRelationData const* relation_data_;
    Vertical focus_;
    unsigned int sample_size_;
    unsigned long long population_size_;
    AgreeSetSample(ColumnLayoutRelationData const* relation_data, Vertical focus,
                   unsigned int sample_size, unsigned long long population_size);

    template<typename T>
    static std::unique_ptr<T> CreateFor(ColumnLayoutRelationData* relation_data, int sample_size);

    template<typename T>
    static std::unique_ptr<T> CreateFocusedFor(ColumnLayoutRelationData const* relation,
                                               Vertical const& restriction_vertical,
                                               PositionListIndex const* restriction_pli,
                                               unsigned int sample_size, CustomRandom& random);
private:
    static double std_dev_smoothing_;

    double RatioToRelationRatio(double ratio) const {
        return ratio * population_size_ / relation_data_->GetNumTuplePairs();
    }
    double ObservationsToRelationRatio(double num_observations) const {
        return RatioToRelationRatio(num_observations / sample_size_);
    }
    static double CalculateNonNegativeFraction(double a, double b);

    ConfidenceInterval EstimateGivenNumHits(unsigned long long num_hits, double confidence) const;
    // Inverse cumulative distribution function (aka the probit function)
    double ProbitFunction(double quantile) const;
};

} // namespace util

//include template implementation
#include "AgreeSetSample_impl.h"

