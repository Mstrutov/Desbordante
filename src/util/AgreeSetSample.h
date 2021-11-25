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

    virtual unsigned long long getNumAgreeSupersets(Vertical const& agreement) const = 0;
    virtual unsigned long long getNumAgreeSupersets(Vertical const& agreement, Vertical const& disagreement) const = 0;
    virtual std::unique_ptr<std::vector<unsigned long long>> getNumAgreeSupersetsExt(
            Vertical const& agreement, Vertical const& disagreement) const;

    double estimateAgreements(Vertical const& agreement) const;
    ConfidenceInterval estimateAgreements(Vertical const& agreement, double confidence) const;
    ConfidenceInterval estimateMixed(Vertical const& agreement, Vertical const& disagreement, double confidence) const;

    double getSamplingRatio() const { return sampleSize / static_cast<double>(populationSize); }
    bool isExact() const { return populationSize == sampleSize; }

    virtual ~AgreeSetSample() = default;

protected:
    ::ColumnLayoutRelationData const* relationData;
    Vertical focus;
    unsigned int sampleSize;
    unsigned long long populationSize;
    AgreeSetSample(ColumnLayoutRelationData const* relationData, Vertical  focus, unsigned int sampleSize, unsigned long long populationSize);

    template<typename T>
    static std::unique_ptr<T> createFor(ColumnLayoutRelationData* relationData, int sampleSize);

    template<typename T>
    static std::unique_ptr<T> createFocusedFor(ColumnLayoutRelationData const* relation,
                                          Vertical const& restrictionVertical,
                                          PositionListIndex const* restrictionPli,
                                          unsigned int sampleSize, CustomRandom& random);
private:
    static double stdDevSmoothing;

    double ratioToRelationRatio(double ratio) const {
        return ratio * populationSize / relationData->getNumTuplePairs(); }
    double observationsToRelationRatio(double numObservations) const {
        return ratioToRelationRatio(numObservations / sampleSize); }
    static double calculateNonNegativeFraction(double a, double b);

    ConfidenceInterval estimateGivenNumHits(unsigned long long numHits, double confidence) const;
    // Inverse cumulative distribution function (aka the probit function)
    double probitFunction(double quantile) const;
};

} // namespace util

//include template implementation
#include "AgreeSetSample_impl.h"

