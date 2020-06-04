#include "AgreeSetSample.h"

#include <algorithm>
#include <exception>
#include <limits>
#include <memory>
#include <random>
#include <utility>

using namespace std;

double AgreeSetSample::stdDevSmoothing = 1;

AgreeSetSample::AgreeSetSample(shared_ptr<ColumnLayoutRelationData> relationData, shared_ptr<Vertical> focus, int sampleSize,
                               long populationSize):
        relationData(std::move(relationData)),
        focus(std::move(focus)),
        sampleSize(sampleSize),
        populationSize(populationSize){
}



double AgreeSetSample::calculateNonNegativeFraction(double a, double b) {
    if (a == 0)
        return 0;
    return max(std::numeric_limits<double>::min(), a / b);
}

//TODO: is this manipulation with virtual and pure virtual methods even legal?
std::shared_ptr<std::vector<long long>>
AgreeSetSample::getNumAgreeSupersetsExt(std::shared_ptr<Vertical> agreement, std::shared_ptr<Vertical> disagreement) {
    return std::make_shared<std::vector<long long>> (std::vector<long long> {this->getNumAgreeSupersets(agreement), this->getNumAgreeSupersets(agreement, disagreement)});
}

double AgreeSetSample::estimateAgreements(std::shared_ptr<Vertical> agreement) {
    if (!agreement->contains(*this->focus)) {
        throw std::runtime_error("An agreement in estimateAgreemnts should contain the focus");
    }

    if (populationSize == 0) {
        return 0;
    }
    return observationsToRelationRatio(this->getNumAgreeSupersets(agreement));
}

ConfidenceInterval AgreeSetSample::estimateAgreements(std::shared_ptr<Vertical> agreement, double confidence) {
    if (!agreement->contains(*this->focus)) {
        throw std::runtime_error("An agreement in estimateAgreemnts with confidence should contain the focus");
    }
    if (populationSize == 0) {
        return ConfidenceInterval(0, 0, 0);
    }

    //Counting the sampled tuples agreeing as requested - calling virtual method
    long long numHits = this->getNumAgreeSupersets(agreement);

    return estimateGivenNumHits(numHits, confidence);
}

ConfidenceInterval AgreeSetSample::estimateMixed(std::shared_ptr<Vertical> agreement, std::shared_ptr<Vertical> disagreement, double confidence) {
    if (!agreement->contains(*this->focus)) {
        throw std::runtime_error("An agreement in estimateMixed should contain the focus");
    }
    if (populationSize == 0) {
        return ConfidenceInterval(0, 0, 0);
    }

    //Counting the sampled tuples agreeing as requested - calling virtual method
    long long numHits = this->getNumAgreeSupersets(agreement, disagreement);

    return estimateGivenNumHits(numHits, confidence);
}

ConfidenceInterval AgreeSetSample::estimateGivenNumHits(long long numHits, double confidence) {
    double sampleRatio = numHits / static_cast<double>(sampleSize);
    double relationRatio = ratioToRelationRatio(sampleRatio);
    if (this->isExact() || confidence == -1) {       //TODO: check all the stuff with confidence interval and what confidence is; possibility to use boost::optional
        return ConfidenceInterval(relationRatio);
    }

    normal_distribution normalDistribution;
    double z = probitFunction((confidence + 1) / 2);
    double smoothedSampleRatio = (numHits + stdDevSmoothing / 2) / (sampleSize + stdDevSmoothing);
    double stdDevPositiveTuples = sqrt(smoothedSampleRatio * (1 - smoothedSampleRatio) / sampleSize);
    double minRatio = max(sampleRatio - z * stdDevPositiveTuples, calculateNonNegativeFraction(numHits, relationData->getNumTuplePairs()));
    double maxRatio = sampleRatio + z * stdDevPositiveTuples;

    return ConfidenceInterval(ratioToRelationRatio(minRatio), relationRatio, ratioToRelationRatio(maxRatio));
}

// Inverse cumulative distribution function of normal distribution
//taken from https://www.quantstart.com/articles/Statistical-Distributions-in-C/
double AgreeSetSample::probitFunction(double quantile) const {
    // This is the Beasley-Springer-Moro algorithm which can
    // be found in Glasserman [2004].
    static double a[4] = {   2.50662823884,
                             -18.61500062529,
                             41.39119773534,
                             -25.44106049637};

    static double b[4] = {  -8.47351093090,
                            23.08336743743,
                            -21.06224101826,
                            3.13082909833};

    static double c[9] = {0.3374754822726147,
                          0.9761690190917186,
                          0.1607979714918209,
                          0.0276438810333863,
                          0.0038405729373609,
                          0.0003951896511919,
                          0.0000321767881768,
                          0.0000002888167364,
                          0.0000003960315187};

    if (quantile >= 0.5 && quantile <= 0.92) {
        double num = 0.0;
        double denom = 1.0;

        for (int i=0; i<4; i++) {
            num += a[i] * pow((quantile - 0.5), 2*i + 1);
            denom += b[i] * pow((quantile - 0.5), 2*i);
        }
        return num/denom;

    } else if (quantile > 0.92 && quantile < 1) {
        double num = 0.0;

        for (int i=0; i<9; i++) {
            num += c[i] * pow((log(-log(1-quantile))), i);
        }
        return num;

    } else {
        return -1.0*probitFunction(1-quantile);
    }
}
