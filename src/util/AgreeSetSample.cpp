#include "AgreeSetSample.h"

#include <algorithm>
#include <exception>
#include <limits>
#include <memory>
#include <random>
#include <utility>

namespace util {

using namespace std;

double AgreeSetSample::std_dev_smoothing_ = 1;

AgreeSetSample::AgreeSetSample(ColumnLayoutRelationData const* relation_data, Vertical focus, unsigned int sample_size,
                               unsigned long long population_size):
    relation_data_(relation_data),
    focus_(std::move(focus)),
    sample_size_(sample_size),
    population_size_(population_size){
}



double AgreeSetSample::CalculateNonNegativeFraction(double a, double b) {
    // TODO: checking 0/b, comparing double to 0.
    if (a == 0)
        return 0;
    return max(std::numeric_limits<double>::min(), a / b);
}

std::unique_ptr<std::vector<unsigned long long>>
AgreeSetSample::GetNumAgreeSupersetsExt(Vertical const& agreement, Vertical const& disagreement) const {
    return std::make_unique<std::vector<unsigned long long>> (
            std::vector<unsigned long long> {
                this->GetNumAgreeSupersets(agreement), this->GetNumAgreeSupersets(agreement, disagreement)
            }
        );
}

double AgreeSetSample::EstimateAgreements(Vertical const& agreement) const {
    if (!agreement.Contains(this->focus_)) {
        throw std::runtime_error("An agreement in estimateAgreemnts should contain the focus");
    }

    if (population_size_ == 0) {
        return 0;
    }
    return ObservationsToRelationRatio(this->GetNumAgreeSupersets(agreement));
}

ConfidenceInterval AgreeSetSample::EstimateAgreements(Vertical const& agreement, double confidence) const {
    if (!agreement.Contains(this->focus_)) {
        throw std::runtime_error("An agreement in estimateAgreemnts with confidence should contain the focus");
    }
    if (population_size_ == 0) {
        return ConfidenceInterval(0, 0, 0);
    }

    //Counting the sampled tuples agreeing as requested - calling virtual method
    long long num_hits = this->GetNumAgreeSupersets(agreement);

    return EstimateGivenNumHits(num_hits, confidence);
}

ConfidenceInterval AgreeSetSample::EstimateMixed(
        Vertical const& agreement, Vertical const& disagreement, double confidence) const {
    if (!agreement.Contains(this->focus_)) {
        throw std::runtime_error("An agreement in EstimateMixed should contain the focus");
    }
    if (population_size_ == 0) {
        return ConfidenceInterval(0, 0, 0);
    }

    //Counting the sampled tuples agreeing as requested - calling virtual method
    long long num_hits = this->GetNumAgreeSupersets(agreement, disagreement);

    return EstimateGivenNumHits(num_hits, confidence);
}

ConfidenceInterval AgreeSetSample::EstimateGivenNumHits(unsigned long long num_hits, double confidence) const {
    double sample_ratio = num_hits / static_cast<double>(sample_size_);
    double relation_ratio = RatioToRelationRatio(sample_ratio);
    if (this->IsExact() || confidence == -1) {
        //TODO: check all the stuff with confidence interval and what confidence is; possibility to use boost::optional
        return ConfidenceInterval(relation_ratio);
    }

    normal_distribution normal_distribution;
    double z = ProbitFunction((confidence + 1) / 2);
    double smoothed_sample_ratio = (num_hits + std_dev_smoothing_ / 2) / (sample_size_ + std_dev_smoothing_);
    double std_dev_positive_tuples = sqrt(smoothed_sample_ratio * (1 - smoothed_sample_ratio) / sample_size_);
    double min_ratio = max(sample_ratio - z * std_dev_positive_tuples,
                           CalculateNonNegativeFraction(num_hits, relation_data_->GetNumTuplePairs()));
    double max_ratio = sample_ratio + z * std_dev_positive_tuples;

    return ConfidenceInterval(RatioToRelationRatio(min_ratio), relation_ratio, RatioToRelationRatio(max_ratio));
}

// Inverse cumulative distribution function of normal distribution
// taken from https://www.quantstart.com/articles/Statistical-Distributions-in-C/
double AgreeSetSample::ProbitFunction(double quantile) const {
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
        return -1.0* ProbitFunction(1 - quantile);
    }
}

} // namespace util

