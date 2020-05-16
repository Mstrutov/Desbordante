#pragma once
class ConfidenceInterval {
private:
    double min_;
    double mean_;
    double max_;
public:
    //TODO: !(min <= mean <= max) -> exception
    ConfidenceInterval(double min, double mean, double max) :
        min_(min), mean_(mean), max_(max) {}

    ConfidenceInterval(double value) : ConfidenceInterval(value, value, value) {}

    double getMin() const { return min_; }
    double getMax() const { return max_; }
    double getMean() const { return mean_; }

    //TODO: assert min_ == max_
    double get() const { return mean_; }

    ConfidenceInterval multiply(double scalar) { return ConfidenceInterval(min_ * scalar, mean_ * scalar, max_ * scalar); }
    bool isPoint() { return min_ == max_; }
};
