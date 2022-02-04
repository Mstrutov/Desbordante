#pragma once
#include <boost/format.hpp>
#include <string>

namespace util {

class ConfidenceInterval {
private:
    double min_;
    double mean_;
    double max_;
public:
    //TODO: !(min <= mean <= max) -> exception
    ConfidenceInterval(double min, double mean, double max) :
        min_(min), mean_(mean), max_(max) {}

    explicit ConfidenceInterval(double value) : ConfidenceInterval(value, value, value) {}

    double GetMin() const { return min_; }
    double GetMax() const { return max_; }
    double GetMean() const { return mean_; }

    //TODO: assert min_ == max_
    double Get() const { return mean_; }

    ConfidenceInterval Multiply(double scalar) { return ConfidenceInterval(min_ * scalar, mean_ * scalar, max_ * scalar); }
    bool IsPoint() const { return min_ == max_; }
    explicit operator std::string() const { return (boost::format("error=(%f, %f, %f)") % min_ % mean_ % max_).str(); }

    friend std::ostream& operator<<(std::ostream&, ConfidenceInterval const&);
};

} // namespace util

