#pragma once
class VerticalInfo {
private:

public:
    double error_;
    bool is_dependency_;
    mutable bool is_extremal_;
    VerticalInfo(bool is_dependency, bool is_extremal, double error = 0) :
        error_(error), is_dependency_(is_dependency), is_extremal_(is_extremal) {}

    VerticalInfo(VerticalInfo const& other) = default;

    static VerticalInfo ForMinimalDependency() { return VerticalInfo(true, true); }
    static VerticalInfo ForNonDependency() { return VerticalInfo(false, false); }
    static VerticalInfo ForMaximalNonDependency() { return VerticalInfo(false, true); }

    bool IsPruningSupersets() const { return is_dependency_ || is_extremal_; }
    bool IsPruningSubsets() const {return !is_dependency_ || is_extremal_; }
};
