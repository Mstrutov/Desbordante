#pragma once
class VerticalInfo {
private:

public:
    double error_;
    bool isDependency_;
    bool isExtremal_;
    VerticalInfo(bool isDependency, bool isExtremal, double error = 0) :
        error_(error), isDependency_(isDependency), isExtremal_(isExtremal) {}

    static VerticalInfo forMinimalDependency() { return VerticalInfo(true, true); }
    static VerticalInfo forNonDependency() { return VerticalInfo(false, false); }
    static VerticalInfo forMaximalNonDependency() { return VerticalInfo(false, true); }

    bool isPruningSupersets() { return isDependency_ || isExtremal_; }
    bool isPruningSubsets() {return !isDependency_ || isExtremal_; }
};
