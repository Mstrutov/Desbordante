#pragma once
class VerticalInfo {
private:
    bool isDependency_;
    bool isExtremal_;
    double error_;

public:
    VerticalInfo(bool isDependency, bool isExtremal, double error = 0) :
        isDependency_(isDependency), isExtremal_(isExtremal), error_(error) {}

    static VerticalInfo forMinimalDependency() { return VerticalInfo(true, true); }
    static VerticalInfo forNonDependency() { return VerticalInfo(false, false); }
    static VerticalInfo forMaximalNonDependency() { return VerticalInfo(false, true); }

    bool isPruningSupersets() { return isDependency_ || isExtremal_; }
    bool isPruningSubsets() {return !isDependency_ || isExtremal_; }
};
