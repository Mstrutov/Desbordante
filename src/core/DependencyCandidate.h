#pragma once
#include <functional>
#include "model/Vertical.h"
#include "util/ConfidenceInterval.h"

class DependencyCandidate {
private:
    std::shared_ptr<Vertical> vertical_;
    ConfidenceInterval error_;
    bool isExact_;

public:
    DependencyCandidate(std::shared_ptr<Vertical> vertical, ConfidenceInterval error, bool isExact) :
        vertical_(vertical), error_(error), isExact_(isExact) {}

    bool operator<(DependencyCandidate const& other) const;
    
    //TODO: implement if used
    //bool operator==(DependencyCandidate const& other) const;

    bool isExact() { return isExact_ && error_.isPoint(); }

    //TODO: check if these comparators are used in a right way
    //using comparator = std::function<bool (DependencyCandidate const &, DependencyCandidate const &)>;
    static bool arityComparator(DependencyCandidate const &, DependencyCandidate const &);
    static bool minErrorComparator(DependencyCandidate const &, DependencyCandidate const &);
    //static bool maxErrorComparator(DependencyCandidate const &, DependencyCandidate const &);
    static bool fullArityErrorComparator(DependencyCandidate const &, DependencyCandidate const &);
    static bool fullErrorArityComparator(DependencyCandidate const &, DependencyCandidate const &);
};
