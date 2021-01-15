#pragma once
#include <functional>
#include "Vertical.h"
#include "ConfidenceInterval.h"

class DependencyCandidate {
private:
    bool isExact_;

public:
    ConfidenceInterval error_;
    std::shared_ptr<Vertical> vertical_;

    DependencyCandidate(std::shared_ptr<Vertical> vertical, ConfidenceInterval error, bool isExact) :
        vertical_(vertical), error_(error), isExact_(isExact) {}
    bool operator<(DependencyCandidate const& other) const;
    
    //TODO: implement if used
    //bool operator==(DependencyCandidate const& other) const;

    bool isExact() const { return isExact_ && error_.isPoint(); }

    //TODO: check if these comparators are used in a right way
    //using comparator = std::function<bool (DependencyCandidate const &, DependencyCandidate const &)>;
    static bool arityComparator(DependencyCandidate const &, DependencyCandidate const &);
    static bool minErrorComparator(DependencyCandidate const &, DependencyCandidate const &);
    //static bool maxErrorComparator(DependencyCandidate const &, DependencyCandidate const &);
    static bool fullArityErrorComparator(DependencyCandidate const &, DependencyCandidate const &);
    static bool fullErrorArityComparator(DependencyCandidate const &, DependencyCandidate const &);
    explicit operator std::string() const
        { return "Candidate " + static_cast<std::string>(*vertical_) + static_cast<std::string>(error_); }

    friend std::ostream& operator<< (std::ostream&, DependencyCandidate const&);
};
