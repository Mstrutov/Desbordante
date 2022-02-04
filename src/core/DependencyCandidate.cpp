#include <boost/dynamic_bitset.hpp>

#include "DependencyCandidate.h"

// TODO: these methods are used in priority_queues, where operator> is needed. (>) !<=> (>=) due to strict weak ordering
bool DependencyCandidate::ArityComparator(DependencyCandidate const &dc1, DependencyCandidate const &dc2) {
    if (dc1.vertical_.GetArity() > dc2.vertical_.GetArity())
        return true;
    else if (dc1.vertical_.GetArity() == dc2.vertical_.GetArity())
        return (dc1.error_.GetMean() > dc2.error_.GetMean());
    return false;
}

bool DependencyCandidate::MinErrorComparator(DependencyCandidate const &dc1, DependencyCandidate const &dc2) {
    if (dc1.error_.GetMin() > dc2.error_.GetMin())
        return true;
    else if (dc1.error_.GetMin() == dc2.error_.GetMin())
        return (dc1.vertical_.GetArity() > dc2.vertical_.GetArity());
    return false;
}

bool DependencyCandidate::FullErrorArityComparator(DependencyCandidate const &dc1, DependencyCandidate const &dc2) {
    return dc1 < dc2;
}

//TODO: perhaps in Java code they meant to implement error -> arity -> lexicographical? Check usages
bool DependencyCandidate::FullArityErrorComparator(DependencyCandidate const &dc1, DependencyCandidate const &dc2) {
    if (dc1.error_.GetMean() < dc2.error_.GetMean())
        return true;
    else if (dc1.error_.GetMean() == dc2.error_.GetMean()) {
        if (dc1.vertical_.GetArity() < dc2.vertical_.GetArity())
            return true;
        else if (dc1.vertical_.GetArity() == dc2.vertical_.GetArity()) {
            boost::dynamic_bitset<> dc1_cols = dc1.vertical_.GetColumnIndices();
            boost::dynamic_bitset<> dc2_cols = dc2.vertical_.GetColumnIndices();

            for (size_t a = dc1_cols.find_first(), b = dc2_cols.find_first();
                 a < dc1_cols.size();
                 a = dc1_cols.find_next(a), b = dc2_cols.find_next(b))
                if (a != b)
                    return (a < b);
        }
    }
    return false;
}

bool DependencyCandidate::operator<(DependencyCandidate const & other) const {
    if (error_.GetMean() < other.error_.GetMean())
        return true;
    else if (error_.GetMean() == other.error_.GetMean()) {
        if (vertical_.GetArity() < other.vertical_.GetArity())
            return true;
        else if (vertical_.GetArity() == other.vertical_.GetArity()) {
            boost::dynamic_bitset<> dc1_cols = vertical_.GetColumnIndices();
            boost::dynamic_bitset<> dc2_cols = other.vertical_.GetColumnIndices();

            for (size_t a = dc1_cols.find_first(), b = dc2_cols.find_first();
                 a < dc1_cols.size();
                 a = dc1_cols.find_next(a), b = dc2_cols.find_next(b))
                if (a != b)
                    return (a < b);
        }
    }
    return false;
}

std::ostream& operator<< (std::ostream& ofs, DependencyCandidate const& dependency_candidate) {
    return ofs << static_cast<std::string>(dependency_candidate);
}

/*bool DependencyCandidate::operator==(DependencyCandidate const & other) const {
    return (*vertical_ == *(other.vertical_) && )
}*/