#include <boost/dynamic_bitset.hpp>

#include "DependencyCandidate.h"

// TODO: these methods are used in priority_queues, where operator> is needed. (>) !<=> (>=) due to strict weak ordering
bool DependencyCandidate::arityComparator(DependencyCandidate const & dc1, DependencyCandidate const & dc2) {
    if (dc1.vertical_->getArity() > dc2.vertical_->getArity())
        return true;
    else if (dc1.vertical_->getArity() == dc2.vertical_->getArity())
        return (dc1.error_.getMean() > dc2.error_.getMean());
    return false;
}

bool DependencyCandidate::minErrorComparator(DependencyCandidate const & dc1, DependencyCandidate const & dc2) {
    if (dc1.error_.getMin() > dc2.error_.getMin())
        return true;
    else if (dc1.error_.getMin() == dc2.error_.getMin())
        return (dc1.vertical_->getArity() > dc2.vertical_->getArity());
    return false;
}

bool DependencyCandidate::fullErrorArityComparator(DependencyCandidate const & dc1, DependencyCandidate const & dc2) {
    return dc1 < dc2;
}

//TODO: perhaps in Java code they meant to implement error -> arity -> lexicographical? Check usages
bool DependencyCandidate::fullArityErrorComparator(DependencyCandidate const & dc1, DependencyCandidate const & dc2) {
    if (dc1.error_.getMean() < dc2.error_.getMean())
        return true;
    else if (dc1.error_.getMean() == dc2.error_.getMean()) {
        if (dc1.vertical_->getArity() < dc2.vertical_->getArity())
            return true;
        else if (dc1.vertical_->getArity() == dc2.vertical_->getArity()) {
            boost::dynamic_bitset<> dc1Cols = dc1.vertical_->getColumnIndices();
            boost::dynamic_bitset<> dc2Cols = dc2.vertical_->getColumnIndices();

            for (size_t a = dc1Cols.find_first(), b = dc2Cols.find_first();
                 a < dc1Cols.size();
                 a = dc1Cols.find_next(a), b = dc2Cols.find_next(b))
                if (a != b)
                    return (a < b);
        }
    }
    return false;
}

bool DependencyCandidate::operator<(DependencyCandidate const & other) const {
    if (error_.getMean() < other.error_.getMean())
        return true;
    else if (error_.getMean() == other.error_.getMean()) {
        if (vertical_->getArity() < other.vertical_->getArity())
            return true;
        else if (vertical_->getArity() == other.vertical_->getArity()) {
            boost::dynamic_bitset<> dc1Cols = vertical_->getColumnIndices();
            boost::dynamic_bitset<> dc2Cols = other.vertical_->getColumnIndices();

            for (size_t a = dc1Cols.find_first(), b = dc2Cols.find_first();
                 a < dc1Cols.size();
                 a = dc1Cols.find_next(a), b = dc2Cols.find_next(b))
                if (a != b)
                    return (a < b);
        }
    }
    return false;
}

std::ostream& operator<< (std::ostream& ofs, DependencyCandidate const& dependencyCandidate) {
    return ofs << static_cast<std::string>(dependencyCandidate);
}

/*bool DependencyCandidate::operator==(DependencyCandidate const & other) const {
    return (*vertical_ == *(other.vertical_) && )
}*/