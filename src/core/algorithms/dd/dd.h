#pragma once

#include <list>
#include <sstream>
#include <string>
#include <vector>

namespace model {

using DF = std::vector<std::pair<double, double>>;

struct DD {
    DF lhs;
    DF rhs;

    DD() = default;

    DD(DF const& left, DF const& right) : lhs(left), rhs(right) {}
};

struct DFConstraint {
    std::string column_name;
    double lower_bound;
    double upper_bound;

    DFConstraint() = default;

    DFConstraint(std::string name, double lower, double upper)
        : column_name(name), lower_bound(lower), upper_bound(upper) {}

    bool operator<(DFConstraint const& other) const {
        return column_name < other.column_name ||
               (column_name == other.column_name && lower_bound < other.lower_bound) ||
               (column_name == other.column_name && lower_bound == other.lower_bound &&
                upper_bound < other.upper_bound);
    }

    std::string ToString() const {
        std::stringstream s;
        s << column_name << " [" << lower_bound << ", " << upper_bound << "]";
        return s.str();
    }
};

struct DDString {
    std::list<DFConstraint> left;
    std::list<DFConstraint> right;

    std::string ToString() const {
        return DFToString(left) + " -> " + DFToString(right);
    }

    std::string DFToString(std::list<DFConstraint> const& df) const {
        std::stringstream s;
        bool has_constraints = false;
        for (auto constraint : df) {
            if (has_constraints)
                s << " ; ";
            else
                has_constraints = true;
            s << constraint.ToString();
        }
        return s.str();
    }
};

}  // namespace model
