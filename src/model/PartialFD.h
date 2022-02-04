#pragma once

#include <boost/lexical_cast.hpp>
#include "Column.h"
#include "Vertical.h"

class PartialFD {
public:
    double error_;
    Vertical lhs_;
    Column rhs_;
    double score_;

    PartialFD(Vertical lhs, Column rhs, double error, double score)
        : error_(error), lhs_(std::move(lhs)), rhs_(std::move(rhs)), score_(score) {}


    std::string ToIndicesString() const { return lhs_.ToIndicesString() + " - " + rhs_.ToIndicesString(); }
    std::string ToString() const {
        return lhs_.ToString() + "~>" + rhs_.ToString()
            + boost::lexical_cast<std::string>(error_)
            + boost::lexical_cast<std::string>(score_); }

    double GetError() const { return error_; }
    int GetArity() const { return lhs_.GetColumns().size(); }
};