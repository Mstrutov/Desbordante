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


    std::string toIndicesString() const { return lhs_.toIndicesString() + " - " + rhs_.toIndicesString(); }
    std::string toString() const {
        return lhs_.toString() + "~>" + rhs_.toString()
            + boost::lexical_cast<std::string>(error_)
            + boost::lexical_cast<std::string>(score_); }

    double getError() const { return error_; }
    int getArity() const { return lhs_.getColumns().size(); }
};