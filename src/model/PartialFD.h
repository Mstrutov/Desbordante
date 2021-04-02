#pragma once

#include <boost/lexical_cast.hpp>
#include "Column.h"
#include "Vertical.h"

class PartialFD {
public:
    double error_;
    std::shared_ptr<Vertical> lhs_;
    std::shared_ptr<Column> rhs_;
    double score_;

    PartialFD(std::shared_ptr<Vertical> lhs, std::shared_ptr<Column> rhs, double error, double score) : error_(error), lhs_(std::move(lhs)), rhs_(std::move(rhs)), score_(score) {}

    std::string toString() const { return lhs_->toString() + "~>" + rhs_->toString() + boost::lexical_cast<std::string>(error_) + boost::lexical_cast<std::string>(score_); }
    std::string toIndicesString() const { return lhs_->toIndicesString() + " - " + rhs_->toIndicesString(); }

    double getError() const { return error_; }
    int getArity() const { return lhs_->getColumns().size(); }
};