#pragma once

#include <boost/lexical_cast.hpp>
#include "Column.h"
#include "Vertical.h"

class PartialKey {
public:
    double error_;
    std::shared_ptr<Vertical> vertical_;
    double score_;

    PartialKey(std::shared_ptr<Vertical> vertical, double error, double score) : error_(error), vertical_(std::move(vertical)), score_(score) {}

    std::string toString() const { return vertical_->toString() + "~>" + boost::lexical_cast<std::string>(error_) + boost::lexical_cast<std::string>(score_); }
    //double getError() const { return error_; }
    //int getArity() const { return lhs_->getColumns().size(); }
};