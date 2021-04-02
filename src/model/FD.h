#pragma once

#include <string>

#include "Column.h"
#include "Vertical.h"

class FD {
private:
    Vertical lhs_;
    Column rhs_;

public:
    FD(Vertical const& lhs, Column const& rhs) : lhs_(lhs), rhs_(rhs) {}

    std::string toJSONString() const {
        return "{lhs: " + lhs_.toIndicesString() + ", rhs: " + rhs_.toIndicesString() + "}";
    }

    Vertical const& getLhs() const { return lhs_; }
    Column const& getRhs() const { return rhs_; }

    // unsigned int fletcher16() const;
};