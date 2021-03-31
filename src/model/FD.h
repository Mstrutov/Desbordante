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

    // unsigned int fletcher16() const;
};