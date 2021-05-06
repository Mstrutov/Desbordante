//
// Created by alex on 01.05.2021.
//

#pragma once

#include "Vertical.h"

struct custom_comparator {
    bool operator()(shared_ptr<Vertical> const &v1, shared_ptr<Vertical> const &v2) const {
        return *v1 == *v2;
    }
};