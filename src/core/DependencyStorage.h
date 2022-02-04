#pragma once

#include <list>
#include "model/PartialFD.h"
#include "model/PartialKey.h"

class DependencyStorage {
public:
    std::list<PartialFD> discovered_fds_;
    std::list<PartialKey> discovered_uccs_;
};