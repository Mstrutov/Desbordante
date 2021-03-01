#pragma once

#include <list>
#include "model/PartialFD.h"
#include "model/PartialKey.h"

class DependencyStorage {
public:
    std::list<PartialFD> discoveredFDs_;
    std::list<PartialKey> discoveredUCCs_;
};