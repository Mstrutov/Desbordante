#pragma once

#include "dc/predicate_provider.h"

namespace model {

class PredicateBuilder {
public:
    PredicateBuilder() {
        PredicateProvider::CreateInstance();
    }
};
}  // namespace model
