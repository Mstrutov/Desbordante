#pragma once

#include "Type.h"

namespace model {

class IMetrizableType : public Type {
public:
    using Type::Type;

    virtual double Dist(std::byte const* l, std::byte const* r) const = 0;
};

}  // namespace model
