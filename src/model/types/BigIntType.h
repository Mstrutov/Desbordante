#pragma once

#include "StringType.h"

namespace model {

class BigIntType : public StringType {
public:
    BigIntType() : StringType(TypeId::kBigInt){};
};

}  // namespace model
