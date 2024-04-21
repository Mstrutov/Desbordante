#pragma once
#include "builtin.h"
#include "date_type.h"
#include "string_type.h"

namespace model {
inline void Destruct(TypeId type_id, std::byte const* value) {
    if (type_id == +TypeId::kString || type_id == +TypeId::kBigInt) {
        StringType::Destruct(value);
    }
    if (type_id == +TypeId::kDate) {
        DateType::Destruct(value);
    }
}
}  // namespace model