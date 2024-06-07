#include "config/indices/od_context.h"

#include "config/names_and_descriptions.h"
#include "indices/type.h"

namespace config {
extern CommonOption<IndicesType> const kODContextOpt{names::kODContext, descriptions::kDODContext,
                                                     IndicesType({})};
}  // namespace config
