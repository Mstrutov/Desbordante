#include "ascending_od/option.h"

#include "ascending_od/type.h"
#include "config/names_and_descriptions.h"

namespace config {
extern CommonOption<AscendingODFlagType> const kAscendingODOpt{names::kAscendingOD,
                                                               descriptions::kDAscendingOD, true};
}  // namespace config
