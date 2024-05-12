#include "config/custom_random/option.h"

#include <utility>

#include "config/names_and_descriptions.h"

namespace config {
extern CommonOption<CustomRandomFlagType> const kCustomRandomFlagOpt{
        names::kCustomRandom, descriptions::kDCustomRandom, std::make_pair(false, 47)
        // names::kCustomRandom, descriptions::kDCustomRandom, {false, 47}
};
}  // namespace config
