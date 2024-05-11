#include "config/names_and_descriptions.h"
#include "config/custom_random/option.h"
#include <utility>

namespace config {
extern CommonOption<CustomRandomFlagType> const CustomRandomFlagOpt {
    names::kCustomRandom, descriptions::kDCustomRandom, std::make_pair(false, 47)
    // names::kCustomRandom, descriptions::kDCustomRandom, {false, 47}
};
}  // namespace config
