#include "config/tabular_data/crud_operations/delete/option.h"

#include "config/names_and_descriptions.h"

namespace config {
using names::kDeleteStatements, descriptions::kDDeleteStatements;
extern CommonOption<std::vector<size_t>> const kDeleteStatementsOpt = {
        kDeleteStatements, kDDeleteStatements, std::nullopt, nullptr, nullptr};
}  // namespace config
