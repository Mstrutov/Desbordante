#include "config/tabular_data/input_table/option.h"

#include "config/names_and_descriptions.h"

namespace config {
using names::kTable, descriptions::kDTable;
extern CommonOption<InputTable> const TableOpt = {
        kTable, kDTable, std::nullopt, nullptr, [](InputTable const& table) {
            // Shouldn't happen normally.
            if (table == nullptr) throw ConfigurationError("Table must be a valid pointer!");
        }};
}  // namespace config
