#include "relational_schema.h"

#include <memory>
#include <utility>

#include "vertical_map.h"

namespace {
auto MakeColumns(RelationalSchema* schema, std::vector<std::string> column_names) {
    std::vector<Column> columns;
    std::size_t const number_of_columns = column_names.size();
    columns.reserve(number_of_columns);
    for (model::ColumnIndex i = 0; i < number_of_columns; ++i) {
        columns.emplace_back(schema, column_names[i], i);
    }
    return columns;
}
}  // namespace

RelationalSchema::RelationalSchema(std::string name, std::vector<std::string> column_names)
    : columns_(MakeColumns(this, std::move(column_names))),
      name_(std::move(name)),
      empty_vertical_(std::make_unique<Vertical>(this, boost::dynamic_bitset<>(columns_.size()))) {}

std::unique_ptr<RelationalSchema> RelationalSchema::CreateFrom(model::IDatasetStream& table) {
    std::vector<std::string> column_names;
    std::size_t const number_of_columns = table.GetNumberOfColumns();
    column_names.reserve(number_of_columns);
    for (model::ColumnIndex i = 0; i < number_of_columns; ++i) {
        column_names.push_back(table.GetColumnName(i));
    }
    return std::make_unique<RelationalSchema>(table.GetRelationName(), std::move(column_names));
}

Vertical RelationalSchema::GetVertical(boost::dynamic_bitset<> indices) const {
    return Vertical(this, std::move(indices));
}

Column const& RelationalSchema::GetColumn(std::string const& col_name) const {
    auto found_entry_iterator =
            std::find_if(columns_.begin(), columns_.end(),
                         [&col_name](auto& column) { return column.name_ == col_name; });
    if (found_entry_iterator != columns_.end()) return *found_entry_iterator;

    throw std::invalid_argument("Couldn't match column name \'" + col_name +
                                "\' to any of the schema's column names");
}

Column const& RelationalSchema::GetColumn(size_t index) const {
    return columns_.at(index);
}

size_t RelationalSchema::GetNumColumns() const {
    return columns_.size();
}

// TODO: critical part - consider optimization
// TODO: list -> vector as list doesn't have RAIterators therefore can't be sorted
std::unordered_set<Vertical> RelationalSchema::CalculateHittingSet(
        std::vector<Vertical> verticals,
        boost::optional<std::function<bool(Vertical const&)>> pruning_function) const {
    std::sort(verticals.begin(), verticals.end(), [](auto& vertical1, auto& vertical2) {
        return vertical1.GetArity() < vertical2.GetArity();
    });
    model::VerticalMap<Vertical> consolidated_verticals(this);

    model::VerticalMap<Vertical> hitting_set(this);
    hitting_set.Put(*empty_vertical_, Vertical::EmptyVertical(this));

    for (auto& vertical : verticals) {
        if (consolidated_verticals.GetAnySubsetEntry(vertical).second != nullptr) {
            continue;
        }
        // TODO: костыль, тк VerticalMap хранит unique_ptr - лишнее копирование
        consolidated_verticals.Put(vertical, std::make_unique<Vertical>(vertical));

        auto invalid_hitting_set_members = hitting_set.GetSubsetKeys(vertical.Invert());
        std::sort(invalid_hitting_set_members.begin(), invalid_hitting_set_members.end(),
                  [](auto& vertical1, auto& vertical2) {
                      return vertical1.GetArity() < vertical2.GetArity();
                  });

        for (auto& invalid_hitting_set_member : invalid_hitting_set_members) {
            hitting_set.Remove(invalid_hitting_set_member);
        }

        for (auto& invalid_member : invalid_hitting_set_members) {
            for (size_t corrective_column_index = vertical.GetColumnIndices().find_first();
                 corrective_column_index != boost::dynamic_bitset<>::npos;
                 corrective_column_index =
                         vertical.GetColumnIndices().find_next(corrective_column_index)) {
                auto corrective_column = GetColumn(corrective_column_index);
                auto corrected_member =
                        invalid_member.Union(static_cast<Vertical>(corrective_column));

                if (hitting_set.GetAnySubsetEntry(corrected_member).second == nullptr) {
                    if (pruning_function) {
                        bool is_pruned = (*pruning_function)(corrected_member);
                        if (is_pruned) {
                            continue;
                        }
                    }
                    hitting_set.Put(corrected_member, std::make_unique<Vertical>(corrected_member));
                }
            }
        }
        if (hitting_set.IsEmpty()) break;
    }
    return hitting_set.KeySet();
}

RelationalSchema::~RelationalSchema() = default;
