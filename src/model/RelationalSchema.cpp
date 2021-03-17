#include <memory>
#include <utility>

#include "RelationalSchema.h"
#include "Vertical.h"
#include "VerticalMap.h"

RelationalSchema::RelationalSchema(std::string name, bool isNullEqNull) :
        columns(),
        name(std::move(name)),
        isNullEqNull(isNullEqNull),
        emptyVertical() {
    init();
}

void RelationalSchema::init() {
    emptyVertical = Vertical::emptyVertical(this);
}

//TODO: В оригинале тут что-то непонятное + приходится пересоздавать emptyVertical -- тут
//должен быть unique_ptr, тк создаём в остальных случаях новую вершину и выдаём наружу с овнершипом
Vertical RelationalSchema::getVertical(boost::dynamic_bitset<> indices) const {
    if (indices.empty()) return *Vertical::emptyVertical(this);

    if (indices.count() == 1){
        return Vertical(this, std::move(indices));
    }
    return Vertical(this, std::move(indices));
}

Column const* RelationalSchema::getColumn(const std::string &colName) const {
    auto foundEntryIterator = std::find_if(columns.begin(), columns.end(),
                                           [&colName](auto& column) { return column->name == colName; });
    if (foundEntryIterator != columns.end()) return foundEntryIterator->get();

    throw std::invalid_argument("Couldn't match column name \'"
        + colName
        + "\' to any of the schema's column names");
}

Column const* RelationalSchema::getColumn(int index) const {
    return columns.at(index).get();
}

void RelationalSchema::appendColumn(const std::string& colName) {
    columns.push_back(std::make_unique<Column>(this, colName, columns.size()));
}

void RelationalSchema::appendColumn(Column column) {
    columns.push_back(std::make_unique<Column>(std::move(column)));
}

int RelationalSchema::getNumColumns() const {
    return columns.size();
}

bool RelationalSchema::isNullEqualNull() const { return isNullEqNull; }

// TODO: critical part - consider optimization
// TODO: list -> vector as list doesn't have RAIterators therefore can't be sorted
std::unordered_set<Vertical> RelationalSchema::calculateHittingSet(
        std::vector<Vertical> verticals, boost::optional<std::function<bool (Vertical const&)>> pruningFunction) const {
    std::sort(verticals.begin(), verticals.end(),
              [](auto vertical1, auto vertical2) { return vertical1->getArity() < vertical2->getArity(); });
    VerticalMap<Vertical> consolidatedVerticals(this);

    VerticalMap<Vertical> hittingSet(this);
    hittingSet.put(*emptyVertical, Vertical::emptyVertical(this));

    for (auto& vertical : verticals) {
        if (consolidatedVerticals.getAnySubsetEntry(vertical).second != nullptr) {
            continue;
        }
        // TODO: костыль, тк VerticalMap хранит unique_ptr - лишнее копирование
        consolidatedVerticals.put(vertical, std::make_unique<Vertical>(vertical));

        auto invalidHittingSetMembers = hittingSet.getSubsetKeys(vertical.invert());
        std::sort(invalidHittingSetMembers.begin(), invalidHittingSetMembers.end(),
                [](auto& vertical1, auto& vertical2) { return vertical1.getArity() < vertical2.getArity(); });

        for (auto& invalidHittingSetMember : invalidHittingSetMembers) {
            hittingSet.remove(invalidHittingSetMember);
        }

        for (auto& invalidMember : invalidHittingSetMembers) {
            for (size_t correctiveColumnIndex = vertical.getColumnIndices().find_first();
                 correctiveColumnIndex != boost::dynamic_bitset<>::npos;
                 correctiveColumnIndex = vertical.getColumnIndices().find_next(correctiveColumnIndex)) {

                auto correctiveColumn = *getColumn(correctiveColumnIndex);
                auto correctedMember = invalidMember.Union(static_cast<Vertical>(correctiveColumn));

                if (hittingSet.getAnySubsetEntry(correctedMember).second == nullptr) {
                    if (pruningFunction) {
                        bool isPruned = (*pruningFunction)(correctedMember);
                        if (isPruned) {
                            continue;
                        }
                    }
                    hittingSet.put(correctedMember, std::make_unique<Vertical>(correctedMember));
                }
            }
        }
        if (hittingSet.isEmpty()) break;
    }
    return hittingSet.keySet();
}