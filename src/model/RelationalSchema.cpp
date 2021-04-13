#include <utility>

#include "RelationalSchema.h"
#include "Vertical.h"
#include "VerticalMap.h"

using namespace std;

RelationalSchema::RelationalSchema(string name, bool isNullEqNull) :
        columns(),
        name(std::move(name)),
        isNullEqNull(isNullEqNull),
        emptyVertical() {
}
// good practice is using std::make_shared instead
shared_ptr<RelationalSchema> RelationalSchema::create(string name, bool isNullEqNull) {
    auto schema = shared_ptr<RelationalSchema>(new RelationalSchema(std::move(name), isNullEqNull));
    schema->init();
    return schema;
}

// this is hard to comprehend
void RelationalSchema::init() {
    emptyVertical.reset(new Vertical(std::move(Vertical::emptyVertical(shared_from_this()))));
}

//TODO: В оригинале тут что-то непонятное
std::shared_ptr<Vertical> RelationalSchema::getVertical(dynamic_bitset<> indices) {
    if (indices.empty()) return this->emptyVertical;

    if (indices.count() == 1){
        return std::make_unique<Vertical>(static_cast<Vertical>(*this->columns[indices.find_first()]));          //TODO: TEMPORAL KOSTYL'
    }
    return std::make_unique<Vertical>(shared_from_this(), indices);
}

string RelationalSchema::getName() { return name; }

vector<shared_ptr<Column>> RelationalSchema::getColumns() { return columns; }

//TODO: assert'ы пофиксить на нормальные эксепшены
shared_ptr<Column> RelationalSchema::getColumn(const string &colName) {
    for (auto &column : columns){
        if (column->name == colName)
            return column;
    }
    assert(0);
}

shared_ptr<Column> RelationalSchema::getColumn(int index) {
    return columns[index];
}

void RelationalSchema::appendColumn(const string& colName) {
    columns.push_back(make_shared<Column>(shared_from_this(), colName, columns.size()));
}

// if you have nothing else to do: push_back through move semantics
void RelationalSchema::appendColumn(shared_ptr<Column> column) {
    columns.push_back(column);
}

int RelationalSchema::getNumColumns() const {
    return columns.size();
}

bool RelationalSchema::isNullEqualNull() { return isNullEqNull; }

// TODO: critical part - consider optimization
// TODO: list -> vector as list doesn't have RAIterators therefore can't be sorted
std::unordered_set<std::shared_ptr<Vertical>> RelationalSchema::calculateHittingSet(std::list<std::shared_ptr<Vertical>>&& verticals, boost::optional<std::function<bool (Vertical const&)>> pruningFunction) {
    using std::shared_ptr;
    //auto arityComparator = [](auto vertical1, auto vertical2) { return vertical1->getArity() < vertical2->getArity(); };
    verticals.sort([](auto vertical1, auto vertical2) { return vertical1->getArity() < vertical2->getArity(); });
    VerticalMap<shared_ptr<Vertical>> consolidatedVerticals(shared_from_this());

    VerticalMap<shared_ptr<Vertical>> hittingSet(shared_from_this());
    hittingSet.put(*emptyVertical, emptyVertical);

    for (auto vertical_ptr : verticals) {
        if (consolidatedVerticals.getAnySubsetEntry(*vertical_ptr).second != nullptr) {
            continue;
        }
        consolidatedVerticals.put(*vertical_ptr, vertical_ptr);

        auto invalidHittingSetMembers = hittingSet.getSubsetKeys(*vertical_ptr->invert());
        std::sort(invalidHittingSetMembers.begin(), invalidHittingSetMembers.end(),
                [](auto vertical1, auto vertical2) { return vertical1->getArity() < vertical2->getArity(); });

        for (auto& invalidHittingSetMember : invalidHittingSetMembers) {
            hittingSet.remove(*invalidHittingSetMember);
        }

        for (auto& invalidMember : invalidHittingSetMembers) {
            for (size_t correctiveColumnIndex = vertical_ptr->getColumnIndices().find_first();
                 correctiveColumnIndex != boost::dynamic_bitset<>::npos;
                 correctiveColumnIndex = vertical_ptr->getColumnIndices().find_next(correctiveColumnIndex)) {

                auto correctiveColumn = *getColumn(correctiveColumnIndex);
                auto correctedMember = invalidMember->Union(static_cast<Vertical>(correctiveColumn));

                if (hittingSet.getAnySubsetEntry(*correctedMember).second == nullptr) {
                    if (pruningFunction) {
                        bool isPruned = (*pruningFunction)(*correctedMember);
                        if (isPruned) {
                            continue;
                        }
                    }
                    hittingSet.put(*correctedMember, correctedMember);
                }
            }
        }
        if (hittingSet.isEmpty()) break;
    }
    return hittingSet.keySet();
}