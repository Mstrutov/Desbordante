#include "TypedColumnData.h"

#include "CreateType.h"

namespace model {

TypedColumnDataFactory::TypeMap TypedColumnDataFactory::CreateTypeMap() const {
    TypeMap type_map;
    auto const match = [&type_map](std::string const& val, size_t const row) {
        bool matched = false;
        for (auto const& [type_id, regex] : type_id_to_regex_) {
            if (std::regex_match(val, regex)) {
                type_map[type_id].insert(row);
                matched = true;
                break;
            }
        }
        if (!matched) {
            type_map[TypeId::kString].insert(row);
        }
    };

    for (size_t i = 0; i != unparsed_.size(); ++i) {
        match(unparsed_[i], i);
    }

    if (type_map.count(TypeId::kBigInt) && type_map.count(TypeId::kInt)) {
        std::unordered_set<unsigned>& big_ints = type_map[TypeId::kBigInt];
        std::unordered_set<unsigned> ints = std::move(type_map.extract(TypeId::kInt).mapped());
        big_ints.insert(ints.begin(), ints.end());
    }

    return type_map;
}

std::vector<TypeId> TypedColumnDataFactory::GetTypesLayout(TypeMap const& tm) const {
    std::vector<TypeId> types_layout(unparsed_.size(), TypeId::kString);

    for (auto const& [type_id, indices] : tm) {
        for (unsigned const index : indices) {
            types_layout[index] = type_id;
        }
    }

    return types_layout;
}

TypedColumnDataFactory::TypeIdToType TypedColumnDataFactory::MapTypeIdsToTypes(
    TypeMap const& tm) const {
    std::unordered_map<TypeId, std::unique_ptr<Type>> type_id_to_type;
    for (auto const& [type_id, indices] : tm) {
        type_id_to_type.emplace(type_id, CreateType(type_id, is_null_equal_null_));
    }
    return type_id_to_type;
}

size_t TypedColumnDataFactory::CalculateMixedBufSize(MixedType const* mixed,
                                                     TypeIdToType const& type_id_to_type,
                                                     TypeMap const& type_map) const noexcept {
    size_t result = 0;
    for (auto const& [type_id, indices] : type_map) {
        result += (mixed->GetMixedValueSize(type_id_to_type.at(type_id).get())) * indices.size();
    }
    return result;
}

TypedColumnData TypedColumnDataFactory::CreateMixedFromTypeMap(std::unique_ptr<Type const> type,
                                                               TypeMap type_map) {
    assert(type->GetTypeId() == +TypeId::kMixed);
    MixedType const* mixed_type = static_cast<MixedType const*>(type.get());
    std::vector<std::byte const*> data;
    data.reserve(unparsed_.size());

    std::unordered_set<unsigned> const& nulls = type_map[TypeId::kNull];
    std::unordered_set<unsigned> const& empties = type_map[TypeId::kEmpty];
    unsigned int const rows_num = unparsed_.size();
    unsigned int const nulls_num = nulls.size();
    unsigned int const empties_num = empties.size();

    TypeIdToType type_id_to_type = MapTypeIdsToTypes(type_map);
    size_t const buf_size = CalculateMixedBufSize(mixed_type, type_id_to_type, type_map);
    std::unique_ptr<std::byte[]> buf(new std::byte[buf_size]);
    std::vector<TypeId> types_layout = GetTypesLayout(type_map);
    type_map.clear(); /* type_map is no longer needed, so saving space */

    std::byte* next = buf.get();
    for (size_t i = 0; i != types_layout.size(); ++i) {
        TypeId const type_id = types_layout[i];
        Type const* concrete_type = type_id_to_type.at(type_id).get();
        size_t const value_size = mixed_type->GetMixedValueSize(concrete_type);

        assert(next + value_size <= buf.get() + buf_size);

        mixed_type->ValueFromStr(next, std::move(unparsed_[i]), concrete_type);

        data.push_back(next);
        next += value_size;
    }

    return TypedColumnData(column_, std::move(type), rows_num, nulls_num, empties_num,
                           std::move(buf), std::move(data), {}, {});
}

TypedColumnData TypedColumnDataFactory::CreateConcreteFromTypeMap(std::unique_ptr<Type const> type,
                                                                  TypeMap type_map) {
    TypeId const type_id = type->GetTypeId();

    if (type_id == +TypeId::kMixed) {
        /* For mixed type use CreateMixedFromTypeMap. */
        assert(0);
    }

    std::unordered_set<unsigned> nulls = std::move(type_map[TypeId::kNull]);
    std::unordered_set<unsigned> empties = std::move(type_map[TypeId::kEmpty]);
    unsigned int const rows_num = unparsed_.size();
    unsigned int const nulls_num = nulls.size();
    unsigned int const empties_num = empties.size();
    assert(rows_num >= nulls_num + empties_num);

    std::vector<std::byte const*> data(unparsed_.size());

    if (type_id == +TypeId::kUndefined) {
        return TypedColumnData(column_, std::move(type), rows_num, nulls_num, empties_num, nullptr,
                               std::move(data), std::move(nulls), std::move(empties));
    }

    std::unique_ptr<std::byte[]> buf(type->Allocate(rows_num - nulls_num - empties_num));

    unsigned buf_index = 0;
    size_t const value_size = type->GetSize();
    for (unsigned i : type_map.at(type_id)) {
        assert(buf_index <= type->GetSize() * type_map.at(type_id).size());
        std::byte* next = buf.get() + buf_index;
        type->ValueFromStr(next, std::move(unparsed_[i]));
        data[i] = next;
        buf_index += value_size;
    }

    return TypedColumnData(column_, std::move(type), rows_num, nulls_num, empties_num,
                           std::move(buf), std::move(data), std::move(nulls), std::move(empties));
}

TypedColumnData TypedColumnDataFactory::CreateFromTypeMap(std::unique_ptr<Type const> type,
                                                          TypeMap type_map) {
    if (type->GetTypeId() == +TypeId::kMixed) {
        return CreateMixedFromTypeMap(std::move(type), std::move(type_map));
    } else {
        return CreateConcreteFromTypeMap(std::move(type), std::move(type_map));
    }
}

TypedColumnData TypedColumnDataFactory::CreateFrom() {
    TypeMap type_map = CreateTypeMap();

    decltype(type_map)::node_type null_node = type_map.extract(TypeId::kNull);
    decltype(type_map)::node_type empty_node = type_map.extract(TypeId::kEmpty);

    TypeId type_id = TypeId::kMixed;

    if (type_map.empty()) {
        type_id = TypeId::kUndefined;
    } else if (type_map.size() == 1) {
        type_id = type_map.begin()->first;
    }

    type_map.insert(std::move(null_node));
    type_map.insert(std::move(empty_node));

    return CreateFromTypeMap(CreateType(type_id, is_null_equal_null_), std::move(type_map));
}

}  // namespace model
