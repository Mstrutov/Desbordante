#pragma once

#include <regex>
#include <string>
#include <vector>

#include "AbstractColumnData.h"
#include "RelationData.h"
#include "types/Types.h"

namespace model {

class TypedColumnData : public model::AbstractColumnData {
private:
    std::unique_ptr<Type const> type_;
    unsigned int rows_num_;
    unsigned int nulls_num_;
    unsigned int empties_num_;
    std::unique_ptr<std::byte[]> buffer_;
    std::vector<std::byte const*> data_;
    /* For non-mixed type only */
    std::unordered_set<unsigned> nulls_;
    std::unordered_set<unsigned> empties_;

    TypedColumnData(Column const* column, std::unique_ptr<Type const> type,
                    unsigned int const rows_num, unsigned int nulls_num, unsigned int empties_num,
                    std::unique_ptr<std::byte[]> buffer, std::vector<std::byte const*> data,
                    std::unordered_set<unsigned> nulls,
                    std::unordered_set<unsigned> empties) noexcept
        : AbstractColumnData(column),
          type_(std::move(type)),
          rows_num_(rows_num),
          nulls_num_(nulls_num),
          empties_num_(empties_num),
          buffer_(std::move(buffer)),
          data_(std::move(data)),
          nulls_(std::move(nulls)),
          empties_(std::move(empties)) {}

    friend class TypedColumnDataFactory;

public:
    TypedColumnData(TypedColumnData const& other) = delete;
    TypedColumnData& operator=(TypedColumnData const& other) = delete;
    TypedColumnData(TypedColumnData&& other) noexcept = default;
    TypedColumnData& operator=(TypedColumnData&& other) noexcept = default;

    ~TypedColumnData() {
        if (type_ == nullptr) {
            return;
        }

        MixedType const* mixed = GetIfMixed();
        for (size_t i = 0; i != data_.size(); ++i) {
            if (GetValueTypeId(i) == +TypeId::kString || GetValueTypeId(i) == +TypeId::kBigInt) {
                std::byte const* value = (mixed) ? mixed->RetrieveValue(data_[i]) : data_[i];
                StringType::Destruct(value);
            }
        }
    }

    TypeId GetTypeId() const noexcept {
        return type_->GetTypeId();
    }

    Type const& GetType() const noexcept {
        return *type_;
    }

    std::vector<std::byte const*> const& GetData() const noexcept {
        return data_;
    }

    unsigned int GetNumNulls() const noexcept {
        return nulls_num_;
    }

    unsigned int GetNumEmpties() const noexcept {
        return empties_num_;
    }

    unsigned int GetNumRows() const noexcept {
        return rows_num_;
    }

    bool IsNull(unsigned int index) const noexcept {
        MixedType const* mixed = GetIfMixed();
        if (mixed != nullptr) {
            return mixed->RetrieveTypeId(data_[index]) == +TypeId::kNull;
        } else {
            return nulls_.find(index) != nulls_.end();
        }
    }

    bool IsEmpty(unsigned int index) const noexcept {
        MixedType const* mixed = GetIfMixed();
        if (mixed != nullptr) {
            return mixed->RetrieveTypeId(data_[index]) == +TypeId::kEmpty;
        } else {
            return empties_.find(index) != empties_.end();
        }
    }

    TypeId GetValueTypeId(unsigned int index) const noexcept {
        TypeId const type_id = type_->GetTypeId();
        if (type_id == +TypeId::kMixed) {
            return static_cast<MixedType const*>(type_.get())->RetrieveTypeId(data_[index]);
        }

        if (IsNull(index)) {
            return TypeId::kNull;
        }
        if (IsEmpty(index)) {
            return TypeId::kEmpty;
        }

        return type_id;
    }

    bool IsNumeric() const noexcept {
        TypeId type_id = GetTypeId();
        return type_id == +TypeId::kInt || /* type_id == +ColumnTypeId::kBigInt || */
               type_id == +TypeId::kDouble;
    }

    bool IsMixed() const noexcept {
        return GetTypeId() == +TypeId::kMixed;
    }

    MixedType const* GetIfMixed() const noexcept {
        return dynamic_cast<MixedType const*>(type_.get());
    }

    std::string ToString() const final {
        return "Data for column " + column_->ToString() + " of type " + GetTypeId()._to_string();
    }
};

class TypedColumnDataFactory {
private:
    using TypeMap = std::unordered_map<TypeId, std::unordered_set<unsigned>>;
    using TypeIdToType = std::unordered_map<TypeId, std::unique_ptr<Type>>;

    Column const* column_;
    std::vector<std::string> unparsed_;
    bool is_null_equal_null_;

    inline static const std::unordered_map<TypeId, std::regex> type_id_to_regex_ = {
        {TypeId::kInt, std::regex(R"(^(\+|-)?\d{1,19}$)")},
        {TypeId::kBigInt, std::regex(R"(^(\+|-)?\d{20,}$)")},
        {TypeId::kDouble, std::regex(R"(^(\+|-)?\d+\.\d*$)")},
        {TypeId::kNull, std::regex(Null::kValue.data())},
        {TypeId::kEmpty, std::regex(R"(^$)")}};

    size_t CalculateMixedBufSize(MixedType const* mixed, TypeIdToType const& tid_to_map,
                                 TypeMap const& type_map) const noexcept;
    std::vector<TypeId> GetTypesLayout(TypeMap const& tm) const;
    TypeIdToType MapTypeIdsToTypes(TypeMap const& tm) const;
    TypeMap CreateTypeMap() const;
    TypedColumnData CreateMixedFromTypeMap(std::unique_ptr<Type const> type, TypeMap type_map);
    TypedColumnData CreateConcreteFromTypeMap(std::unique_ptr<Type const> type, TypeMap type_map);
    TypedColumnData CreateFromTypeMap(std::unique_ptr<Type const> type, TypeMap type_map);
    TypedColumnData CreateFrom();

    TypedColumnDataFactory(Column const* col, std::vector<std::string> unparsed,
                           bool is_null_equal_null)
        : column_(col), unparsed_(std::move(unparsed)), is_null_equal_null_(is_null_equal_null) {}

public:
    static TypedColumnData CreateFrom(Column const* col, std::vector<std::string> unparsed,
                                      bool is_null_equal_null) {
        TypedColumnDataFactory f(col, std::move(unparsed), is_null_equal_null);
        return f.CreateFrom();
    }
};

}  // namespace model
