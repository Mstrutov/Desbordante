#pragma once

#include <sstream>

#include "BuiltIn.h"

namespace model {

class Type {
private:
    TypeId type_id_;

public:
    explicit Type(TypeId const type_id) noexcept : type_id_(type_id) {}
    virtual ~Type() = default;

    /* Operations on the type itself */

    friend inline bool operator==(Type const& l, Type const& r) {
        return l.type_id_ == r.type_id_;
    }
    friend inline bool operator!=(Type const& l, Type const& r) {
        return !(l == r);
    }

    [[nodiscard]] TypeId GetTypeId() const noexcept {
        return type_id_;
    }

    [[nodiscard]] bool IsNumeric() const noexcept {
        return type_id_ == +TypeId::kInt /* || type_id_ == +TypeId::kBigInt */ ||
               type_id_ == +TypeId::kDouble;
    }

    [[nodiscard]] std::string ToString() const {
        return type_id_._to_string();
    }

    /* Operations on values of current type */

    class Hasher {
        Type const* type_;

    public:
        explicit Hasher(Type const* type) noexcept : type_(type) {}
        size_t operator()(std::byte const* key) const {
            return type_->Hash(key);
        }
    };

    class Comparator {
        Type const* type_;

    public:
        explicit Comparator(Type const* type) noexcept : type_(type) {}
        size_t operator()(std::byte const* a, std::byte const* b) const {
            return (type_->Compare(a, b) == CompareResult::kEqual);
        }
    };

    Comparator GetComparator() const noexcept {
        return Comparator(this);
    }

    Hasher GetHasher() const noexcept {
        return Hasher(this);
    }

    void Print(std::byte const* value, std::ostream& os) const {
        os << ValueToString(value);
    }

    [[nodiscard]] virtual std::byte* Clone(std::byte const* value) const {
        size_t size = GetSize();
        auto* new_value = new std::byte[size];
        std::memcpy(new_value, value, size);
        return new_value;
    }

    [[nodiscard]] std::byte* Allocate(size_t count = 1) const {
        return new std::byte[GetSize() * count];
    }

    virtual void Free(std::byte const* val) const noexcept {
        delete[] val;
    }
    [[nodiscard]] virtual std::string ValueToString(std::byte const* value) const = 0;
    [[nodiscard]] virtual CompareResult Compare(std::byte const* l, std::byte const* r) const = 0;
    [[nodiscard]] virtual size_t Hash(std::byte const* value) const = 0;
    [[nodiscard]] virtual size_t GetSize() const = 0;
    virtual void ValueFromStr(std::byte* dest, std::string s) const = 0;

    template <typename Container>
    void DeallocateContainer(Container const& c) {
        for (std::byte const* val : c) {
            Free(val);
        }
    }

    template <typename T>
    [[nodiscard]] static T const& GetValue(std::byte const* buf) {
        return *reinterpret_cast<T const*>(buf);
    }
    template <typename T>
    [[nodiscard]] static T& GetValue(std::byte* buf) {
        return *reinterpret_cast<T*>(buf);
    }
};

}  // namespace model
