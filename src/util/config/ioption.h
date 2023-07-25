#pragma once

#include <memory>
#include <string_view>
#include <typeindex>

#include "boost/any.hpp"

namespace util::config {

class IOption {
public:
    virtual std::vector<std::string_view> Set(boost::any const& value) = 0;
    virtual void Unset() = 0;
    [[nodiscard]] virtual bool IsSet() const = 0;
    [[nodiscard]] virtual std::string_view GetName() const = 0;
    [[nodiscard]] virtual std::string_view GetDescription() const = 0;
    [[nodiscard]] virtual std::type_index GetTypeIndex() const = 0;
    [[nodiscard]] virtual std::unique_ptr<IOption> MoveToHeap() = 0;
    virtual ~IOption() = default;
};

}  // namespace util::config
