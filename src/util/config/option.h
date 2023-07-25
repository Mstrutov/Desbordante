#pragma once

#include <cassert>
#include <functional>
#include <typeindex>
#include <typeinfo>

#include <boost/any.hpp>

#include "util/config/ioption.h"

namespace util::config {

// This class is responsible for configuration values. It is provided with a
// pointer to a field of an algorithm object, option name, and option
// description.
// It can also normalize the provided value and check it for correctness.
// Some options may have a default value, which can be provided as well.
// However, some defaults may depend on internal state of an algorithms, so a
// function that returns that default value may be passed instead.
template <typename T>
class Option : public IOption {
public:
    using DefaultFunc = std::function<T()>;
    using ValueCheckFunc = std::function<void(T const &)>;
    using CondCheckFunc = std::function<bool(T const &val)>;
    using OptCondVector = std::vector<std::pair<CondCheckFunc, std::vector<std::string_view>>>;
    using NormalizeFunc = std::function<void(T &)>;

    Option(T *value_ptr, std::string_view name, std::string_view description,
           DefaultFunc default_func = nullptr)
        : value_ptr_(value_ptr),
          name_(name),
          description_(description),
          default_func_(std::move(default_func)) {}

    Option(T *value_ptr, std::string_view name, std::string_view description, T default_value)
        : Option(value_ptr, name, description,
                 [default_value = std::move(default_value)]() { return default_value; }) {}

    std::vector<std::string_view> Set(boost::any const &value) override;

    void Unset() override {
        is_set_ = false;
    }

    [[nodiscard]] std::string_view GetName() const override {
        return name_;
    }

    [[nodiscard]] std::string_view GetDescription() const override {
        return description_;
    }

    [[nodiscard]] std::type_index GetTypeIndex() const override {
        return typeid(T);
    }

    [[nodiscard]] bool IsSet() const override {
        return is_set_;
    }

    std::unique_ptr<IOption> MoveToHeap() override {
        return std::make_unique<Option>(std::move(*this));
    }

    Option &&SetValueCheck(ValueCheckFunc value_check) {
        assert(!value_check_);
        value_check_ = std::move(value_check);
        return std::move(*this);
    }

    // Some options may become required depending on this option's value and/or
    // their correctness may depend on the value. If that is the case, a vector
    // of {predicate, option names} pairs may be provided to this method.
    // If after a predicate in one of the pairs holds for the provided value,
    // then the corresponding option names will be returned, which will then
    // become required options of the algorithm.
    // Order matters: predicates are checked first-to-last, and `Set` returns
    // only the names of the first pair where the predicate holds. An empty
    // predicate is equivalent to an always-true predicate, and thus must
    // always be last.
    Option &&SetConditionalOpts(OptCondVector opt_cond) {
        assert(opt_cond_.empty());
        assert(!opt_cond.empty());
        opt_cond_ = std::move(opt_cond);
        return std::move(*this);
    }

    Option &&SetNormalizeFunc(NormalizeFunc normalize_func) {
        assert(!normalize_func_);
        normalize_func_ = std::move(normalize_func);
        return std::move(*this);
    }

private:
    T ConvertValue(boost::any const &value) const;

    bool is_set_ = false;
    T *value_ptr_;
    std::string_view name_;
    std::string_view description_;
    DefaultFunc default_func_;
    ValueCheckFunc value_check_{};
    OptCondVector opt_cond_{};
    NormalizeFunc normalize_func_{};
};

template <typename T>
std::vector<std::string_view> Option<T>::Set(boost::any const &value) {
    assert(!is_set_);
    T converted_value = ConvertValue(value);
    if (normalize_func_) normalize_func_(converted_value);
    if (value_check_) value_check_(converted_value);

    assert(value_ptr_ != nullptr);
    is_set_ = true;
    std::vector<std::string_view> new_options{};
    for (auto const &[cond, opts] : opt_cond_) {
        if (!cond || cond(converted_value)) {
            new_options = opts;
            break;
        }
    }
    *value_ptr_ = std::move(converted_value);
    return new_options;
}

template <typename T>
T Option<T>::ConvertValue(boost::any const &value) const {
    std::string const no_value_no_default =
            std::string("No value was provided to an option without a default value (") +
            GetName().data() + ")";
    if (value.empty()) {
        if (!default_func_) throw std::invalid_argument(no_value_no_default);
        return default_func_();
    } else {
        if (value.type() != typeid(T)) throw std::invalid_argument("Incorrect option value type");
        return boost::any_cast<T>(value);
    }
}

}  // namespace util::config
