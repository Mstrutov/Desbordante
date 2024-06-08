#pragma once

#include <memory>
#include <stdexcept>

namespace model {

class PredicateBuilder;

template <typename Derived>
class BaseProvider {
protected:
    static std::shared_ptr<Derived> instance_;

    BaseProvider() = default;

    // FIXME: this is wrong, there should be one class representing the algorithm,
    // which will create/clear singletons and call PredicateBuilder, PliShardBuiler,
    // SomeOtherBuilder in order..
    friend PredicateBuilder;

    static void CreateInstance() {
        if (instance_) {
            throw std::runtime_error(Derived::ClassName() +
                                     " instance is already created. "
                                     "Perhaps PredicateBuilder has not been deleted");
        }
        instance_ = std::shared_ptr<Derived>(new Derived());
    }

    static void ClearInstance() {
        if (instance_) {
            Derived::Clear();
            instance_.reset();
        }
    }

public:
    BaseProvider(BaseProvider const &) = delete;
    BaseProvider &operator=(BaseProvider const &) = delete;

    static std::shared_ptr<Derived> GetInstance() {
        if (!instance_) {
            throw std::runtime_error(Derived::ClassName() +
                                     " instance is not created. "
                                     "Perhaps PredicateBuilder has not been created");
        }
        return instance_;
    }
};

template <typename Derived>
std::shared_ptr<Derived> BaseProvider<Derived>::instance_ = nullptr;

}  // namespace model
