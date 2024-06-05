#pragma once

#include <array>
#include <cstddef>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

#include <pybind11/pybind11.h>

#include "algorithms/algorithm.h"

namespace python_bindings {

namespace detail {
template <typename AlgoType>
std::string MakeDocString() {
    auto algo = AlgoType();
    std::stringstream docstring;
    docstring << "Options:\n";
    for (std::string_view option_name : algo.GetPossibleOptions()) {
        docstring << option_name << ": " << algo.GetDescription(option_name) << "\n";
    }
    return docstring.str();
}

template <typename AlgorithmType, typename BaseType>
auto RegisterAlgorithm(pybind11::module_ module, auto&& name) {
    namespace py = pybind11;
    // py::multiple_inheritance only needed for Pyro. May incur an unnecessary performance cost for
    // other classes, but whatever pybind11 manages is probably not a performance-critical part, so
    // it should be fine.
    auto cls = py::class_<AlgorithmType, BaseType>(module, std::move(name),
                                                   py::multiple_inheritance());

    cls.doc() = MakeDocString<AlgorithmType>();
    cls.def(py::init<>());
    return cls;
}

template <typename T>
struct MemberPointerClassHelper;

template <typename Class, typename Type>
struct MemberPointerClassHelper<Type Class::*> {
    using type = Class;
};

template <typename MemberPointer>
using MemberPointerClass = typename MemberPointerClassHelper<MemberPointer>::type;
}  // namespace detail

template <typename Base, typename Default, typename... Others>
void BindAlgos(pybind11::module_& primitive_module,
               std::array<char const*, sizeof...(Others) + 1> algo_names) {
    auto algos_module = primitive_module.def_submodule("algorithms");
    auto arr_iter = algo_names.begin();
    auto default_algorithm = detail::RegisterAlgorithm<Default, Base>(algos_module, *arr_iter++);
    (detail::RegisterAlgorithm<Others, Base>(algos_module, *arr_iter++), ...);
    algos_module.attr("Default") = default_algorithm;
}

template <typename Default, typename... Others>
void BindPrimitive(
        pybind11::module_& module, auto result_method, char const* base_name,
        char const* base_result_method_name,
        std::array<char const*, sizeof...(Others) + 1> algo_names,
        // If using the default, make sure your primitive's objects will not become invalid after
        // the algorithm's object is destroyed (for example, if they contain a plain pointer to
        // schema). If they cannot work without the algorithm's object, rework them.
        pybind11::return_value_policy result_rv_policy = pybind11::return_value_policy::copy) {
    namespace py = pybind11;
    using algos::Algorithm;

    using ResultMethodType = std::decay_t<decltype(result_method)>;
    static_assert(std::is_member_pointer_v<ResultMethodType>);
    using Base = detail::MemberPointerClass<ResultMethodType>;
    py::class_<Base, Algorithm>(module, base_name)
            .def(base_result_method_name, result_method, result_rv_policy);
    BindAlgos<Base, Default, Others...>(module, std::move(algo_names));
}

template <typename AlgorithmType>
auto BindNoBaseAndSetDefault(pybind11::module_& module, char const* algo_name) {
    namespace py = pybind11;
    using algos::Algorithm;

    auto algos_module = module.def_submodule("algorithms");
    auto default_algorithm =
            detail::RegisterAlgorithm<AlgorithmType, Algorithm>(algos_module, algo_name);
    algos_module.attr("Default") = default_algorithm;
    return std::pair{default_algorithm, algos_module};
}

template <typename AlgorithmType>
auto BindPrimitiveNoBase(pybind11::module_& module, char const* algo_name) {
    return BindNoBaseAndSetDefault<AlgorithmType>(module, algo_name).first;
}

}  // namespace python_bindings
