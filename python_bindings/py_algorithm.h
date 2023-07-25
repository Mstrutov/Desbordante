#pragma once

#include <memory>
#include <type_traits>
#include <unordered_map>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "algorithms/algorithm.h"
#include "util/config/tabular_data/input_table_type.h"

namespace python_bindings {

// Currently, only two scenarios are supported. Either
// (scenario 1)
//  1. create the algorithm
//  2. set options using `set_option`
//  3. call `load_data` with no arguments
//  4. set options using `set_option`
//  5. call `execute` with no arguments
// or
// (scenario 2)
//  1. create the algorithm
//  2. call `load_data` with option arguments
//  3. call `execute` with option arguments
// Using the class in any other way from Python is undefined behaviour.
// scenario 2 is intended for use from Python's interactive shell
// scenario 1 is intended for use in scripts where parameters are set one-by-one
class PyAlgorithmBase {
private:
    bool set_option_used_on_stage_ = false;

    void LoadProvidedData(pybind11::kwargs const& kwargs, util::config::InputTable table);

protected:
    std::unique_ptr<algos::Algorithm> algorithm_;

    explicit PyAlgorithmBase(std::unique_ptr<algos::Algorithm> ptr) : algorithm_(std::move(ptr)) {}

    void Configure(pybind11::kwargs const& kwargs, util::config::InputTable table = nullptr);

public:
    // (scenario 1)
    void SetOption(std::string_view option_name, pybind11::handle option_value);

    [[nodiscard]] std::unordered_set<std::string_view> GetNeededOptions() const;

    [[nodiscard]] std::unordered_set<std::string_view> GetPossibleOptions() const;

    [[nodiscard]] std::string_view GetDescription(std::string_view option_name) const;

    [[nodiscard]] pybind11::frozenset GetOptionType(std::string_view option_name) const;

    // For pandas dataframes
    // (scenario 2)
    void LoadData(pybind11::handle dataframe, std::string name, pybind11::kwargs const& kwargs);

    // (scenario 2)
    void LoadData(std::string_view path, char separator, bool has_header,
                  pybind11::kwargs const& kwargs);

    // (scenario 1)
    void LoadData();

    pybind11::int_ Execute(pybind11::kwargs const& kwargs);
};

template <typename AlgorithmType, typename Base>
class PyAlgorithm : public Base {
    static_assert(std::is_base_of_v<PyAlgorithmBase, Base>);

public:
    PyAlgorithm() : Base(std::make_unique<AlgorithmType>()) {}
};

}  // namespace python_bindings
