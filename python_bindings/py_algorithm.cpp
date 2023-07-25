#include "py_algorithm.h"

#include <string>
#include <utility>

#include <boost/any.hpp>

#include "algorithms/algo_factory.h"
#include "create_dataframe_reader.h"
#include "dataframe_reader.h"
#include "get_py_type.h"
#include "parser/csv_parser.h"
#include "py_to_any.h"
#include "util/config/names.h"
#include "util/config/tabular_data/input_table_type.h"

namespace python_bindings {

using util::config::InputTable;
using util::config::names::kTable;

static const auto void_index = std::type_index{typeid(void)};

namespace py = pybind11;

void PyAlgorithmBase::Configure(py::kwargs const& kwargs, InputTable table) {
    algorithm_->ResetConfiguration();
    set_option_used_on_stage_ = false;
    algos::ConfigureFromFunction(
            *algorithm_,
            [this, &kwargs, table = std::move(table)](std::string_view option_name) -> boost::any {
                if (option_name == kTable) {
                    // If this is reached, it is assumed the algorithm is currently at the data
                    // loading stage, and thus a `LoadData` overload of PyAlgorithm has been called,
                    // which means the table pointer is not nullptr.
                    assert(table != nullptr);
                    return table;
                }
                std::type_index type_index = algorithm_->GetTypeIndex(option_name);
                assert(type_index != void_index);
                return kwargs.contains(option_name)
                               ? PyToAny(type_index, kwargs[py::str{option_name}])
                               : boost::any{};
            });
}

void PyAlgorithmBase::SetOption(std::string_view option_name, py::handle option_value) {
    set_option_used_on_stage_ = true;
    if (option_value.is_none()) {
        algorithm_->SetOption(option_name);
        return;
    }
    algorithm_->SetOption(option_name,
                          PyToAny(algorithm_->GetTypeIndex(option_name), option_value));
}

std::unordered_set<std::string_view> PyAlgorithmBase::GetNeededOptions() const {
    return algorithm_->GetNeededOptions();
}

std::unordered_set<std::string_view> PyAlgorithmBase::GetPossibleOptions() const {
    return algorithm_->GetPossibleOptions();
}

std::string_view PyAlgorithmBase::GetDescription(std::string_view option_name) const {
    return algorithm_->GetDescription(option_name);
}

py::frozenset PyAlgorithmBase::GetOptionType(std::string_view option_name) const {
    auto type_index = algorithm_->GetTypeIndex(option_name);
    if (type_index == void_index)
        throw std::invalid_argument{std::string{"Option named \""} + option_name.data() +
                                    "\" doesn't exist!"};
    return GetPyType(type_index);
}

void PyAlgorithmBase::LoadProvidedData(pybind11::kwargs const& kwargs, InputTable table) {
    Configure(kwargs, std::move(table));
    algorithm_->LoadData();
}

void PyAlgorithmBase::LoadData(std::string_view path, char separator, bool has_header,
                               py::kwargs const& kwargs) {
    LoadProvidedData(kwargs, std::make_shared<CSVParser>(path, separator, has_header));
}

void PyAlgorithmBase::LoadData(py::handle dataframe, std::string name, py::kwargs const& kwargs) {
    LoadProvidedData(kwargs, CreateDataFrameReader(dataframe, std::move(name)));
}

void PyAlgorithmBase::LoadData() {
    algorithm_->LoadData();
    set_option_used_on_stage_ = false;
}

py::int_ PyAlgorithmBase::Execute(py::kwargs const& kwargs) {
    if (!(set_option_used_on_stage_ && kwargs.empty())) {
        Configure(kwargs);
    }
    py::int_ time = algorithm_->Execute();
    set_option_used_on_stage_ = false;
    return time;
}

}  // namespace python_bindings
