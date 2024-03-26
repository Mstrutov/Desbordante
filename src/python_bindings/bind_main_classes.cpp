#include "bind_main_classes.h"

#include <typeindex>
#include <typeinfo>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
//#include <pybind11/stl/filesystem.h>

#include "algorithms/algo_factory.h"
#include "algorithms/algorithm.h"
#include "config/exceptions.h"
#include "config/names.h"
#include "py_util/get_py_type.h"
#include "py_util/opt_to_py.h"
#include "py_util/py_to_any.h"

namespace {
namespace py = pybind11;
using algos::Algorithm;
auto const kVoidIndex = std::type_index{typeid(void)};

void ConfigureAlgo(Algorithm& algorithm, py::kwargs const& kwargs) {
    using python_bindings::PyToAny;
    algos::ConfigureFromFunction(
            algorithm, [&kwargs, &algorithm](std::string_view option_name) -> boost::any {
                std::type_index type_index = algorithm.GetTypeIndex(option_name);
                assert(type_index != kVoidIndex);
                return kwargs.contains(option_name)
                               ? PyToAny(option_name, type_index, kwargs[option_name.data()])
                               : boost::any{};
            });
}
}  // namespace

namespace python_bindings {
void BindMainClasses(py::module_& main_module) {
    using namespace pybind11::literals;

    main_module.doc() =
            "A high-performance data profiling library oriented towards exploratory data analysis";

    py::register_exception<config::ConfigurationError>(main_module, "ConfigurationError",
                                                       PyExc_ValueError);

#define CERTAIN_SCRIPTS_ONLY                                                       \
    "\nThis option is only expected to be used by Python scripts in which it is\n" \
    "easier to set all options one by one. For normal use, you may set the\n"      \
    "algorithms' options using keyword arguments of the load_data and execute\nmethods."
    py::class_<Algorithm>(main_module, "Algorithm")
            .def(
                    "load_data",
                    [](Algorithm& algo, py::kwargs const& kwargs) {
                        ConfigureAlgo(algo, kwargs);
                        algo.LoadData();
                    },
                    "Load data for execution")
            .def("get_possible_options", &Algorithm::GetPossibleOptions,
                 "Get names of options the algorithm may request.")
            .def("get_description", &Algorithm::GetDescription, "option_name"_a,
                 "Get description of an option.")
            .def("get_needed_options", &Algorithm::GetNeededOptions,
                 "Get names of options the algorithm requires to be set at the "
                 "moment." CERTAIN_SCRIPTS_ONLY)
            .def(
                    "get_option_type",
                    [](Algorithm const& algo, std::string_view option_name) {
                        auto type_index = algo.GetTypeIndex(option_name);
                        if (type_index == kVoidIndex)
                            throw config::ConfigurationError{std::string{"Option named \""} +
                                                             option_name.data() +
                                                             "\" doesn't exist!"};
                        return GetPyType(type_index);
                    },
                    "option_name"_a, "Get info about the option's type.")
            .def(
                    "set_option",
                    [](Algorithm& algorithm, std::string_view option_name,
                       py::handle option_value) {
                        if (option_value.is_none()) {
                            algorithm.SetOption(option_name);
                            return;
                        }
                        algorithm.SetOption(
                                option_name,
                                PyToAny(option_name, algorithm.GetTypeIndex(option_name),
                                        option_value));
                    },
                    "option_name"_a, "option_value"_a = py::none(),
                    "Set option value. Passing None means setting the default "
                    "value." CERTAIN_SCRIPTS_ONLY)
            .def(
                    "get_opts",
                    [](Algorithm& algorithm) {
                        auto opt_value_info = algorithm.GetOptValues();
                        std::unordered_map<std::string_view, pybind11::object> res;
                        for (auto const& [name, value_info] : opt_value_info) {
                            if (name == config::names::kTable) {
                                continue;
                            }
                            res[name] = OptToPy(value_info.type, value_info.value);
                        }
                        return res;
                    },
                    "Get option values represented as the closest Python type")
            .def(
                    "execute",
                    [](Algorithm& algo, py::kwargs const& kwargs) {
                        ConfigureAlgo(algo, kwargs);
                        algo.Execute();
                    },
                    "Process data.");
#undef CERTAIN_SCRIPTS_ONLY
}
}  // namespace python_bindings
