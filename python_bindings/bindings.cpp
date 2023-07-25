// This file exposes Desbordante algorithm classes to Python.
// This Python library heavily uses trampoline classes, which share a lot of
// methods. In theory, it all can be done with fewer classes, but pybind11
// requires each Python class to have a corresponding C++ class, which is why
// all algorithms use the PyAlgorithm template.
// The classes which algorithms inherit from ("bases") like `FDAlgorithm` are
// also exposed to Python to allow for `isinstance` and `issubclass` checks.
// This file defines several macros that help reduce code repetition. These are:
// DEFINE_ALGORITHM: register a Python class `type` inherited from `base`. All
// algorithm trampoline classes are named using the scheme `Py[name]`, and
// [name] is also used as the corresponding Python class's name.
// DEFINE_ALGORITHM_BASE: register an algorithm base class. All algorithm base
// trampoline class are named using the scheme `Py[name]Base`, where [name] is
// the name of the corresponding Python class. An algorithm base is expected to
// implement methods that get execution results of the algorithms. The main base
// `PyAlgorithmBase` is defined separately before any usages of this macro, and
// it exposes methods common to all algorithms.
// DEFINE_{FD,AR}_ALGORITHM: register a class inherited from the corresponding
// base class. These base classes implement methods for getting results. Do note
// that some algorithms are inherited directly from `Algorithm`, as they don't
// belong to the "normal" algorithm types like FDAlgorithm. One example of such
// a class is `FDVerifier`.
// If you want to add your own algorithm, you need to either use the one of the
// already existing bases if applicable or inherit from
// PyAlgorithm<YourAlgorithmType, PyAlgorithmBase> and implement your own result
// getter methods. To be able to use the aforementioned macros, you have to name
// your class using the scheme Py[PythonClassName]. If you want to return custom
// result classes, you need to register them using
// `py::class<YourResultClass>(module, "PythonClassName")`.
// Consult pybind11's documentation for details.

#include <filesystem>

#include <easylogging++.h>
#include <pybind11/pybind11.h>

#include "algorithms/algorithms.h"
#include "model/ar.h"
#include "py_ac_algorithm.h"
#include "py_ar_algorithm.h"
#include "py_data_stats.h"
#include "py_fd_algorithm.h"
#include "py_fd_verifier.h"
#include "py_metric_verifier.h"
#include "py_ucc_algorithm.h"

INITIALIZE_EASYLOGGINGPP

#define DEFINE_ALGORITHM(type, base) \
    py::class_<Py##type, Py##base##Base>(module, #type).def(py::init<>())
#define DEFINE_ALGORITHM_BASE(base) py::class_<Py##base##Base, PyAlgorithmBase>(module, #base)
#define DEFINE_FD_ALGORITHM(type) DEFINE_ALGORITHM(type, FdAlgorithm)
#define DEFINE_AR_ALGORITHM(type) DEFINE_ALGORITHM(type, ArAlgorithm)
#define DEFINE_UCC_ALGORITHM(type) DEFINE_ALGORITHM(type, UCCAlgorithm)

namespace python_bindings {

namespace py = pybind11;
// AR mining algorithms
using PyApriori = PyArAlgorithm<algos::Apriori>;
// FD mining algorithms
using PyTane = PyFDAlgorithm<algos::Tane>;
using PyPyro = PyFDAlgorithm<algos::Pyro>;
using PyFUN = PyFDAlgorithm<algos::FUN>;
using PyFdMine = PyFDAlgorithm<algos::Fd_mine>;
using PyFastFDs = PyFDAlgorithm<algos::FastFDs>;
using PyHyFD = PyFDAlgorithm<algos::hyfd::HyFD>;
using PyFDep = PyFDAlgorithm<algos::FDep>;
using PyDFD = PyFDAlgorithm<algos::DFD>;
using PyDepminer = PyFDAlgorithm<algos::Depminer>;
using PyAid = PyFDAlgorithm<algos::Aid>;
// UCC mining algorithms
using PyHyUCC = PyUCCAlgorithm<algos::HyUCC>;

using FDHighlight = algos::fd_verifier::Highlight;
using algos::ACException;
using algos::RangesCollection;
using model::ARStrings;

PYBIND11_MODULE(desbordante, module) {
    using namespace pybind11::literals;

    if (std::filesystem::exists("logging.conf")) {
        el::Loggers::configureFromGlobal("logging.conf");
    } else {
        el::Configurations conf;
        conf.set(el::Level::Global, el::ConfigurationType::Enabled, "false");
        el::Loggers::reconfigureAllLoggers(conf);
    }

    module.doc() = "A data profiling library";

    py::class_<ARStrings>(module, "AssociativeRule")
            .def("__str__", &ARStrings::ToString)
            .def_readonly("left", &ARStrings::left)
            .def_readonly("right", &ARStrings::right)
            .def_readonly("confidence", &ARStrings::confidence);

    py::class_<PyFD>(module, "FD")
            .def("__str__", &PyFD::ToString)
            .def("__repr__", &PyFD::ToString)
            .def_property_readonly("lhs_indices", &PyFD::GetLhs)
            .def_property_readonly("rhs_index", &PyFD::GetRhs);

    py::class_<PyUCC>(module, "UCC")
            .def("__str__", &PyUCC::ToString)
            .def("__repr__", &PyUCC::ToString)
            .def_property_readonly("indices", &PyUCC::GetUCC);

    py::class_<FDHighlight>(module, "FDHighlight")
            .def_property_readonly("cluster", &FDHighlight::GetCluster)
            .def_property_readonly("num_distinct_rhs_values", &FDHighlight::GetNumDistinctRhsValues)
            .def_property_readonly("most_frequent_rhs_value_proportion",
                                   &FDHighlight::GetMostFrequentRhsValueProportion);

    py::class_<PyRangesCollection>(module, "ACRanges")
            .def_property_readonly("column_indices", &PyRangesCollection::GetColumnIndices)
            .def_property_readonly("ranges", &PyRangesCollection::GetRanges);

    py::class_<ACException>(module, "ACException")
            .def_readonly("row_index", &ACException::row_i)
            .def_readonly("column_pairs", &ACException::column_pairs);

    py::class_<PyAlgorithmBase>(module, "Algorithm")
            .def("load_data",
                 py::overload_cast<std::string_view, char, bool, py::kwargs const&>(
                         &PyAlgorithmBase::LoadData),
                 "path"_a, "separator"_a, "has_header"_a, "Load data from CSV")
            .def("load_data",
                 py::overload_cast<pybind11::handle, std::string, py::kwargs const&>(
                         &PyAlgorithmBase::LoadData),
                 "df"_a, "name"_a = "Pandas dataframe", "Load data from pandas dataframe")
            .def("load_data", py::overload_cast<>(&PyAlgorithmBase::LoadData),
                 "Load data after all options have been set by SetOption calls")
            .def("get_needed_options", &PyAlgorithmBase::GetNeededOptions,
                 "Get names of options the algorithm needs")
            .def("get_possible_options", &PyAlgorithmBase::GetPossibleOptions,
                 "Get names of options the algorithm may request")
            .def("get_description", &PyAlgorithmBase::GetDescription,
                 "Get description of an option of this algorithm")
            .def("set_option", &PyAlgorithmBase::SetOption, "option_name"_a,
                 "option_value"_a = pybind11::none(), "Set option value")
            .def("get_option_type", &PyAlgorithmBase::GetOptionType, "option_name"_a,
                 "Get info about the option's type")
            .def("execute", &PyAlgorithmBase::Execute, "Process data");

    DEFINE_ALGORITHM_BASE(ArAlgorithm).def("get_ars", &PyArAlgorithmBase::GetARs);
    DEFINE_ALGORITHM_BASE(FdAlgorithm).def("get_fds", &PyFdAlgorithmBase::GetFDs);
    DEFINE_ALGORITHM_BASE(UCCAlgorithm).def("get_uccs", &PyUCCAlgorithmBase::GetUCCs);

    DEFINE_ALGORITHM(FDVerifier, Algorithm)
            .def("fd_holds", &PyFDVerifier::FDHolds)
            .def("get_error", &PyFDVerifier::GetError)
            .def("get_num_error_clusters", &PyFDVerifier::GetNumErrorClusters)
            .def("get_num_error_rows", &PyFDVerifier::GetNumErrorRows)
            .def("get_highlights", &PyFDVerifier::GetHighlights);

    DEFINE_ALGORITHM(DataStats, Algorithm).def("get_result_string", &PyDataStats::GetResultString);

    DEFINE_ALGORITHM(MetricVerifier, Algorithm).def("mfd_holds", &PyMetricVerifier::MfdHolds);

    DEFINE_AR_ALGORITHM(Apriori);

    DEFINE_ALGORITHM(ACAlgorithm, Algorithm)
            .def("get_ac_ranges", &PyACAlgorithm::GetACRanges)
            .def("get_ac_exceptions", &PyACAlgorithm::GetACExceptions);

    DEFINE_FD_ALGORITHM(Aid);
    DEFINE_FD_ALGORITHM(Depminer);
    DEFINE_FD_ALGORITHM(DFD);
    DEFINE_FD_ALGORITHM(FastFDs);
    DEFINE_FD_ALGORITHM(FDep);
    DEFINE_FD_ALGORITHM(FdMine);
    DEFINE_FD_ALGORITHM(FUN);
    DEFINE_FD_ALGORITHM(HyFD);
    DEFINE_FD_ALGORITHM(Pyro);
    DEFINE_FD_ALGORITHM(Tane);

    DEFINE_UCC_ALGORITHM(HyUCC);
}
#undef DEFINE_FD_ALGORITHM
#undef DEFINE_AR_ALGORITHM
#undef DEFINE_UCC_ALGORITHM
#undef DEFINE_ALGORITHM_BASE
#undef DEFINE_ALGORITHM

}  // namespace python_bindings
