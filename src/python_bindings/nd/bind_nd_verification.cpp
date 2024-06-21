#include "bind_nd_verification.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "algorithms/nd/nd_verifier/nd_verifier.h"
#include "algorithms/nd/nd_verifier/util/highlight.h"
#include "py_util/bind_primitive.h"

namespace {
namespace py = pybind11;
}  // namespace

namespace python_bindings {
void BindNdVerification(pybind11::module_& main_module) {
    using namespace algos;
    using namespace algos::nd_verifier;
    using namespace algos::nd_verifier::util;

    auto nd_verification_module = main_module.def_submodule("nd_verification");
    py::class_<Highlight>(nd_verification_module, "Highlight")
            .def("__str__", &Highlight::ToValuesString)
            .def("to_short_string", &Highlight::ToIndicesString)
            .def("to_long_string", &Highlight::ToValuesString)
            .def("get_occurences_indices", &Highlight::GetOccurencesIndices)
            .def_property_readonly("occurences_number", &Highlight::GetOccurencesNumber)
            .def_property_readonly("lhs_value", &Highlight::GetLhsValue)
            .def("get_rhs_values", &Highlight::GetRhsValues)
            .def_property_readonly("distinct_rhs_values_number",
                                   &Highlight::GetDistinctRhsValuesNumber)
            .def("get_most_frequent_rhs_value_indices", &Highlight::GetMostFrequentRhsValueIndices)
            .def("get_most_frequent_rhs_values", &Highlight::GetMostFrequentRhsValues);

    BindPrimitiveNoBase<NDVerifier>(nd_verification_module, "NDVerifier")
            .def_property_readonly("nd_holds", &NDVerifier::NDHolds)
            .def_property_readonly("highlights", &NDVerifier::GetHighlights)
            .def_property_readonly("global_min_weight", &NDVerifier::GetGlobalMinWeight)
            .def_property_readonly("real_weight", &NDVerifier::GetRealWeight)
            .def("get_lhs_frequencies", &NDVerifier::GetLhsFrequencies)
            .def("get_rhs_frequencies", &NDVerifier::GetRhsFrequencies);
}
}  // namespace python_bindings
