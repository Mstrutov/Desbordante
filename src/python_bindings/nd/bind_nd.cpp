#include "nd/bind_nd.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "algorithms/nd/nd.h"
#include "py_util/bind_primitive.h"

namespace py = pybind11;

namespace python_bindings {
void BindNd(py::module_& main_module) {
    using namespace model;

    auto nd_module = main_module.def_submodule("nd");
    py::class_<ND>(nd_module, "ND")
            .def("__str__", &ND::ToLongString)
            .def("to_short_string", &ND::ToShortString)
            .def("to_long_string", &ND::ToLongString)
            .def("get_lhs", &ND::GetLhs)
            .def("get_rhs", &ND::GetRhs)
            .def("get_lhs_indices", &ND::GetLhsIndices)
            .def("get_rhs_indices", &ND::GetRhsIndices)
            .def("get_weight", &ND::GetWeight);
}
}  // namespace python_bindings
