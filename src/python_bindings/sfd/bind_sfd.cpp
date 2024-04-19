#include "bind_sfd.h"

#include <algorithm>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "algorithms/fd/sfd/cords.h"
#include "algorithms/fd/sfd/correlation.h"
#include "py_util/bind_primitive.h"

namespace {
namespace py = pybind11;
}  // namespace

namespace python_bindings {

void BindSFD(py::module_& main_module) {
    using namespace algos;
    auto sfd_module = main_module.def_submodule("sfd");
    py::class_<Correlation>(sfd_module, "Correlation")
            .def("__str__", &Correlation::ToString)
            .def("to_string", &Correlation::ToString)
            .def("GetLhsIndex", &Correlation::GetLhsIndex)
            .def("GetRhsIndex", &Correlation::GetRhsInex)
            .def("GetLhsName", &Correlation::GetLhsName)
            .def("GetRhsName", &Correlation::GetRhsName);

    auto sfd_algorithms_module = sfd_module.def_submodule("algorithms");
    auto cls = py::class_<Cords, FDAlgorithm>(sfd_algorithms_module, "SFDAlgorithm")
                       .def(py::init<>())
                       .def("get_correlations", &Cords::GetCorrelations,
                            py::return_value_policy::reference_internal);
    sfd_algorithms_module.attr("Default") = cls;
}
}  // namespace python_bindings
