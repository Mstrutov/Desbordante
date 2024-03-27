#include "bind_sfd.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "algorithms/sfd/sfd_algorithm.h"
#include "py_util/bind_primitive.h"

namespace {
namespace py = pybind11;
}  // namespace

namespace python_bindings {
void BindSFD(py::module_& main_module) {
    using namespace algos;
    auto sfd_module = main_module.def_submodule("sfd");
    BindPrimitiveNoBase<SFDAlgorithm>(sfd_module, "CordsAlgorithm")
            .def("get_sfd", &SFDAlgorithm::GetSFD, py::return_value_policy::reference_internal)
            .def("get_correlations", &SFDAlgorithm::GetCorrelations,
                 py::return_value_policy::reference_internal)
            .def("get_trivial", &SFDAlgorithm::GetTrivialColumns,
                 py::return_value_policy::reference_internal)
            .def("get_soft_keys", &SFDAlgorithm::GetSoftKeys,
                 py::return_value_policy::reference_internal)
            .def("get_unrelated", &SFDAlgorithm::GetUnrelated,
                 py::return_value_policy::reference_internal);
}
}  // namespace python_bindings
