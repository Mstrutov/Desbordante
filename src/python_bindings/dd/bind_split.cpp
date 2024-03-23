#include "bind_split.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "algorithms/dd/dd.h"
#include "algorithms/dd/mining_algorithms.h"
#include "py_util/bind_primitive.h"

namespace {
namespace py = pybind11;
}  // namespace

namespace python_bindings {
void BindSplit(py::module_& main_module) {
    using namespace algos;

    auto dd_module = main_module.def_submodule("dd");
    py::class_<model::DDString>(dd_module, "DD")
            .def("__str__", &model::DDString::ToString)
            .def("__repr__", &model::DDString::ToString);

    BindPrimitiveNoBase<dd::Split>(dd_module, "Split").def("get_dds", &dd::Split::GetDDStringList);
}
}  // namespace python_bindings
