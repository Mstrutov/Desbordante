#include "bind_gfd_verification.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "algorithms/gfd/mining_algorithms.h"
#include "py_util/bind_primitive.h"

namespace {
namespace py = pybind11;
}  // namespace

namespace python_bindings {
void BindGfdVerification(pybind11::module_& main_module) {
    using namespace algos;

    auto gfd_module = main_module.def_submodule("gfd_verification");

    BindPrimitive<GfdValidation, EGfdValidation, NaiveGfdValidation>(
            gfd_module, &GfdHandler::StringGfdList, "GfdAlgorithm", "get_gfds",
            {"GfdValid", "EGfdValid", "NaiveGfdValid"});
}
}  // namespace python_bindings
