#pragma once

#include <pybind11/pybind11.h>

namespace python_bindings {
void BindNd(pybind11::module_& main_module);
}  // namespace python_bindings
