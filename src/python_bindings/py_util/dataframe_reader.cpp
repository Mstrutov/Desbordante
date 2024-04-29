#include "dataframe_reader.h"

#include <string>
#include <utility>
#include <vector>

#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

#include "model/types/builtin.h"

namespace {
namespace py = pybind11;
static std::vector<std::string> GetColumnNames(py::handle dataframe) {
    std::vector<std::string> names;
    py::list name_lst = dataframe.attr("columns").attr("to_list")();
    for (py::handle element : name_lst) {
        names.emplace_back(py::str(element));
    }
    return names;
}
}  // namespace

namespace python_bindings {

DataframeReaderBase::DataframeReaderBase(py::handle dataframe, std::string name)
    : dataframe_(py::reinterpret_borrow<py::object>(dataframe)),
      df_iter_(dataframe_.attr("itertuples")(false, py::none{})),
      name_(std::move(name)),
      column_names_(::GetColumnNames(dataframe_)) {}

void DataframeReaderBase::Reset() {
    df_iter_ = dataframe_.attr("itertuples")(false, py::none{});
}

std::string DataframeReaderBase::GetRelationName() const {
    return name_;
}

std::string DataframeReaderBase::GetColumnName(size_t index) const {
    return column_names_.at(index);
}

size_t DataframeReaderBase::GetNumberOfColumns() const {
    return column_names_.size();
}

bool DataframeReaderBase::HasNextRow() const {
    return df_iter_ != py::iterator::sentinel();
}

std::vector<std::string> StringDataframeReader::GetNextRow() {
    return py::cast<std::vector<std::string>>(*df_iter_++);
}

std::vector<std::string> ArbitraryDataframeReader::GetNextRow() {
    std::vector<std::string> strings{};
    auto tuple_row = py::reinterpret_borrow<py::object>(*df_iter_);
    ++df_iter_;
    for (py::handle el : tuple_row) {
        // When reading from .csv files, pandas treats some values as nulls,
        // and empty values are among those values, which may cause some
        // confusion here, since in Desbordante only the literal "NULL" string
        // is interpreted as the null value.
        strings.emplace_back(is_null_(el) ? model::Null::kValue : py::str(el));
    }
    return strings;
}

}  // namespace python_bindings
