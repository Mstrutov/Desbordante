#include <functional>
#include <unordered_map>

#include <boost/any.hpp>
#include <easylogging++.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "algorithms/algebraic_constraints/bin_operation_enum.h"
#include "algorithms/metric/enums.h"
#include "association_rules/ar_algorithm_enums.h"
#include "config/exceptions.h"
#include "config/tabular_data/input_table_type.h"
#include "parser/csv_parser/csv_parser.h"
#include "py_util/create_dataframe_reader.h"
#include "util/enum_to_available_values.h"

namespace {

namespace py = pybind11;
using ConvFunc = std::function<boost::any(std::string_view, py::handle)>;

template <typename T>
T CastAndReplaceCastError(std::string_view option_name, py::handle value) {
    try {
        return py::cast<T>(value);
    } catch (py::cast_error& e) {
        throw config::ConfigurationError(
                std::string("Unable to cast Python object to C++ type for option \"") +
                option_name.data() + '"');
    }
}

config::InputTable CreateCsvParser(std::string_view option_name, py::tuple const& arguments) {
    if (py::len(arguments) != 3) {
        throw config::ConfigurationError("Cannot create a CSV parser from passed tuple.");
    }

    auto csv_parser = std::make_shared<CSVParser>(
            CastAndReplaceCastError<std::string>(option_name, arguments[0]),
            CastAndReplaceCastError<char>(option_name, arguments[1]),
            CastAndReplaceCastError<bool>(option_name, arguments[2]));
    csv_parser->ValidateSeparator();
    return csv_parser;
}

template <typename Type>
std::pair<std::type_index, ConvFunc> const NormalConvPair{
        std::type_index(typeid(Type)), [](std::string_view option_name, py::handle value) {
            return CastAndReplaceCastError<Type>(option_name, value);
        }};

template <typename EnumType>
std::pair<std::type_index, ConvFunc> const EnumConvPair{
        std::type_index(typeid(EnumType)), [](std::string_view option_name, py::handle value) {
            auto string = CastAndReplaceCastError<std::string>(option_name, value);
            better_enums::optional<EnumType> enum_holder =
                    EnumType::_from_string_nocase_nothrow(string.data());
            if (enum_holder) return *enum_holder;

            std::stringstream error_message;
            error_message << "Incorrect value for option \"" << option_name
                          << "\". Possible values: " << util::EnumToAvailableValues<EnumType>();
            throw config::ConfigurationError(error_message.str());
        }};

template <typename EnumType>
std::pair<std::type_index, ConvFunc> const CharEnumConvPair{
        std::type_index(typeid(EnumType)), [](std::string_view option_name, py::handle value) {
            using EnumValueType = typename EnumType::_integral;
            // May be applicable to other types.
            static_assert(std::is_same_v<EnumValueType, char>);
            auto char_value = CastAndReplaceCastError<char>(option_name, value);
            better_enums::optional<EnumType> enum_holder =
                    EnumType::_from_integral_nothrow(char_value);
            if (enum_holder) return *enum_holder;

            std::stringstream error_message;
            error_message << "Incorrect value for option \"" << option_name
                          << "\". Possible values: ";

            error_message << '[';
            for (auto const& val : EnumType::_values()) {
                error_message << val._to_integral() << '|';
            }
            error_message.seekp(-1, std::stringstream::cur);
            error_message << ']';

            throw config::ConfigurationError(error_message.str());
        }};

boost::any InputTableToAny(std::string_view option_name, py::handle obj) {
    if (py::isinstance<py::tuple>(obj)) {
        return CreateCsvParser(option_name, py::cast<py::tuple>(obj));
    }
    return python_bindings::CreateDataFrameReader(obj);
}

std::unordered_map<std::type_index, ConvFunc> const converters{
        NormalConvPair<bool>,
        NormalConvPair<double>,
        NormalConvPair<unsigned int>,
        NormalConvPair<long double>,
        NormalConvPair<std::vector<unsigned int>>,
        NormalConvPair<ushort>,
        NormalConvPair<int>,
        NormalConvPair<size_t>,
        EnumConvPair<algos::metric::Metric>,
        EnumConvPair<algos::metric::MetricAlgo>,
        EnumConvPair<algos::InputFormat>,
        CharEnumConvPair<algos::Binop>,
        {typeid(config::InputTable), InputTableToAny},
};

}  // namespace

namespace python_bindings {

boost::any PyToAny(std::string_view option_name, std::type_index index, py::handle obj) {
    return converters.at(index)(option_name, obj);
}

}  // namespace python_bindings
