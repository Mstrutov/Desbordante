#include <algorithm>
#include <filesystem>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/program_options.hpp>
#include <easylogging++.h>

#include "algorithms/algo_factory.h"
#include "algorithms/ar_algorithm_enums.h"
#include "algorithms/metric/enums.h"
#include "algorithms/options/descriptions.h"
#include "algorithms/options/names.h"
#include "algorithms/create_primitive.h"
#include "algorithms/cfd/enums.h"

namespace po = boost::program_options;
namespace onam = algos::config::names;
namespace desc = algos::config::descriptions;

using algos::EnumToAvailableValues;

INITIALIZE_EASYLOGGINGPP

namespace algos::metric {

void validate(boost::any& v, const std::vector<std::string>& values, Metric*, int) {
    const std::string& s = po::validators::get_single_string(values);
    try {
        v = boost::any(Metric::_from_string_nocase(s.c_str()));
    } catch (std::runtime_error&e) {
        throw po::validation_error(po::validation_error::invalid_option_value);
    }
}

void validate(boost::any& v, const std::vector<std::string>& values, MetricAlgo*, int) {
    const std::string& s = po::validators::get_single_string(values);
    try {
        v = boost::any(MetricAlgo::_from_string_nocase(s.c_str()));
    } catch (std::runtime_error &e) {
        throw po::validation_error(po::validation_error::invalid_option_value);
    }
}

}  // namespace algos::metric

namespace algos {

void validate(boost::any& v, const std::vector<std::string>& values, PrimitiveType*, int) {
    const std::string& s = po::validators::get_single_string(values);
    try {
        v = boost::any(PrimitiveType::_from_string_nocase(s.c_str()));
    } catch (std::runtime_error &e) {
        throw po::validation_error(po::validation_error::invalid_option_value);
    }
}

void validate(boost::any& v, const std::vector<std::string>& values, InputFormat*, int) {
    const std::string& s = po::validators::get_single_string(values);
    try {
        v = boost::any(InputFormat::_from_string_nocase(s.c_str()));
    } catch (std::runtime_error &e) {
        throw po::validation_error(po::validation_error::invalid_option_value);
    }
}

void validate(boost::any& v, const std::vector<std::string>& values, CfdAlgo*, int) {
    const std::string& s = po::validators::get_single_string(values);
    try {
        v = boost::any(CfdAlgo::_from_string_nocase(s.c_str()));
    } catch (std::runtime_error& e) {
        throw po::validation_error(po::validation_error::invalid_option_value);
    }
}
}  // namespace algos

int main(int argc, char const* argv[]) {
    std::string primitive;
    /*Options for algebraic constraints algorithm*/
    char bin_operation = '+';
    double fuzziness = 0.15;
    double p_fuzz = 0.9;
    double weight = 0.05;
    size_t bumps_limit = 5;
    size_t iterations_limit = 10;
    std::string pairing_rule = "trivial";

    std::string const separator_arg = std::string{onam::kSeparator} + ",s";
    std::string const prim_desc = "algorithm to use for data profiling\n" +
                                  EnumToAvailableValues<algos::PrimitiveType>() + " + [ac]";
    constexpr auto help_opt = "help";
    constexpr auto prim_opt = "primitive";
    std::string const separator_opt = std::string(onam::kSeparator) + ",s";

    po::options_description info_options("Desbordante information options");
    info_options.add_options()
        (help_opt, "print the help message and exit")
        // --version, if needed, goes here too
        ;

    po::options_description general_options("General options");
    general_options.add_options()
            (prim_opt, po::value<std::string>(&primitive)->required(), prim_desc.c_str())
            (onam::kData, po::value<std::string>()->required(), desc::kDData)
            (separator_opt.c_str(), po::value<char>()->default_value(','), desc::kDSeparator)
            (onam::kHasHeader, po::value<bool>()->default_value(true), desc::kDHasHeader)
            (onam::kEqualNulls, po::value<bool>(), desc::kDEqualNulls)
            (onam::kThreads, po::value<ushort>(), desc::kDThreads)
            ;

    po::options_description fd_options("FD options");
    fd_options.add_options()
            (onam::kError, po::value<double>(), desc::kDError)
            (onam::kMaximumLhs, po::value<unsigned int>(), desc::kDMaximumLhs)
            (onam::kSeed, po::value<int>(), desc::kDSeed)
            ;

    po::options_description typo_options("Typo mining options");
    typo_options.add_options()
            (onam::kRatio, po::value<double>(), desc::kDRatio)
            (onam::kRadius, po::value<double>(), desc::kDRadius)
            (onam::kApproximateAlgorithm, po::value<algos::PrimitiveType>(),
             desc::kDApproximateAlgorithm)
            (onam::kPreciseAlgorithm, po::value<algos::PrimitiveType>(), desc::kDPreciseAlgorithm)
            ;

    po::options_description ar_options("AR options");
    ar_options.add_options()
            (onam::kMinimumSupport, po::value<double>(), desc::kDMinimumSupport)
            (onam::kMinimumConfidence, po::value<double>(), desc::kDMinimumConfidence)
            (onam::kInputFormat, po::value<algos::InputFormat>(), desc::kDInputFormat)
            ;

    po::options_description cfd_options("CFD options");
    cfd_options.add_options()
            (onam::kCfdMinimumSupport, po::value<unsigned>(),desc::kDCfdMinimumSupport)
            (onam::kCfdMinimumConfidence, po::value<double>(), desc::kDCfdMinimumConfidence)
            (onam::kCfdMaximumLhs, po::value<unsigned int>()->default_value(0), desc::kDCfdMaximumLhs)
            (onam::kCfdColumnsNumber, po::value<unsigned int>()->default_value(0), desc::kDCfdColumnsNumber)
            (onam::kCfdTuplesNumber, po::value<unsigned int>()->default_value(0), desc::kDCfdTuplesNumber)
            (onam::kCfdAlgo, po::value<algos::CfdAlgo>(), desc::kDCfdAlgo)
        ;

    po::options_description ar_singular_options("AR \"singular\" input format options");
    ar_singular_options.add_options()
            (onam::kTIdColumnIndex, po::value<unsigned>(), desc::kDTIdColumnIndex)
            (onam::kItemColumnIndex, po::value<unsigned>(), desc::kDItemColumnIndex)
            ;

    po::options_description ar_tabular_options("AR \"tabular\" input format options");
    ar_tabular_options.add_options()
            (onam::kFirstColumnTId, po::bool_switch(), desc::kDFirstColumnTId)
            ;

    ar_options.add(ar_singular_options).add(ar_tabular_options);

    po::options_description mfd_options("MFD options");
    mfd_options.add_options()
            (onam::kMetric, po::value<algos::metric::Metric>(), desc::kDMetric)
            (onam::kMetricAlgorithm, po::value<algos::metric::MetricAlgo>(), desc::kDMetricAlgorithm)
            (onam::kLhsIndices, po::value<std::vector<unsigned int>>()->multitoken(),
             desc::kDLhsIndices)
            (onam::kRhsIndices, po::value<std::vector<unsigned int>>()->multitoken(),
             desc::kDRhsIndices)
            (onam::kParameter, po::value<long double>(), desc::kDParameter)
            (onam::kDistFromNullIsInfinity, po::bool_switch(), desc::kDDistFromNullIsInfinity)
            ;

    po::options_description cosine_options("Cosine metric options");
    cosine_options.add_options()
            (onam::kQGramLength, po::value<unsigned int>(), desc::kDQGramLength)
            ;

    mfd_options.add(cosine_options);

    po::options_description ac_options("AC options");
    ac_options.add_options()
        (onam::kBinaryOperation, po::value<char>(&bin_operation)->default_value(bin_operation),
         "one of availible operations: /, *, +, - ")
        (onam::kFuzziness, po::value<double>(&fuzziness)->default_value(0.15),
         "fraction of exceptional records")
        (onam::kFuzzinessProbability, po::value<double>(&p_fuzz)->default_value(p_fuzz),
         "probability, the fraction of exceptional records that lie outside the "
         "bump intervals is at most Fuzziness")
        (onam::kWeight, po::value<double>(&weight)->default_value(weight),
         "value between 0 and 1. Closer to 0 - many short intervals. "
         "Closer to 1 - small number of long intervals")
        (onam::kBumpsLimit, po::value<size_t>(&bumps_limit)->default_value(bumps_limit),
         "max considered intervals amount. Pass 0 to remove limit")
        (onam::kIterationsLimit, po::value<size_t>(&iterations_limit)->default_value(iterations_limit),
         "limit for iterations of sampling")
        (onam::kPairingRule, po::value<std::string>(&pairing_rule)->default_value(pairing_rule),
         "one of available pairing rules: trivial")
        ;

    po::options_description all_options("Allowed options");
    all_options.add(info_options).add(general_options).add(fd_options)
        .add(mfd_options).add(ar_options).add(ac_options).add(typo_options).add(cfd_options);

    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, all_options), vm);
    } catch (po::error &e) {
        std::cout << e.what() << std::endl;
        return 1;
    }
    if (vm.count(help_opt))
    {
        std::cout << all_options << std::endl;
        return 0;
    }
    try {
        po::notify(vm);
    } catch (po::error &e) {
        std::cout << e.what() << std::endl;
        return 1;
    }

    el::Loggers::configureFromGlobal("logging.conf");

    std::unique_ptr<algos::Primitive> algorithm_instance;
    try {
        algorithm_instance = algos::CreatePrimitive(primitive, vm);
    } catch (std::exception& e) {
        std::cout << e.what() << std::endl;
        return 1;
    }
    try {
        unsigned long long elapsed_time = algorithm_instance->Execute();
        std::cout << "> ELAPSED TIME: " << elapsed_time << std::endl;
    } catch (std::runtime_error& e) {
        std::cout << e.what() << std::endl;
        return 1;
    }

    return 0;
}
