#include "borgelt.h"

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/process.hpp>
#include <easylogging++.h>

#include "descriptions.h"
#include "names.h"

namespace {

    struct Config {
        using MinSupType = double;
        double minconf_;
        bool first_column_tid_;
        MinSupType minsup_;
        std::string sep_;
        std::filesystem::path path_;
        std::string algo_;
        std::filesystem::path output_file;
    };

    std::string PrepareCall(Config const& config) {
        std::filesystem::path path = config.path_;
        std::string sep = {config.sep_};
        std::string minsup = std::to_string(static_cast<int>(config.minsup_ * 100));
        std::string minconf = std::to_string(static_cast<int>(config.minconf_ * 100));

        std::filesystem::path const kOutputTmp = config.output_file;
        std::filesystem::path const kExternalExecutable = config.algo_;
        if (!std::filesystem::exists(kExternalExecutable)) {
            throw std::runtime_error("No executable named " + kExternalExecutable.string() +
                                     " found in " + std::filesystem::current_path().string());
        }
        boost::format call = boost::format{"./" + config.algo_ + " -tr -f%1% -s%2% -c%3% %4% %5%"} %
                             sep % minsup % minconf % path.generic_string() % kOutputTmp.string();
        return call.str();
    }

    std::string ExtractOutput(boost::process::ipstream& pipe_stream) {
        std::stringstream output;

        std::string line;
        while (pipe_stream && std::getline(pipe_stream, line) && !line.empty()) {
            output << line << '\n';
        }

        return output.str();
    }

    unsigned long long ExecutePrepared(std::string const& call) {
        auto const start = std::chrono::system_clock::now();

        boost::process::ipstream pipe_stream;
        std::error_code error_code;
        boost::process::system(call, boost::process::std_err > pipe_stream, error_code);

        auto const output = ExtractOutput(pipe_stream);
        if (error_code) {
            throw std::runtime_error(output);
        }
        LOG(INFO) << output;

        return (std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::system_clock::now() - start))
                .count();
    }

    unsigned long long ExecuteExternal(Config const& config) {
        auto const call = PrepareCall(config);
        LOG(INFO) << "Calling `" << call << "`";

        return ExecutePrepared(call);
    }

    using ItemNamePosMap = std::unordered_map<std::string, size_t>;

    auto ToItemNamePosMap(std::vector<std::string> const& itemNamesVector) {
        ItemNamePosMap itemNameId;
        for (size_t i = 0; i < itemNamesVector.size(); i++) {
            itemNameId[itemNamesVector[i]] = i;
        }
        return itemNameId;
    }

    /**
     * Generates compressed, i.e. with item ids instead of item names, list of association rules by
     * parsing external tool output file, which is expected to look like a collection of lines
     * <item> <- (<item> )+\(<support>, <confidence>\)
     *
     * @param kOutputTmp path to external tool output file
     * @param item_name_id map of item names to item ids
     * @return compressed list of association rules
     */
    auto Parse(std::filesystem::path const& kOutputTmp, ItemNamePosMap const& item_name_id) {
        std::list<model::ArIDs> ar_collection;

        std::ifstream ifs(kOutputTmp);
        std::string line;
        while (std::getline(ifs, line)) {
            size_t sep_pos = line.find(" <- ");

            std::vector<unsigned> rhs;
            {
                std::string full_rhs = line.substr(0, sep_pos);
                std::vector<std::string> splitRhs;
                boost::split(
                        splitRhs, full_rhs, [](auto const& ch) { return ch == ' '; },
                        boost::token_compress_on);
                std::transform(splitRhs.begin(), splitRhs.end(), std::back_inserter(rhs),
                               [&item_name_id](auto const& el) { return item_name_id.at(el); });
            }

            size_t par_pos = line.find('(');
            std::vector<unsigned> lhs;
            {
                unsigned const kSidesSeparatorLength = 4;
                std::string full_lhs = line.substr(sep_pos + kSidesSeparatorLength,
                                                   par_pos - 1 - sep_pos - kSidesSeparatorLength);
                std::vector<std::string> splitLhs;
                boost::split(
                        splitLhs, full_lhs, [](auto const& ch) { return ch == ' '; },
                        boost::token_compress_on);
                auto const end = std::remove_if(splitLhs.begin(), splitLhs.end(),
                                                [](auto const& el) { return el.empty(); });
                std::transform(splitLhs.begin(), end, std::back_inserter(lhs),
                               [&item_name_id](auto const& el) { return item_name_id.at(el); });
            }

            double confi;
            {
                std::string pars = line.substr(par_pos + 1, line.size() - 1 - par_pos - 1);
                std::vector<std::string> splitPars;
                boost::split(
                        splitPars, pars, [](auto const& ch) { return ch == ','; },
                        boost::token_compress_on);
                confi = std::stod(splitPars.at(1).substr(1));
            }

            ar_collection.emplace_back(std::move(lhs), std::move(rhs), confi);
        }

        return ar_collection;
    }
}

namespace algos {

decltype(Borgelt::AlgoOpt) Borgelt::AlgoOpt{
        {config::names::kARAlgorithm, config::descriptions::kDARAlgorithm}};

void Borgelt::ConvertAndFillResults() {
    auto const& itemNamesVector = GetItemNamesVector();

    auto const item_name_id = ToItemNamePosMap(itemNamesVector);

    ar_collection_ = Parse(kOutputFile, item_name_id);

    std::filesystem::remove(kOutputFile);
}

void Borgelt::FitInternal(model::IDatasetStream& data_stream) {
    if (std::type_index(typeid(data_stream)) != std::type_index(typeid(CSVParser))) {
        throw std::runtime_error(
                "Integrated Borgelt algorithms currently are not supported via Python interface");
    }
    CSVParser& data_parser = dynamic_cast<CSVParser&>(data_stream);
    sep_ = data_parser.GetSeparator();
    path_ = data_parser.GetPath();

    switch (input_format_) {
        case InputFormat::singular:
            throw std::runtime_error(
                    "Integrated Borgelt algorithms currently support only tabular format");
        case InputFormat::tabular:
            transactional_data_ =
                    model::TransactionalData::CreateFromTabular(data_stream, first_column_tid_);
            break;
        default:
            assert(0);
    }
    if (transactional_data_->GetNumTransactions() == 0) {
        throw std::runtime_error("Got an empty dataset: AR mining is meaningless.");
    }
}

void Borgelt::MakeExecuteOptsAvailable() {
    ARAlgorithm::MakeExecuteOptsAvailable();
    MakeOptionsAvailable(config::GetOptionNames(AlgoOpt));
}

unsigned long long Borgelt::ExecuteInternal() {
    Config config{minconf_, first_column_tid_, minsup_, sep_,
                  path_, algo_, kOutputFile};

    auto const elapsed = ExecuteExternal(config);

    ConvertAndFillResults();

    return elapsed;
}

Borgelt::~Borgelt() {
    std::filesystem::remove(kOutputFile);
}

}  // namespace algos