#include "algorithms/typo_miner.h"

#include <typeindex>
#include <typeinfo>

#include "util/config/equal_nulls/option.h"
#include "util/config/error/option.h"
#include "util/config/names_and_descriptions.h"
#include "util/config/option_using.h"
#include "util/config/tabular_data/input_table/option.h"

namespace algos {

static const std::type_index void_index = typeid(void);
static const std::type_info& error_type = typeid(util::config::ErrorType);

static std::unique_ptr<FDAlgorithm> GetAlgorithm(AlgorithmType algorithm) {
    if (IsDerived<PliBasedFDAlgorithm>(algorithm)) {
        return CreateAlgorithmInstance<PliBasedFDAlgorithm>(algorithm, true);
    }
    return CreateAlgorithmInstance<FDAlgorithm>(algorithm);
}

util::config::Configuration::FuncTuple TypoMiner::MakeConfigFunctions() {
    using util::config::names::kPliRelation;
    using OptionNameSet = util::config::Configuration::OptionNameSet;

    // To whomever it may concern: if you are making another pipeline, please
    // figure out a way to generalize this.
    auto add_external_needed_opts = [this](OptionNameSet& needed_options) {
        if (GetCurrentStage() == +util::config::ConfigurationStage::load_data) {
            for (FDAlgorithm* algo : {precise_algo_.get(), approx_algo_.get()}) {
                // PLI-based algorithm data are loaded by TypoMiner itself on this stage.
                if (dynamic_cast<PliBasedFDAlgorithm*>(algo) != nullptr) continue;
                OptionNameSet needed_by_algo = algo->GetNeededOptions();
                needed_options.insert(needed_by_algo.begin(), needed_by_algo.end());
            }
            return;
        }
        OptionNameSet precise_options = precise_algo_->GetNeededOptions();
        OptionNameSet approx_options = approx_algo_->GetNeededOptions();
        needed_options.insert(precise_options.begin(), precise_options.end());
        needed_options.insert(approx_options.begin(), approx_options.end());
    };
    auto set_external_opt = [this](std::string_view option_name, boost::any const& value) {
        if (option_name == kPliRelation) {
            throw std::logic_error("You do not manage this option.");
        }
        if (option_name == util::config::ErrorOpt.GetName()) {
            if (value.type() != error_type) {
                throw std::invalid_argument("Incorrect error type.");
            }
            if (value.empty()) {
                throw std::invalid_argument("Must specify error value when mining typos.");
            }
            auto error = boost::any_cast<util::config::ErrorType>(value);
            if (error == 0.0) {
                throw std::invalid_argument("Typo mining with error 0 is meaningless");
            }
            // Assumes if both have an option called `config::ErrorOpt.GetName()`,
            // then these options share semantics.
            return TrySetOption(option_name, util::config::ErrorType{0.0}, value);
        }
        return TrySetOption(option_name, value, value);
    };
    auto unset_external_opt = [this](std::string_view option_name) {
        if (option_name == kPliRelation) {
            throw std::logic_error("You do not manage this option.");
        }
        precise_algo_->UnsetOption(option_name);
        approx_algo_->UnsetOption(option_name);
    };
    auto get_external_type_index = [this](std::string_view option_name) {
        if (option_name == kPliRelation) {
            throw std::logic_error("You do not manage this option.");
        }
        std::type_index precise_index = precise_algo_->GetTypeIndex(option_name);
        std::type_index approx_index = approx_algo_->GetTypeIndex(option_name);
        if (precise_index != void_index) {
            if (approx_index != void_index && !precise_algo_->NeedsOption(option_name))
                return approx_index;
            return precise_index;
        }
        return approx_index;
    };
    auto reset_external_config = [this]() {
        precise_algo_->ResetConfiguration();
        approx_algo_->ResetConfiguration();
    };

    return {add_external_needed_opts, set_external_opt, unset_external_opt, get_external_type_index,
            reset_external_config};
}

TypoMiner::TypoMiner(AlgorithmType precise, AlgorithmType approx)
    : TypoMiner(GetAlgorithm(precise), GetAlgorithm(approx)) {}

TypoMiner::TypoMiner(std::unique_ptr<FDAlgorithm> precise_algo,
                     std::unique_ptr<FDAlgorithm> approx_algo)
        : Algorithm({/*"Precise fd algorithm execution", "Approximate fd algoritm execution",
                     "Extracting fds with non-zero error"*/}, MakeConfigFunctions()),
          precise_algo_(std::move(precise_algo)),
          approx_algo_(std::move(approx_algo)) {
    ValidateAlgorithms();
    RegisterOptions();
}

void TypoMiner::RegisterOptions() {
    DESBORDANTE_OPTION_USING;

    auto radius_check = [](double radius) {
        if (!(radius == -1 || radius >= 0)) {
            throw std::invalid_argument("Radius should be greater or equal to zero or equal to -1");
        }
    };
    auto ratio_default = [this]() {
        return (relation_->GetNumRows() <= 1) ? 1 : (2.0 / relation_->GetNumRows());
    };
    auto ratio_check = [](double ratio) {
        if (!(ratio >= 0 && ratio <= 1)) {
            throw std::invalid_argument("Ratio should be between 0 and 1");
        }
    };

    RegisterInitialLoadOption(util::config::TableOpt(&input_table_));
    RegisterInitialLoadOption(util::config::EqualNullsOpt(&is_null_equal_null_));
    RegisterInitialExecOption(
            Option{&radius_, kRadius, kDRadius, -1.0}.SetValueCheck(radius_check));
    RegisterInitialExecOption(
            Option{&ratio_, kRatio, kDRatio, {ratio_default}}.SetValueCheck(ratio_check));
}

void TypoMiner::ValidateAlgorithms() {
    using util::config::ErrorType;
    using util::config::names::kError;

    if (!approx_algo_->IsInitialAtStage(kError, util::config::ConfigurationStage::execute)) {
        throw std::logic_error("Approximate algorithm must have an error option.");
    }

    std::type_index approx_error_index = approx_algo_->GetTypeIndex(kError);
    std::type_index precise_error_index = precise_algo_->GetTypeIndex(kError);
    if (approx_error_index != error_type) {
        throw std::logic_error("Unexpected error option type in the approximate algorithm.");
    }
    if (precise_error_index != void_index && precise_error_index != error_type) {
        throw std::logic_error("Unexpected error option type in the precise algorithm.");
    }
}

void TypoMiner::ResetState() {
    approx_fds_.clear();
}

std::pair<bool, std::string> TypoMiner::TrySetOption(std::string_view option_name,
                                                     boost::any const& value_precise,
                                                     boost::any const& value_approx) {
    int successes = 0;
    std::string error_text;
    for (auto [algo, value] : {std::make_pair(precise_algo_.get(), value_precise),
                               std::make_pair(approx_algo_.get(), value_approx)}) {
        if (!algo->OptionSettable(option_name)) continue;
        auto [success_algo, error_text_algo] = algo->SetOptionNoThrow(option_name, value);
        if (success_algo) ++successes;
        if (successes) continue;
        error_text = error_text_algo;
    }
    if (successes) error_text = "";
    // Options with the same names must have certain semantics. The
    // situation where both algorithms have options with the same name but
    // one throws an exception for the pair of values while the other does
    // not is undefined behaviour. The situation where both algorithms may
    // need options with the same name, but the name is returned in one of
    // the GetNeededOptions calls and then in another is undefined behaviour
    // likewise.
    assert(!(successes == 1 && (precise_algo_->OptionSettable(option_name) ==
                                approx_algo_->OptionSettable(option_name))));
    return {successes, error_text};
}

void TypoMiner::LoadDataInternal() {
    relation_ = ColumnLayoutRelationData::CreateFrom(*input_table_, is_null_equal_null_);
    input_table_->Reset();
    typed_relation_ =
            model::ColumnLayoutTypedRelationData::CreateFrom(*input_table_, is_null_equal_null_);
    for (Algorithm* algo : {precise_algo_.get(), approx_algo_.get()}) {
        if (dynamic_cast<PliBasedFDAlgorithm*>(algo) == nullptr) {
            input_table_->Reset();
        } else {
            algo->SetOption(util::config::names::kPliRelation, relation_);
        }
        algo->LoadData();
    }
}

unsigned long long TypoMiner::ExecuteInternal() {
    auto const start_time = std::chrono::system_clock::now();

    precise_algo_->Execute();
    approx_algo_->Execute();

    std::list<FD>& precise_fds = precise_algo_->FdList();
    std::list<FD>& approx_fds = approx_algo_->FdList();

    precise_fds.sort(FDLess);
    approx_fds.sort(FDLess);

    std::set_difference(approx_fds.begin(), approx_fds.end(), precise_fds.begin(),
                        precise_fds.end(), std::back_inserter(approx_fds_), FDLess);

    auto const elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now() - start_time);

    return elapsed_milliseconds.count();
}

auto TypoMiner::MakeTuplesByIndicesComparator(std::map<int, unsigned> const& frequency_map) {
    return [&frequency_map](int const l, int const r) {
        return frequency_map.at(l) < frequency_map.at(r);
    };
}

std::vector<util::PLI::Cluster> TypoMiner::FindClustersWithTypos(FD const& typos_fd,
                                                                 bool const sort_clusters) const {
    std::vector<util::PLI::Cluster> clusters;
    std::shared_ptr<util::PLI const> intersection_pli;
    std::vector<Column const*> const lhs_columns = typos_fd.GetLhs().GetColumns();
    std::vector<int> const& probing_table =
        relation_->GetColumnData(typos_fd.GetRhs().GetIndex()).GetProbingTable();
    auto const sort_cluster = [this, &typos_fd](util::PLI::Cluster& cluster) {
        std::map<int, unsigned> const frequency_map =
            CreateFrequencyMap(typos_fd.GetRhs(), cluster);
        std::stable_sort(cluster.begin(), cluster.end(),
                         MakeTuplesByIndicesComparator(frequency_map));
    };

    if (lhs_columns.empty()) {
        /* If an approximate fd []->rhs holds then it is implied that rhs contains
         * only equal values with some deviations, which are considered errors,
         * so there is only one 'cluster with typos' containing all rows */
        util::PLI::Cluster cluster(relation_->GetNumRows());
        std::iota(cluster.begin(), cluster.end(), 0);
        if (sort_clusters) {
            sort_cluster(cluster);
        }
        return {cluster};
    }

    for (Column const* col : lhs_columns) {
        ColumnData const& col_data = relation_->GetColumnData(col->GetIndex());
        std::shared_ptr<util::PLI const> pli = col_data.GetPliOwnership();

        if (intersection_pli == nullptr) {
            intersection_pli = pli;
        } else {
            intersection_pli = intersection_pli->Intersect(pli.get());
        }
    }

    for (util::PLI::Cluster const& cluster : intersection_pli->GetIndex()) {
        int cluster_rhs_value = -1;

        /* Check if fd has wrong rhs values in this cluster */
        for (int const tuple_index : cluster) {
            int const probing_table_value = probing_table[tuple_index];

            if (cluster_rhs_value == -1) {
                cluster_rhs_value = probing_table_value;
            } else if (cluster_rhs_value != probing_table_value) {
                cluster_rhs_value = -1;
                break;
            }
        }

        if (cluster_rhs_value == -1 ||
            (ColumnData::IsValueSingleton(cluster_rhs_value) && cluster.size() != 1)) {
            /* Actually intersection_pli is owned only by this method most of the time
             * (when lhs_columns.size() != 1), so we can here move cluster when
             * lhs_columns.size() != 1. But that leads to more cumbersome and complex code.
             * So I decided to leave it as it is until we know for sure that this place causes
             * performance problems.
             */
            clusters.push_back(cluster);

            if (sort_clusters) {
                sort_cluster(clusters.back());
            }
        }
    }

    return clusters;
}

std::vector<TypoMiner::SquashedElement> TypoMiner::SquashCluster(
        FD const& squash_on, util::PLI::Cluster const& cluster) const {
    std::vector<SquashedElement> squashed;
    std::vector<int> const& probing_table =
            relation_->GetColumnData(squash_on.GetRhs().GetIndex()).GetProbingTable();

    if (cluster.empty()) {
        return squashed;
    }

    auto prev = cluster.begin();
    squashed.push_back({.tuple_index = *prev, .amount = 1});

    for (auto it = std::next(cluster.cbegin()); it != cluster.cend(); ++it) {
        if (probing_table[*it] != util::PLI::singleton_value_id_ &&
            probing_table[*it] == probing_table[*prev]) {
            squashed.back().amount++;
        } else {
            squashed.push_back({.tuple_index = *it, .amount = 1});
        }
        prev = it;
    }

    return squashed;
}

void TypoMiner::SortCluster(FD const& sort_on, util::PLI::Cluster& cluster) const {
    std::map<int, unsigned> const frequency_map = CreateFrequencyMap(sort_on.GetRhs(), cluster);

    std::stable_sort(cluster.begin(), cluster.end(), MakeTuplesByIndicesComparator(frequency_map));
}

void TypoMiner::RestoreLineOrder(util::PLI::Cluster& cluster) const {
    std::sort(cluster.begin(), cluster.end());
}

void TypoMiner::RestoreLineOrder(std::vector<TypoMiner::SquashedElement>& squashed_cluster) const {
    std::sort(squashed_cluster.begin(), squashed_cluster.end(),
              [](const TypoMiner::SquashedElement& lhs, const TypoMiner::SquashedElement& rhs) {
                  return lhs.tuple_index < rhs.tuple_index;
              });
}

std::vector<util::PLI::Cluster::value_type> TypoMiner::FindLinesWithTypos(
        FD const& typos_fd, util::PLI::Cluster const& cluster, double new_radius,
        double new_ratio) {
    SetRadius(new_radius);
    SetRatio(new_ratio);
    return FindLinesWithTypos(typos_fd, cluster);
}

std::vector<util::PLI::Cluster::value_type> TypoMiner::FindLinesWithTypos(
    FD const& typos_fd, util::PLI::Cluster const& cluster) const {
    Column const& col = typos_fd.GetRhs();
    model::TypedColumnData const& col_data = typed_relation_->GetColumnData(col.GetIndex());
    std::vector<int> const& probing_table =
        relation_->GetColumnData(col.GetIndex()).GetProbingTable();
    model::Type const& type = col_data.GetType();

    unsigned most_freq_index = GetMostFrequentValueIndex(col, cluster);
    int most_freq_value = probing_table[most_freq_index];

    if (ColumnData::IsValueSingleton(most_freq_value) || col_data.IsMixed() ||
        !type.IsMetrizable() || col_data.IsEmpty(most_freq_index) ||
        col_data.IsNullOrEmpty(most_freq_index)) {
        if (ratio_ == 1) {
            return cluster;
        }
        return {};
    }

    std::vector<util::PLI::Cluster::value_type> typos;
    unsigned long num_of_close_values = 0;
    std::vector<std::byte const*> const& data = col_data.GetData();

    for (util::PLI::Cluster::value_type tuple_index : cluster) {
        /* Temporary ignoring NULL or empty values. Maybe it should be decided by some parameter
         * if NULL (empty) value is close to any non-NULL (non-empty) or not (as in metric
         * dependecies).
         */
        if (most_freq_value == probing_table[tuple_index] || col_data.IsNullOrEmpty(tuple_index)) {
            continue;
        }
        if (radius_ == -1 || ValuesAreClose(data[most_freq_index], data[tuple_index], type)) {
            num_of_close_values++;
            typos.push_back(tuple_index);
        }
    }

    if (double(num_of_close_values) / col_data.GetNumRows() > ratio_) {
        return {};
    }

    return typos;
}

std::vector<TypoMiner::ClusterTyposPair> TypoMiner::FindClustersAndLinesWithTypos(
        const FD& typos_fd, const bool sort_clusters) const {
    std::vector<ClusterTyposPair> result;
    std::vector<util::PLI::Cluster> clusters = FindClustersWithTypos(typos_fd, sort_clusters);

    result.reserve(clusters.size());

    for (auto& cluster : clusters) {
        TyposVec typos = FindLinesWithTypos(typos_fd, cluster);
        if (!typos.empty()) {
            result.emplace_back(std::move(cluster), std::move(typos));
        }
    }

    return result;
}

unsigned TypoMiner::GetMostFrequentValueIndex(Column const& cluster_col,
                                              util::PLI::Cluster const& cluster) const {
    assert(!cluster.empty());
    std::vector<int> const& probing_table =
        relation_->GetColumnData(cluster_col.GetIndex()).GetProbingTable();
    std::unordered_map<int, unsigned> frequencies =
            util::PLI::CreateFrequencies(cluster, probing_table);

    unsigned most_frequent_index = cluster.size();
    unsigned largest_frequency = 0;
    for (int const tuple_index : cluster) {
        int const probing_table_value = probing_table[tuple_index];
        unsigned const frequency = (ColumnData::IsValueSingleton(probing_table_value))
                                   ? 1
                                   : frequencies.at(probing_table_value);
        if (frequency > largest_frequency) {
            largest_frequency = frequency;
            most_frequent_index = tuple_index;
        }
    }

    return most_frequent_index;
}

std::map<int, unsigned> TypoMiner::CreateFrequencyMap(Column const& cluster_col,
                                                      util::PLI::Cluster const& cluster) const {
    std::map<int, unsigned> frequency_map;
    std::vector<int> const& probing_table =
        relation_->GetColumnData(cluster_col.GetIndex()).GetProbingTable();
    std::unordered_map<int, unsigned> frequencies =
            util::PLI::CreateFrequencies(cluster, probing_table);

    for (int const tuple_index : cluster) {
        int const probing_table_value = probing_table[tuple_index];
        unsigned const value = (ColumnData::IsValueSingleton(probing_table_value))
                                   ? 1
                                   : frequencies.at(probing_table_value);
        frequency_map[tuple_index] = value;
    }

    return frequency_map;
}

bool TypoMiner::FDLess(FD const& l, FD const& r) {
    if (l.GetLhs() < r.GetLhs()) {
        return true;
    } else if (l.GetLhs() == r.GetLhs()) {
        return l.GetRhs() < r.GetRhs();
    }

    return false;
}

}  // namespace algos
