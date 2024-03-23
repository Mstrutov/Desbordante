#include "algorithms/dd//split/split.h"

#include <cassert>
#include <chrono>
#include <cstddef>
#include <deque>
#include <list>
#include <regex>
#include <set>
#include <vector>

#include <easylogging++.h>

#include "config/names_and_descriptions.h"
#include "config/option_using.h"
#include "config/tabular_data/input_table/option.h"
#include "model/types/numeric_type.h"
#include "util/levenshtein_distance.h"

namespace algos::dd {

Split::Split() : Algorithm({}) {
    using namespace config::names;
    RegisterOptions();
    MakeOptionsAvailable({config::TableOpt.GetName()});
}

void Split::RegisterOptions() {
    DESBORDANTE_OPTION_USING;

    config::InputTable default_table;

    RegisterOption(config::TableOpt(&input_table_));
    RegisterOption(Option{&difference_table_, kDifferenceTable, kDDifferenceTable, default_table});
    RegisterOption(Option{&num_rows_, kNumRows, kDNumRows, 0U});
    RegisterOption(Option{&num_columns_, kNumColumns, kDNUmColumns, 0U});
}

void Split::MakeExecuteOptsAvailable() {
    using namespace config::names;

    MakeOptionsAvailable({kDifferenceTable, kNumRows, kNumColumns});
}

void Split::LoadDataInternal() {
    relation_ = ColumnLayoutRelationData::CreateFrom(*input_table_, false);  // nulls are
                                                                             // ignored
    input_table_->Reset();
    typed_relation_ = model::ColumnLayoutTypedRelationData::CreateFrom(*input_table_,
                                                                       false);  // nulls are ignored
}

void Split::SetLimits() {
    if (num_rows_ == 0) num_rows_ = typed_relation_->GetNumRows();
    if (num_columns_ == 0) num_columns_ = typed_relation_->GetNumColumns();
}

void Split::ParseDifferenceTable() {
    has_dif_table_ = (difference_table_.get() != nullptr);

    if (has_dif_table_) {
        difference_typed_relation_ =
                model::ColumnLayoutTypedRelationData::CreateFrom(*difference_table_,
                                                                 false);  // nulls are ignored
        assert(typed_relation_->GetNumColumns() == difference_typed_relation_->GetNumColumns());
    }
}

unsigned long long Split::ExecuteInternal() {
    SetLimits();
    ParseDifferenceTable();

    auto start_time = std::chrono::system_clock::now();
    LOG(DEBUG) << "Start";

    CalculateAllDistances();

    if (reduce_method_ == +Reduce::IEHybrid) {
        CalculateTuplePairs();
    }

    LOG(DEBUG) << "Calculated distances";
    auto elapsed_milliseconds1 = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_time);
    LOG(DEBUG) << elapsed_milliseconds1.count();

    for (std::size_t index = 0; index < num_columns_; index++)
        LOG(DEBUG) << min_max_dif_[index].first << " " << min_max_dif_[index].second;

    unsigned search_size = ReduceDDs(start_time);

    LOG(DEBUG) << "Reduced dependencies";

    unsigned num_cycles = RemoveRedundantDDs();

    LOG(DEBUG) << "Removed redundant dependencies";
    LOG(DEBUG) << "Cycles: " << num_cycles;

    num_cycles = RemoveTransitiveDDs();

    LOG(DEBUG) << "Removed transitive dependencies";
    LOG(DEBUG) << "Cycles: " << num_cycles;
    LOG(DEBUG) << "Search space size: " << search_size;

    PrintResults();

    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_time);
    LOG(DEBUG) << "Algorithm time: " << elapsed_milliseconds.count();
    return elapsed_milliseconds.count();
}

unsigned Split::ReduceDDs(auto start_time) {
    unsigned search_size = 0;
    std::list<df> search, dfs_y;
    std::list<dd> reduced;

    auto elapsed_milliseconds1 = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_time);

    for (std::size_t index = 0; index < num_columns_; index++) {
        std::list<std::size_t> indices;
        for (std::size_t j = 0; j < num_columns_; j++) {
            if (j != index) indices.emplace_back(j);
        }

        search = SearchSpace(indices);
        dfs_y = SearchSpace(index);
        search_size += search.size() * (dfs_y.size() - 1);
        elapsed_milliseconds1 = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::system_clock::now() - start_time);
        LOG(DEBUG) << elapsed_milliseconds1.count() << " " << search.size() << " " << dfs_y.size();

        for (auto& df_y : dfs_y) {
            unsigned cnt = 0;
            if (df_y[index] != min_max_dif_[index]) {
                switch (reduce_method_) {
                    case +Reduce::Negative:
                        reduced = NegativePruningReduce(df_y, search, cnt);
                        break;
                    case +Reduce::Hybrid:
                        reduced = HybridPruningReduce(df_y, search, cnt);
                        break;
                    case +Reduce::IEHybrid:
                        reduced = InstanceExclusionReduce(tuple_pairs_, search, df_y, cnt);
                        break;
                    default:
                        break;
                }
                for (auto& reduced_df : reduced) results_.push_back(reduced_df);
            }
            LOG(DEBUG) << cnt;
        }
        elapsed_milliseconds1 = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::system_clock::now() - start_time);
        LOG(DEBUG) << "Current time: " << elapsed_milliseconds1.count();
    }
    return search_size;
}

unsigned Split::RemoveRedundantDDs() {
    unsigned num_cycles = 0;
    std::list<dd> results_copy;
    dd ldd, rdd;

    while (true) {
        num_cycles++;
        results_copy.clear();
        std::size_t left_index = 0;
        for (auto& left_dd : results_) {
            bool h = true;
            std::size_t right_index = 0;
            for (auto& right_dd : results_) {
                if (left_index != right_index) {
                    if (Subsume(right_dd.lhs, left_dd.lhs)) {
                        if (Subsume(left_dd.rhs, right_dd.rhs)) {
                            h = false;
                            break;
                        }
                    }
                }
                right_index++;
            }
            if (h) results_copy.push_back(left_dd);
            left_index++;
        }
        if (results_copy.size() == results_.size()) break;
        results_ = results_copy;
    }
    return num_cycles;
}

unsigned Split::RemoveTransitiveDDs() {
    unsigned num_cycles = 0;
    std::list<dd> results_copy;

    while (true) {
        num_cycles++;
        results_copy.clear();
        bool is_removable = false;
        for (auto& dd3 : results_) {
            bool remove = false;
            for (auto& dd1 : results_) {
                for (auto& dd2 : results_) {
                    if (Subsume(dd2.lhs, dd1.rhs) && dd1.lhs == dd3.lhs && dd2.rhs == dd3.rhs) {
                        if (!is_removable) remove = true;
                        is_removable = true;
                        break;
                    }
                }
                if (is_removable) break;
            }
            if (!remove) results_copy.push_back(dd3);
        }
        if (results_copy.size() == results_.size()) break;
        results_ = results_copy;
    }
    return num_cycles;
}

double Split::CalculateDistance(std::size_t column_index,
                                std::pair<std::size_t, std::size_t> tuple_pair) {
    model::TypedColumnData const& column = typed_relation_->GetColumnData(column_index);
    model::TypeId type_id = column.GetTypeId();

    if (type_id == +model::TypeId::kUndefined) {
        throw std::invalid_argument("Column with index \"" + std::to_string(column_index) +
                                    "\" type undefined.");
    }
    if (type_id == +model::TypeId::kMixed) {
        throw std::invalid_argument("Column with index \"" + std::to_string(column_index) +
                                    "\" contains values of different types.");
    }
    if (column.IsNull(tuple_pair.first) || column.IsNull(tuple_pair.second)) {
        throw std::runtime_error("Some of the value coordinates are nulls.");
    }
    if (column.IsEmpty(tuple_pair.first) || column.IsEmpty(tuple_pair.second)) {
        throw std::runtime_error("Some of the value coordinates are empties.");
    }
    double dif = 0;
    if (column.GetType().IsMetrizable()) {
        std::byte const* first_value = column.GetValue(tuple_pair.first);
        std::byte const* second_value = column.GetValue(tuple_pair.second);
        auto const& type = static_cast<model::IMetrizableType const&>(column.GetType());
        dif = type.Dist(first_value, second_value);
    }
    return dif;
}

// must be inline for optimization (gcc 11.4.0)
inline bool Split::CheckDF(df const& dif_func, std::pair<std::size_t, std::size_t> tuple_pair) {
    for (std::size_t column_index = 0; column_index < num_columns_; column_index++) {
        double dif = distances_[column_index][tuple_pair.first][tuple_pair.second];
        if (dif < dif_func[column_index].first || dif > dif_func[column_index].second) {
            return false;
        }
    }
    return true;
}

bool Split::VerifyDD(dd const& dep) {
    for (std::size_t i = 0; i < num_rows_; i++) {
        for (std::size_t j = i + 1; j < num_rows_; j++) {
            if (CheckDF(dep.lhs, {i, j}) && !CheckDF(dep.rhs, {i, j})) return false;
        }
    }
    return true;
}

void Split::CalculateAllDistances() {
    for (std::size_t column_index = 0; column_index < num_columns_; column_index++) {
        distances_.push_back({});

        std::shared_ptr<model::PLI const> pli =
                relation_->GetColumnData(column_index).GetPliOwnership();
        std::deque<model::PLI::Cluster> const& index = pli->GetIndex();
        std::shared_ptr<std::vector<int> const> probing_table = pli->CalculateAndGetProbingTable();
        model::PLI::Cluster const& pt = *probing_table.get();
        for (std::size_t i = 0; i < num_rows_; i++) {
            distances_[column_index].push_back({});
            for (std::size_t j = 0; j < num_rows_; j++) {
                distances_[column_index][i].push_back(0);
            }
        }
        double max_dif = 0, min_dif = 1000000000;
        for (std::size_t i = 0; i < index.size(); i++) {
            for (std::size_t j = 0; j < index.size(); j++) {
                std::size_t first_index = index[i][0];
                std::size_t second_index = index[j][0];
                if (first_index < num_rows_ && second_index < num_rows_) {
                    if (first_index < second_index) {
                        double dif;
                        dif = CalculateDistance(column_index, {first_index, second_index});
                        max_dif = std::max(max_dif, dif);
                        min_dif = std::min(min_dif, dif);
                        distances_[column_index][first_index][second_index] = dif;
                    } else if (first_index == second_index)
                        distances_[column_index][first_index][second_index] = 0;
                    else
                        distances_[column_index][first_index][second_index] =
                                distances_[column_index][second_index][first_index];
                }
            }
        }
        for (std::size_t i = 0; i < num_rows_; i++) {
            for (std::size_t j = 0; j < num_rows_; j++) {
                if (pt[i] != 0 && pt[j] != 0) {
                    distances_[column_index][i][j] =
                            distances_[column_index][index[pt[i] - 1][0]][index[pt[j] - 1][0]];
                    if (i != j && pt[i] == pt[j]) min_dif = 0;
                } else {
                    if (i < j) {
                        double dif;
                        dif = CalculateDistance(column_index, {i, j});
                        max_dif = std::max(max_dif, dif);
                        min_dif = std::min(min_dif, dif);
                        distances_[column_index][i][j] = dif;
                    } else if (i == j)
                        distances_[column_index][i][j] = 0;
                    else
                        distances_[column_index][i][j] = distances_[column_index][j][i];
                }
            }
        }
        min_max_dif_.push_back({min_dif, max_dif});
    }
}

bool Split::IsFeasible(df const& d) {
    for (std::size_t i = 0; i < num_rows_; i++) {
        for (std::size_t j = i + 1; j < num_rows_; j++) {
            if (CheckDF(d, {i, j})) return true;
        }
    }
    return false;
}

std::list<df> Split::SearchSpace(std::size_t index) {
    std::list<df> dfs;
    df d = min_max_dif_;
    dfs.push_back(d);

    if (!has_dif_table_) {
        // differential functions should be put in this exact order for further reducing
        for (int i = 4; i >= 0; i--) {
            if (i >= min_max_dif_[index].first && i < min_max_dif_[index].second) {
                d[index] = {min_max_dif_[index].first, i};
                dfs.push_back(d);
            }
        }
        return dfs;
    }

    std::size_t dif_num_rows = difference_typed_relation_->GetNumRows();

    model::TypedColumnData const& dif_column = difference_typed_relation_->GetColumnData(index);
    model::Type const& type = dif_column.GetType();

    auto pair_compare = [](std::pair<double, double> const& first_pair,
                           std::pair<double, double> const& second_pair) {
        return ((first_pair.second - first_pair.first) >
                (second_pair.second - second_pair.first)) ||
               ((first_pair.second - first_pair.first) ==
                        (second_pair.second - second_pair.first) &&
                (first_pair.first > second_pair.first));
    };

    std::set<std::pair<double, double>, decltype(pair_compare)> limits(pair_compare);

    std::regex df_regex(R"(\[\d{1,19}(\.\d*)?\;\d{1,19}(\.\d*)?\]$)");

    for (std::size_t row_index = 0; row_index < dif_num_rows; row_index++) {
        model::TypeId type_id = dif_column.GetValueTypeId(row_index);
        if (type_id == +model::TypeId::kString) {
            std::string df_str;
            if (dif_column.IsMixed()) {
                model::MixedType const* mixed_type = dif_column.GetIfMixed();
                std::byte const* df_string = dif_column.GetValue(row_index);
                std::unique_ptr<model::Type> t = mixed_type->RetrieveType(df_string);
                std::byte const* df_string_value = mixed_type->RetrieveValue(df_string);
                model::StringType const* str_type = static_cast<model::StringType const*>(t.get());
                df_str = str_type->ValueToString(df_string_value);
            } else {
                model::StringType const& str_type = static_cast<model::StringType const&>(type);
                std::byte const* df_string = dif_column.GetValue(row_index);
                df_str = str_type.ValueToString(df_string);
            }
            double lower_limit, upper_limit;
            if (std::regex_match(df_str, df_regex)) {
                std::size_t sep_index = df_str.find(';');
                std::string lower_str = df_str.substr(1, sep_index - 1);
                std::string upper_str = df_str.substr(sep_index + 1, df_str.size() - 2 - sep_index);
                lower_limit = model::TypeConverter<double>::convert(lower_str);
                upper_limit = model::TypeConverter<double>::convert(upper_str);

                if (upper_limit >= min_max_dif_[index].first &&
                    lower_limit <= min_max_dif_[index].second && lower_limit <= upper_limit) {
                    std::pair<double, double> intersect = {
                            std::max(lower_limit, min_max_dif_[index].first),
                            std::min(upper_limit, min_max_dif_[index].second)};
                    if (intersect != min_max_dif_[index]) limits.insert(intersect);
                }
            }
        }
    }

    // differential functions should be put in this exact order for further reducing
    for (auto limit : limits) {
        d[index] = limit;
        dfs.push_back(d);
    }
    return dfs;
}

std::list<df> Split::SearchSpace(std::list<std::size_t>& indices) {
    if (indices.size() == 1) return SearchSpace(*indices.begin());
    std::size_t last_index = *indices.rbegin();
    std::list<df> d = SearchSpace(last_index);
    indices.pop_back();
    std::list<df> e = SearchSpace(indices);
    std::list<df> f;
    df g = min_max_dif_;
    for (auto const& first_df : e) {
        for (auto const& second_df : d) {
            for (std::size_t k = 0; k < num_columns_; k++) {
                g[k] = {std::max(first_df[k].first, second_df[k].first),
                        std::min(first_df[k].second, second_df[k].second)};
            }
            if (IsFeasible(g))
                f.push_back(g);
            else
                break;
        }
    }
    return f;
}

bool Split::Subsume(df const& df1, df const& df2) {
    for (std::size_t i = 0; i < num_columns_; i++) {
        if (df2[i].first < df1[i].first || df2[i].second > df1[i].second) return false;
    }
    return true;
}

std::list<dd> Split::NegativePruningReduce(df const& rhs, std::list<df> const& search,
                                           unsigned& cnt) {
    if (!search.size()) return {};

    std::list<dd> dds;
    df last_df = *search.rbegin();

    cnt++;
    if (!VerifyDD({last_df, rhs})) {
        std::list<df> remainder;
        for (auto const& search_df : search) {
            if (search_df != last_df) {
                if (!Subsume(search_df, last_df)) {
                    remainder.push_back(search_df);
                }
            }
        }
        return NegativePruningReduce(rhs, remainder, cnt);
    }

    std::list<df> prune;
    std::list<df> remainder;
    for (auto const& search_df : search) {
        if (search_df != last_df) {
            if (Subsume(search_df, last_df)) {
                prune.push_back(search_df);
            } else {
                remainder.push_back(search_df);
            }
        }
    }

    dds = NegativePruningReduce(rhs, prune, cnt);
    if (!dds.size()) dds.push_back({last_df, rhs});
    std::size_t cur_size = dds.size();
    std::list<dd> remaining_dds = NegativePruningReduce(rhs, remainder, cnt);

    for (auto const& remaining_dd : remaining_dds) {
        bool h = true;
        std::size_t index = 0;
        for (auto const& pruning_dd : dds) {
            if (Subsume(pruning_dd.lhs, remaining_dd.lhs)) {
                h = false;
                break;
            }
            index++;
            if (index == cur_size) break;
        }
        if (h) dds.push_back(remaining_dd);
    }

    return dds;
}

std::list<dd> Split::HybridPruningReduce(df const& rhs, std::list<df> const& search,
                                         unsigned& cnt) {
    if (!search.size()) return {};

    std::list<dd> dds;
    df first_df = *search.begin();
    df last_df = *search.rbegin();

    cnt++;
    if (VerifyDD({first_df, rhs})) {
        dds.push_back({first_df, rhs});
        std::list<df> remainder;
        for (auto const& search_df : search) {
            if (search_df != first_df) {
                if (!Subsume(first_df, search_df)) {
                    remainder.push_back(search_df);
                }
            }
        }
        std::list<dd> remaining_dds = HybridPruningReduce(rhs, remainder, cnt);
        for (auto const& remaining_dd : remaining_dds) dds.push_back(remaining_dd);
        return dds;
    }

    cnt++;
    if (!VerifyDD({last_df, rhs})) {
        std::list<df> remainder;
        for (auto const& search_df : search) {
            if (search_df != last_df) {
                if (!Subsume(search_df, last_df)) {
                    remainder.push_back(search_df);
                }
            }
        }
        return HybridPruningReduce(rhs, remainder, cnt);
    }

    std::list<df> prune;
    std::list<df> remainder;
    for (auto const& search_df : search) {
        if (search_df != first_df) {
            if (Subsume(first_df, search_df)) {
                prune.push_back(search_df);
            } else {
                remainder.push_back(search_df);
            }
        }
    }

    dds = HybridPruningReduce(rhs, remainder, cnt);
    std::size_t cur_size = dds.size();
    std::list<dd> pruning_dds = HybridPruningReduce(rhs, prune, cnt);

    for (auto const& pruning_dd : pruning_dds) {
        bool h = true;
        std::size_t index = 0;
        for (auto const& remaining_dd : dds) {
            if (Subsume(remaining_dd.lhs, pruning_dd.lhs)) {
                h = false;
                break;
            }
            index++;
            if (index == cur_size) break;
        }
        if (h) dds.push_back(pruning_dd);
    }

    return dds;
}

std::list<dd> Split::InstanceExclusionReduce(
        std::list<std::pair<std::size_t, std::size_t>> const& tuple_pairs,
        std::list<df> const& search, df const& rhs, unsigned& cnt) {
    if (!search.size()) return {};

    std::list<dd> dds;
    df first_df = *search.begin();
    df last_df = *search.rbegin();
    std::list<std::pair<std::size_t, std::size_t>> remaining_tuple_pairs;

    cnt++;
    for (auto pair : tuple_pairs) {
        if (CheckDF(first_df, pair) && !CheckDF(rhs, pair)) remaining_tuple_pairs.push_back(pair);
    }

    if (!remaining_tuple_pairs.size()) {
        dds.push_back({first_df, rhs});
        std::list<df> remainder;
        for (auto const& search_df : search) {
            if (search_df != first_df) {
                if (!Subsume(first_df, search_df)) {
                    remainder.push_back(search_df);
                }
            }
        }
        std::list<dd> remaining_dds = InstanceExclusionReduce(tuple_pairs, remainder, rhs, cnt);
        for (auto const& remaining_dd : remaining_dds) dds.push_back(remaining_dd);
        return dds;
    }

    std::list<std::pair<std::size_t, std::size_t>> other_remaining_tuple_pairs;

    cnt++;
    for (auto pair : tuple_pairs) {
        if (CheckDF(last_df, pair) && !CheckDF(rhs, pair))
            other_remaining_tuple_pairs.push_back(pair);
    }

    if (other_remaining_tuple_pairs.size()) {
        std::list<df> remainder;
        for (auto const& search_df : search) {
            if (search_df != last_df) {
                if (!Subsume(search_df, last_df)) {
                    remainder.push_back(search_df);
                }
            }
        }
        return InstanceExclusionReduce(tuple_pairs, remainder, rhs, cnt);
    }

    std::list<df> prune;
    std::list<df> remainder;
    for (auto const& search_df : search) {
        if (search_df != first_df) {
            if (Subsume(first_df, search_df)) {
                prune.push_back(search_df);
            } else {
                remainder.push_back(search_df);
            }
        }
    }

    dds = InstanceExclusionReduce(tuple_pairs, remainder, rhs, cnt);
    std::size_t cur_size = dds.size();
    std::list<dd> pruning_dds = InstanceExclusionReduce(remaining_tuple_pairs, prune, rhs, cnt);

    for (auto const& pruning_dd : pruning_dds) {
        bool h = true;
        std::size_t index = 0;
        for (auto const& remaining_dd : dds) {
            if (Subsume(remaining_dd.lhs, pruning_dd.lhs)) {
                h = false;
                break;
            }
            index++;
            if (index == cur_size) break;
        }
        if (h) dds.push_back(pruning_dd);
    }

    return dds;
}

void Split::CalculateTuplePairs() {
    for (std::size_t i = 0; i < num_rows_; i++) {
        for (std::size_t j = i + 1; j < num_rows_; j++) {
            tuple_pairs_.push_back({i, j});
        }
    }
}

void Split::PrintResults() {
    std::list<model::DDString> result_strings = GetDDStringList();
    LOG(DEBUG) << "Minimal cover size: " << result_strings.size();
    for (auto const& result_str : result_strings) {
        LOG(DEBUG) << result_str.ToString();
    }
}

std::list<dd> Split::GetResults() {
    return results_;
}

std::vector<std::pair<double, double>> Split::GetMinMaxDif() {
    return min_max_dif_;
}

model::DDString Split::DDToDDString(dd const& dd) {
    model::DDString dd_string;
    for (std::size_t index = 0; index < num_columns_; index++) {
        if (dd.lhs[index] != min_max_dif_[index]) {
            dd_string.left.emplace_back(input_table_->GetColumnName(index), dd.lhs[index].first,
                                        dd.lhs[index].second);
        }
    }

    for (std::size_t index = 0; index < num_columns_; index++) {
        if (dd.rhs[index] != min_max_dif_[index]) {
            dd_string.right.emplace_back(input_table_->GetColumnName(index), dd.rhs[index].first,
                                         dd.rhs[index].second);
        }
    }
    return dd_string;
}

std::list<model::DDString> Split::GetDDStringList() {
    std::list<model::DDString> dd_strings;
    for (auto const& result_dd : results_) {
        dd_strings.emplace_back(DDToDDString(result_dd));
    }
    return dd_strings;
}

}  // namespace algos::dd
