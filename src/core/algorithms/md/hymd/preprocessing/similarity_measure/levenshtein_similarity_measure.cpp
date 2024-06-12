#include "algorithms/md/hymd/preprocessing/similarity_measure/levenshtein_similarity_measure.h"

#include <atomic>
#include <cstddef>
#include <numeric>
#include <set>

#include "algorithms/md/hymd/lowest_bound.h"
#include "algorithms/md/hymd/utility/make_unique_for_overwrite.h"
#include "config/exceptions.h"
#include "model/types/double_type.h"
#include "model/types/string_type.h"
#include "util/pick_m_highest_bias.h"

namespace {
using namespace algos::hymd;
using SimValPair = std::pair<preprocessing::Similarity, ValueIdentifier>;

struct SimTaskData {
    std::vector<model::md::DecisionBoundary> row_decision_bounds;
    model::md::DecisionBoundary row_lowest = 1.0;
    std::vector<SimValPair> sim_value_id_vec;
    std::size_t valid_records_number = 0;

    indexes::SimilarityMatrixRow sim_matrix_row;
    indexes::MatchingRecsMapping matching_recs_mapping;
};

std::size_t GetLevenshteinBufferSize(auto const& right_string) noexcept {
    return right_string.size() + 1;
}

/* An optimized version of the Levenshtein distance computation algorithm from
 * https://en.wikipedia.org/wiki/Levenshtein_distance, using preallocated buffers
 */
unsigned LevenshteinDistance(auto const* l_ptr, auto const* r_ptr, unsigned* v0, unsigned* v1,
                             std::size_t max_dist, std::size_t bad_value) noexcept {
    std::size_t r_size = r_ptr->size();
    assert(v0 < v1);
    assert(GetLevenshteinBufferSize(*r_ptr) == std::size_t(v1 - v0));
    std::size_t l_size = l_ptr->size();
    if (r_size > l_size) {
        std::swap(l_ptr, r_ptr);
        std::swap(l_size, r_size);
    }
    if (l_size - r_size > max_dist) {
        return bad_value;
    }

    auto const& l = *l_ptr;
    auto const& r = *r_ptr;

    std::iota(v0, v0 + r_size + 1, 0);

    auto compute_arrays = [&](auto* v0, auto* v1, unsigned i) {
        *v1 = i + 1;
        auto const li = l[i];

        for (unsigned j = 0; j != r_size;) {
            unsigned const insert_cost = v1[j] + 1;
            unsigned const substition_cost = v0[j] + (li != r[j]);
            ++j;
            unsigned const deletion_cost = v0[j] + 1;

            v1[j] = std::min({deletion_cost, insert_cost, substition_cost});
        }
    };
    auto loop_to_l_size = [&l_size, &v0, &v1, &compute_arrays]() {
        for (unsigned i = 0; i != l_size; ++i) {
            compute_arrays(v0, v1, i);
            ++i;
            compute_arrays(v1, v0, i);
        }
    };
    if (l_size & 1) {
        --l_size;
        loop_to_l_size();
        compute_arrays(v0, v1, l_size);
        return v1[r_size];
    } else {
        loop_to_l_size();
        return v0[r_size];
    }
}
}  // namespace

namespace algos::hymd::preprocessing::similarity_measure {

indexes::SimilarityMeasureOutput LevenshteinSimilarityMeasure::MakeIndexes(
        std::shared_ptr<DataInfo const> data_info_left,
        std::shared_ptr<DataInfo const> data_info_right,
        std::vector<indexes::PliCluster> const& clusters_right,
        util::WorkerThreadPool& thread_pool) const {
    std::set<model::md::DecisionBoundary> decision_bounds_set;
    indexes::SimilarityMatrix similarity_matrix;
    indexes::SimilarityIndex similarity_index;
    std::size_t const data_left_size = data_info_left->GetElementNumber();
    std::size_t const data_right_size = data_info_right->GetElementNumber();
    Similarity lowest = 1.0;
    auto const& left_nulls = data_info_left->GetNulls();
    auto const& left_empty = data_info_left->GetEmpty();
    auto const& right_nulls = data_info_right->GetNulls();
    auto const& right_empty = data_info_right->GetEmpty();
    std::vector<SimTaskData> task_info(data_left_size);
    auto process_value_id = [&](ValueIdentifier const value_id_left) {
        auto simple_case = [&](SimTaskData& data, auto const& collection) {
            std::size_t const collection_size = collection.size();
            assert(data_right_size != 0);
            if (collection_size != data_right_size) {
                data.row_lowest = kLowestBound;
            }
            data.row_decision_bounds.assign(collection_size, 1.0);
            data.sim_value_id_vec.reserve(collection_size);
            for (ValueIdentifier value_id_right : collection) {
                data.valid_records_number += clusters_right[value_id_right].size();
                data.sim_value_id_vec.emplace_back(1.0, value_id_right);
            }

            if (data.sim_value_id_vec.empty()) return;
            auto& rec_set = data.matching_recs_mapping[1.0];
            rec_set.reserve(data.valid_records_number);
            for (ValueIdentifier value_id_right : collection) {
                data.sim_matrix_row[value_id_right] = 1.0;
                indexes::PliCluster const& cluster = clusters_right[value_id_right];
                rec_set.insert(cluster.begin(), cluster.end());
            }
        };
        if (left_nulls.find(value_id_left) != left_nulls.end()) {
            if (!is_null_equal_null_) task_info[value_id_left].row_lowest = kLowestBound;
            simple_case(task_info[value_id_left], right_nulls);
        } else if (left_empty.find(value_id_left) != left_empty.end()) {
            simple_case(task_info[value_id_left], right_empty);
        } else {
            auto const& string_left =
                    model::Type::GetValue<std::string>(data_info_left->GetAt(value_id_left));
            std::size_t const left_size = string_left.size();

            std::size_t const buf_size = GetLevenshteinBufferSize(string_left);
            auto buf = /* TODO: replace with std::make_unique_for_overwrite when GCC in CI is
                          upgraded */
                    utility::MakeUniqueForOverwrite<unsigned[]>(buf_size * 2);

            auto get_similarity = [this, &string_left, left_size, &data_info_right,
                                   buf1 = buf.get(),
                                   buf2 = buf.get() + buf_size](ValueIdentifier value_id_right) {
                auto const& right_nulls = data_info_right->GetNulls();
                if (right_nulls.find(value_id_right) != right_nulls.end()) return kLowestBound;
                auto const& right_empty = data_info_right->GetEmpty();
                if (right_empty.find(value_id_right) != right_empty.end()) return kLowestBound;

                auto const& string_right =
                        model::Type::GetValue<std::string>(data_info_right->GetAt(value_id_right));
                std::size_t const max_dist = std::max(left_size, string_right.size());
                // Left has to be second since that's what the function uses to determine the buffer
                // size it needs
                Similarity value =
                        static_cast<Similarity>(
                                max_dist - LevenshteinDistance(&string_right, &string_left, buf1,
                                                               buf2, max_dist * (1 - min_sim_), max_dist)) /
                        static_cast<Similarity>(max_dist);
                return value;
            };

            SimTaskData& data = task_info[value_id_left];
            assert(data_right_size > 0);
            for (ValueIdentifier value_id_right = 0; value_id_right < data_right_size;
                 ++value_id_right) {
                Similarity similarity = get_similarity(value_id_right);
                if (similarity < min_sim_) {
                    // Metanome keeps the actual value for some reason.
                    data.row_lowest = kLowestBound /*similarity???*/;
                    continue;
                }
                if (data.row_lowest > similarity) data.row_lowest = similarity;
                data.sim_value_id_vec.emplace_back(similarity, value_id_right);
                data.valid_records_number += clusters_right[value_id_right].size();
                data.row_decision_bounds.push_back(similarity);
            }

            // TODO: move to decision bound indices to turn some logarithmic-time std::map searches
            // to constant-time
            if (data.sim_value_id_vec.empty()) {
                assert(data.row_decision_bounds.empty());
                assert(data.valid_records_number == 0);
                assert(data.row_lowest == kLowestBound);
                return;
            }
            std::sort(
                    data.sim_value_id_vec.begin(), data.sim_value_id_vec.end(),
                    [](SimValPair const& p1, SimValPair const& p2) { return p1.first > p2.first; });
            std::vector<RecordIdentifier> valid_records;
            valid_records.reserve(data.valid_records_number);
            auto val_rec_begin = valid_records.begin();
            Similarity previous_similarity = data.sim_value_id_vec.begin()->first;
            for (auto [similarity, value_id_right] : data.sim_value_id_vec) {
                data.sim_matrix_row[value_id_right] = similarity;
                auto val_rec_end = valid_records.end();
                if (similarity != previous_similarity) {
                    auto& prev_rec_set = data.matching_recs_mapping[previous_similarity];
                    prev_rec_set.reserve(val_rec_end - val_rec_begin);
                    prev_rec_set.insert(val_rec_begin, val_rec_end);
                    previous_similarity = similarity;
                }
                auto const& records = clusters_right[value_id_right];
                valid_records.insert(val_rec_end, records.begin(), records.end());
            }
            auto& last_rec_set = data.matching_recs_mapping[previous_similarity];
            last_rec_set.reserve(data.valid_records_number);
            last_rec_set.insert(val_rec_begin, valid_records.end());
        }
    };
    thread_pool.ExecIndex(process_value_id, data_left_size);
    thread_pool.WorkUntilComplete();
    for (SimTaskData& task : task_info) {
        similarity_index.push_back(std::move(task.matching_recs_mapping));
        similarity_matrix.push_back(std::move(task.sim_matrix_row));
        decision_bounds_set.insert(task.row_decision_bounds.begin(),
                                   task.row_decision_bounds.end());
        if (task.row_lowest < lowest) lowest = task.row_lowest;
    }
    std::vector<model::md::DecisionBoundary> decision_bounds{decision_bounds_set.begin(),
                                                             decision_bounds_set.end()};
    if (size_limit_ == 0) {
        return {std::move(decision_bounds),
                {lowest, std::move(similarity_matrix), std::move(similarity_index)}};
    }
    return {util::PickMHighestBias(decision_bounds, size_limit_),
            {lowest, std::move(similarity_matrix), std::move(similarity_index)}};
}

LevenshteinSimilarityMeasure::LevenshteinSimilarityMeasure(model::md::DecisionBoundary min_sim,
                                                           bool is_null_equal_null,
                                                           std::size_t size_limit)
    : SimilarityMeasure(std::make_unique<model::StringType>(),
                        std::make_unique<model::DoubleType>()),
      is_null_equal_null_(is_null_equal_null),
      min_sim_(min_sim),
      size_limit_(size_limit) {}

LevenshteinSimilarityMeasure::Creator::Creator(model::md::DecisionBoundary min_sim,
                                               bool is_null_equal_null, std::size_t size_limit)
    : SimilarityMeasureCreator(kName),
      min_sim_(min_sim),
      is_null_equal_null_(is_null_equal_null),
      size_limit_(size_limit) {
    if (!(0.0 <= min_sim_ && min_sim_ <= 1.0))
        throw config::ConfigurationError("Minimum similarity out of range");
}

}  // namespace algos::hymd::preprocessing::similarity_measure
