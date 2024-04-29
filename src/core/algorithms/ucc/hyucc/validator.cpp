#include "validator.h"

#include <future>

#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>

#include "fd/hycommon/efficiency_threshold.h"
#include "fd/hycommon/validator_helpers.h"
#include "ucc/hyucc/model/ucc_tree_vertex.h"
#include "util/bitset_utils.h"

namespace {

using model::RawUCC;

// Check out comment for Validator::ValidateAndExtendCandidates, it's appropriate for this function
// too
size_t AddExtendedCandidatesFromInvalid(std::vector<algos::hyucc::LhsPair>& next_level,
                                        algos::hyucc::UCCTree& ucc_tree,
                                        std::vector<RawUCC> const& invalid_uccs,
                                        size_t num_attributes) {
    size_t candidates = 0;
    for (auto const& ucc : invalid_uccs) {
        for (size_t attr = 0; attr < num_attributes; ++attr) {
            if (ucc.test(attr)) {
                continue;
            }

            RawUCC ucc_ext = ucc;
            ucc_ext.set(attr);

            if (ucc_tree.FindUCCOrGeneralization(ucc_ext)) {
                continue;
            }

            algos::hyucc::UCCTreeVertex* child = ucc_tree.AddUCCGetIfNew(ucc_ext);
            if (child == nullptr) {
                continue;
            }
            next_level.emplace_back(child, std::move(ucc_ext));
            candidates++;
        }
    }
    return candidates;
}

}  // namespace

namespace algos::hyucc {

using model::RawUCC;

bool Validator::IsUnique(model::PLI const& pivot_pli, RawUCC const& ucc,
                         hy::IdPairs& comparison_suggestions) {
    std::vector<hy::ClusterId> indices = util::BitsetToIndices<hy::ClusterId>(ucc);
    for (model::PLI::Cluster const& cluster : pivot_pli.GetIndex()) {
        auto cluster_to_record =
                hy::MakeClusterIdentifierToTMap<model::PLI::Cluster::value_type>(cluster.size());
        for (auto const record_id : cluster) {
            std::vector<hy::ClusterId> cluster_id =
                    hy::BuildClustersIdentifier((*compressed_records_)[record_id], indices);
            if (cluster_id.empty()) {
                continue;
            }

            auto const [it, inserted] =
                    cluster_to_record.try_emplace(std::move(cluster_id), record_id);
            if (!inserted) {
                comparison_suggestions.emplace_back(record_id, it->second);
                return false;
            }
        }
    }

    return true;
}

Validator::UCCValidations Validator::GetValidations(LhsPair const& vertex_and_ucc) {
    auto [vertex, ucc] = vertex_and_ucc;
    UCCValidations validations;
    validations.SetCountValidations(1);
    validations.SetCountIntersections(1);

    size_t ucc_attr = ucc.find_first();
    bool is_unique;
    if (current_level_number_ == 1) {
        assert(ucc_attr != boost::dynamic_bitset<>::npos);
        is_unique = (*plis_)[ucc_attr]->AllValuesAreUnique();
    } else {
        ucc.reset(ucc_attr);
        // It's guaranteed that the first cluster has the fewest records because of sorting
        // in the Sampler
        model::PLI const* pivot_pli = (*plis_)[ucc_attr];
        is_unique = IsUnique(*pivot_pli, ucc, validations.ComparisonSuggestions());
        ucc.set(ucc_attr);
    }

    if (!is_unique) {
        vertex->SetIsUCC(false);
        validations.InvalidInstances().push_back(std::move(ucc));
    }

    return validations;
}

Validator::UCCValidations Validator::ValidateAndExtendSeq(
        std::vector<LhsPair> const& current_level) {
    UCCValidations result;
    for (auto const& vertex_and_ucc : current_level) {
        if (!vertex_and_ucc.first->IsUCC()) {
            continue;
        }
        result.Add(GetValidations(vertex_and_ucc));
    }
    return result;
}

Validator::UCCValidations Validator::ValidateAndExtendParallel(
        std::vector<LhsPair> const& current_level) {
    UCCValidations result;
    boost::asio::thread_pool pool(threads_num_);
    std::vector<std::future<UCCValidations>> validation_futures;

    for (auto const& vertex_and_ucc : current_level) {
        if (!vertex_and_ucc.first->IsUCC()) {
            continue;
        }

        std::packaged_task<UCCValidations()> task(
                [this, &vertex_and_ucc]() { return GetValidations(vertex_and_ucc); });
        validation_futures.push_back(task.get_future());
        boost::asio::post(pool, std::move(task));
    }

    pool.join();

    for (auto& future : validation_futures) {
        assert(future.valid());
        result.Add(future.get());
    }

    return result;
}

Validator::UCCValidations Validator::ValidateAndExtend(std::vector<LhsPair> const& current_level) {
    assert(threads_num_ > 0);
    if (threads_num_ > 1) {
        return ValidateAndExtendParallel(current_level);
    } else {
        return ValidateAndExtendSeq(current_level);
    }
}

// This function is similar to hyfd::Validator::ValidateAndExtendCandidates(), but with some key
// differences in the way it initializes the current_level and handles loop conditions. While there
// is some shared structure, attempting to refactor the functions to eliminate duplication
// introduced excessive complexity and reduced readability. As a result, it was decided to keep the
// two functions separate and accept the minor duplication to maintain simplicity and clarity.
hy::IdPairs Validator::ValidateAndExtendCandidates() {
    size_t const num_attributes = plis_->size();
    assert(num_attributes == tree_->GetNumAttributes());

    size_t previous_num_invalid_uccs = 0;
    std::vector<LhsPair> current_level = tree_->GetLevel(current_level_number_);
    hy::IdPairs comparison_suggestions;
    while (!current_level.empty()) {
        UCCValidations result = ValidateAndExtend(current_level);
        comparison_suggestions.insert(comparison_suggestions.end(),
                                      result.ComparisonSuggestions().begin(),
                                      result.ComparisonSuggestions().end());
        std::vector<LhsPair> next_level = hy::CollectCurrentChildren(current_level, num_attributes);

        size_t candidates = AddExtendedCandidatesFromInvalid(
                next_level, *tree_, result.InvalidInstances(), num_attributes);

        LogLevel(current_level, result, candidates, current_level_number_, "UCC");

        size_t const num_invalid_uccs = result.InvalidInstances().size();
        size_t const num_valid_uccs = result.CountValidations() - num_invalid_uccs;
        current_level = std::move(next_level);
        current_level_number_++;

        if (num_invalid_uccs > (long double)hy::kEfficiencyThreshold * num_valid_uccs &&
            previous_num_invalid_uccs < num_invalid_uccs) {
            return comparison_suggestions;
        }
        previous_num_invalid_uccs = num_invalid_uccs;
    }

    return {};
}

}  // namespace algos::hyucc
