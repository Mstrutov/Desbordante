#include "hyucc.h"

#include <chrono>

#include <easylogging++.h>

#include "hycommon/types.h"
#include "inductor.h"
#include "preprocessor.h"
#include "sampler.h"
#include "util/config/thread_number/option.h"
#include "validator.h"

namespace algos {

HyUCC::HyUCC() : UCCAlgorithm({}) {
    RegisterOptions();
}

void HyUCC::RegisterOptions() {
    RegisterInitialExecOption(util::config::ThreadNumberOpt(&threads_num_));
}

void HyUCC::LoadDataInternal() {
    relation_ = ColumnLayoutRelationData::CreateFrom(*input_table_, is_null_equal_null_);

    if (relation_->GetColumnData().empty()) {
        throw std::runtime_error("Got an empty dataset: UCC mining is meaningless.");
    }
}

unsigned long long HyUCC::ExecuteInternal() {
    using namespace hy;
    using namespace hyucc;
    auto const start_time = std::chrono::system_clock::now();

    auto [plis, pli_records, og_mapping] = Preprocess(relation_.get());
    auto const plis_shared = std::make_shared<PLIs>(std::move(plis));
    auto const pli_records_shared = std::make_shared<Rows>(std::move(pli_records));

    hyucc::Sampler sampler(plis_shared, pli_records_shared, threads_num_);

    auto ucc_tree = std::make_unique<UCCTree>(relation_->GetNumColumns());
    Inductor inductor(ucc_tree.get());
    Validator validator(ucc_tree.get(), plis_shared, pli_records_shared, threads_num_);

    IdPairs comparison_suggestions;

    while (true) {
        LOG(DEBUG) << "Sampling...";
        NonUCCList non_uccs = sampler.GetNonUCCs(comparison_suggestions);

        LOG(DEBUG) << "Inducing...";
        inductor.UpdateUCCTree(std::move(non_uccs));

        LOG(DEBUG) << "Validating...";
        comparison_suggestions = validator.ValidateAndExtendCandidates();

        if (comparison_suggestions.empty()) {
            break;
        }
    }

    auto uccs = ucc_tree->FillUCCs();
    RegisterUCCs(std::move(uccs), og_mapping);

    LOG(DEBUG) << "Mined UCCs:";
    for (model::UCC const& ucc : UCCList()) {
        LOG(DEBUG) << ucc.ToString();
    }

    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_time);
    return elapsed_milliseconds.count();
}

void HyUCC::RegisterUCCs(std::vector<boost::dynamic_bitset<>>&& uccs,
                         const std::vector<hy::ClusterId>& og_mapping) {
    const auto* const schema = relation_->GetSchema();
    for (auto&& ucc : uccs) {
        boost::dynamic_bitset<> mapped_ucc =
                hy::RestoreAgreeSet(ucc, og_mapping, schema->GetNumColumns());
        ucc_collection_.Register(schema, std::move(mapped_ucc));
    }
}

}  // namespace algos
