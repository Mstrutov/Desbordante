#include "algorithms/hyfd/hyfd.h"

#include <algorithm>
#include <chrono>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#include <boost/dynamic_bitset.hpp>
#include <easylogging++.h>

#include "algorithms/hycommon/preprocessor.h"
#include "algorithms/hycommon/util/pli_util.h"
#include "algorithms/hyfd/inductor.h"
#include "algorithms/hyfd/sampler.h"
#include "algorithms/hyfd/validator.h"

namespace algos::hyfd {

HyFD::HyFD(bool request_prepared_data) : PliBasedFDAlgorithm({}, request_prepared_data) {}

unsigned long long HyFD::ExecuteInternal() {
    using namespace hy;
    LOG(TRACE) << "Executing";
    auto const start_time = std::chrono::system_clock::now();

    auto [plis, pli_records, og_mapping] = Preprocess(relation_.get());
    auto const plis_shared = std::make_shared<PLIs>(std::move(plis));
    auto const pli_records_shared = std::make_shared<Rows>(std::move(pli_records));

    Sampler sampler(plis_shared, pli_records_shared);

    auto const positive_cover_tree =
            std::make_shared<fd_tree::FDTree>(GetRelation().GetNumColumns());
    Inductor inductor(positive_cover_tree);
    Validator validator(positive_cover_tree, plis_shared, pli_records_shared);

    IdPairs comparison_suggestions;

    while (true) {
        auto non_fds = sampler.GetNonFDs(comparison_suggestions);

        inductor.UpdateFdTree(std::move(non_fds));

        comparison_suggestions = validator.ValidateAndExtendCandidates();

        if (comparison_suggestions.empty()) {
            break;
        }

        LOG(TRACE) << "Cycle done";
    }

    auto fds = positive_cover_tree->FillFDs();
    RegisterFDs(std::move(fds), og_mapping);

    SetProgress(kTotalProgressPercent);

    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_time);
    return elapsed_milliseconds.count();
}

void HyFD::RegisterFDs(std::vector<RawFD>&& fds, const std::vector<hy::ClusterId>& og_mapping) {
    const auto* const schema = GetRelation().GetSchema();
    for (auto&& [lhs, rhs] : fds) {
        boost::dynamic_bitset<> mapped_lhs =
                hy::RestoreAgreeSet(lhs, og_mapping, schema->GetNumColumns());
        Vertical lhs_v(schema, std::move(mapped_lhs));

        auto const mapped_rhs = og_mapping[rhs];
        // todo(strutovsky): make indices unsigned in core structures
        // NOLINTNEXTLINE(*-narrowing-conversions)
        Column rhs_c(schema, schema->GetColumn(mapped_rhs)->GetName(), mapped_rhs);

        RegisterFd(std::move(lhs_v), std::move(rhs_c));
    }
}

}  // namespace algos::hyfd
