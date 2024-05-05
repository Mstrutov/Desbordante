#pragma once

#include <memory>
#include <utility>
#include <vector>

#include "algorithms/fd/hycommon/primitive_validations.h"
#include "algorithms/fd/hyfd/model/fd_tree.h"
#include "algorithms/fd/raw_fd.h"
#include "config/thread_number/type.h"
#include "model/table/position_list_index.h"
#include "types.h"

namespace algos::hyfd {

using LhsPair = fd_tree::LhsPair;

class Validator {
private:
    using FDValidations = hy::PrimitiveValidations<RawFD>;

    std::shared_ptr<fd_tree::FDTree> fds_;

    hy::PLIsPtr plis_;
    hy::RowsPtr compressed_records_;

    unsigned current_level_number_ = 0;

    FDValidations ProcessZeroLevel(LhsPair const& lhsPair);
    FDValidations ProcessFirstLevel(LhsPair const& lhs_pair);
    FDValidations ProcessHigherLevel(LhsPair const& lhs_pair);

    FDValidations GetValidations(LhsPair const& lhsPair);

    FDValidations ValidateAndExtendSeq(std::vector<LhsPair> const& vertices);

    FDValidations ValidateAndExtendPar(std::vector<LhsPair> const& vertices);

    [[nodiscard]] unsigned GetLevelNum() const {
        return current_level_number_;
    }

    config::ThreadNumType threads_num_ = 1;

public:
    Validator(std::shared_ptr<fd_tree::FDTree> fds, hy::PLIsPtr plis,
              hy::RowsPtr compressed_records, config::ThreadNumType threads_num) noexcept
        : fds_(std::move(fds)),
          plis_(std::move(plis)),
          compressed_records_(std::move(compressed_records)),
          threads_num_(threads_num) {}

    hy::IdPairs ValidateAndExtendCandidates();
};

}  // namespace algos::hyfd
