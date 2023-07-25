#include "algorithms/cfd/cfd_discovery.h"

#include <iterator>
#include <thread>

#include "algorithms/cfd/util/cfd_output_util.h"
#include "algorithms/cfd/util/set_util.h"
#include "util/config/equal_nulls/option.h"
#include "util/config/names_and_descriptions.h"
#include "util/config/option_using.h"
#include "util/config/tabular_data/input_table/option.h"

// see algorithms/cfd/LICENSE

namespace algos::cfd {

CFDDiscovery::CFDDiscovery(std::vector<std::string_view> phase_names)
    : Algorithm(std::move(phase_names)) {
    RegisterOptions();
}

CFDDiscovery::CFDDiscovery() : CFDDiscovery({kDefaultPhaseName}) {}

void CFDDiscovery::LoadDataInternal() {
    relation_ = CFDRelationData::CreateFrom(*input_table_, is_null_equal_null_, columns_number_,
                                            tuples_number_);

    if (relation_->GetColumnData().empty()) {
        throw std::runtime_error("Got an empty .csv file: CFD mining is meaningless.");
    }
}

void CFDDiscovery::ResetState() {
    cfd_list_.clear();
    ResetStateCFD();
}

void CFDDiscovery::RegisterOptions() {
    DESBORDANTE_OPTION_USING;

    RegisterInitialLoadOption(util::config::TableOpt(&input_table_));
    RegisterInitialLoadOption(util::config::EqualNullsOpt(&is_null_equal_null_));
    RegisterInitialLoadOption(Option{&columns_number_, kCfdColumnsNumber, kDCfdColumnsNumber, 0u});
    RegisterInitialLoadOption(Option{&tuples_number_, kCfdTuplesNumber, kDCfdTuplesNumber, 0u});
}

int CFDDiscovery::NrCfds() const {
    return (int)cfd_list_.size();
}

CFDList CFDDiscovery::GetCfds() const {
    return cfd_list_;
}

std::string CFDDiscovery::GetCfdString(CFD const& cfd) const {
    return Output::CFDToString(cfd, relation_);
}

std::string CFDDiscovery::GetRelationString(char delim) const {
    return relation_->GetStringFormat(delim);
}

std::string CFDDiscovery::GetRelationString(const SimpleTIdList& subset, char delim) const {
    return relation_->GetStringFormat(subset, delim);
}
}  // namespace algos::cfd
