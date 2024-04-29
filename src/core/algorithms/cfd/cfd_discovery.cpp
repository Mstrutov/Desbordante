#include "cfd_discovery.h"

#include <iterator>
#include <thread>

#include "algorithms/cfd/util/cfd_output_util.h"
#include "algorithms/cfd/util/set_util.h"
#include "config/equal_nulls/option.h"
#include "config/names_and_descriptions.h"
#include "config/option_using.h"
#include "config/tabular_data/input_table/option.h"

// see algorithms/cfd/LICENSE

namespace algos::cfd {

CFDDiscovery::CFDDiscovery(std::vector<std::string_view> phase_names)
    : Algorithm(std::move(phase_names)) {
    using namespace config::names;
    RegisterOptions();
    MakeOptionsAvailable({kTable, kCfdColumnsNumber, kCfdTuplesNumber});
}

CFDDiscovery::CFDDiscovery() : CFDDiscovery({kDefaultPhaseName}) {}

void CFDDiscovery::LoadDataInternal() {
    relation_ = CFDRelationData::CreateFrom(*input_table_, columns_number_, tuples_number_);

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

    RegisterOption(config::kTableOpt(&input_table_));
    RegisterOption(Option{&columns_number_, kCfdColumnsNumber, kDCfdColumnsNumber, 0ul});
    RegisterOption(Option{&tuples_number_, kCfdTuplesNumber, kDCfdTuplesNumber, 0ul});
}

int CFDDiscovery::NrCfds() const {
    return (int)cfd_list_.size();
}

ItemsetCFDList const& CFDDiscovery::GetItemsetCfds() const {
    return cfd_list_;
}

CFDList CFDDiscovery::GetCfds() const {
    CFDList list;

    auto to_cfd = [&rel = *relation_](ItemsetCFD const& dep) -> RawCFD {
        auto& [lhs, rhs] = dep;

        RawCFD::RawItems cfd_lhs;
        for (uint i = 0; i < lhs.size(); i++) {
            Item item = lhs[i];
            AttributeIndex attr = (item == 0) ? static_cast<AttributeIndex>(i)
                                              : Output::ItemToAttrIndex(item, rel);
            std::optional<std::string> pattern_opt = Output::ItemToPatternOpt(item, rel);
            cfd_lhs.push_back({.attribute = attr, .value = pattern_opt});
        }

        RawCFD::RawItem cfd_rhs = {.attribute = Output::ItemToAttrIndex(rhs, rel),
                                   .value = Output::ItemToPatternOpt(rhs, rel)};
        return RawCFD{cfd_lhs, cfd_rhs};
    };

    for (ItemsetCFD const& dep : GetItemsetCfds()) {
        list.push_back(to_cfd(dep));
    }
    return list;
}

std::string CFDDiscovery::GetCfdString(ItemsetCFD const& cfd) const {
    return Output::CFDToString(cfd, relation_);
}

std::string CFDDiscovery::GetRelationString(char delim) const {
    return relation_->GetStringFormat(delim);
}

std::string CFDDiscovery::GetRelationString(SimpleTIdList const& subset, char delim) const {
    return relation_->GetStringFormat(subset, delim);
}
}  // namespace algos::cfd
