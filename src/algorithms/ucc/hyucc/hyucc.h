#pragma once

#include <memory>

#include "hycommon/types.h"
#include "model/column_layout_relation_data.h"
#include "model/idataset_stream.h"
#include "ucc/ucc_algorithm.h"
#include "util/config/thread_number/type.h"

namespace algos {

class HyUCC : public UCCAlgorithm {
private:
    std::unique_ptr<ColumnLayoutRelationData> relation_;
    util::config::ThreadNumType threads_num_ = 1;

    void RegisterOptions();
    void LoadDataInternal() override;
    unsigned long long ExecuteInternal() override;
    void ResetUCCAlgorithmState() override {}
    void RegisterUCCs(std::vector<boost::dynamic_bitset<>>&& uccs,
                      const std::vector<hy::ClusterId>& og_mapping);

public:
    HyUCC();
};

}  // namespace algos
