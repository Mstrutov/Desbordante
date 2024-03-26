#include "order.h"

#include "config/names_and_descriptions.h"
#include "config/tabular_data/input_table/option.h"

namespace algos::order {

Order::Order() : Algorithm({}) {
    RegisterOptions();
    MakeOptionsAvailable({config::kTableOpt.GetName()});
}

void Order::RegisterOptions() {
    using namespace config::names;
    using namespace config::descriptions;
    using config::Option;

    RegisterOption(config::kTableOpt(&input_table_));
}

void Order::LoadDataInternal() {
    typed_relation_ = model::ColumnLayoutTypedRelationData::CreateFrom(*input_table_, false);
}

void Order::ResetState() {}

unsigned long long Order::ExecuteInternal() {}

}  // namespace algos::order
