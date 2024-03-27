#pragma once

#include <utility>

#include "model/table/column_index.h"

namespace algos {

// represents a pair of column indices
using ColumnIndexPair = std::pair<model::ColumnIndex, model::ColumnIndex>;

}  // namespace algos
