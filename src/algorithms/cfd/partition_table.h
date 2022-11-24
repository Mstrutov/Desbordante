#pragma once

// see ./LICENSE

#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <algorithm>

#include "model/cfd_types.h"
namespace algos {

class PartitionTable {
public:
    [[maybe_unused]] static int database_row_number_;

    [[maybe_unused]] static PartitionTidList Intersection(const PartitionTidList& lhs,
                                                          const PartitionTidList& rhs);
    [[maybe_unused]] static std::vector<PartitionTidList> Intersection(
            const PartitionTidList& lhs, const std::vector<PartitionTidList*> rhses);
    static std::vector<PartitionTidList> Intersection(
            const PartitionTidList& lhs, const std::vector<const PartitionTidList*> rhses);
    [[maybe_unused]] static std::vector<std::pair<int, std::vector<int> > > PartitionMap(
            const PartitionTidList& x, const PartitionTidList& xa);
    [[maybe_unused]] static std::vector<int> OtherPartitionMap(const PartitionTidList& x,
                                                               const PartitionTidList& xa);
    [[maybe_unused]] static int PartitionError(const std::vector<std::vector<int> >& x,
                                               const std::vector<std::vector<int> >& xa);
    static int PartitionError(const PartitionTidList&, const PartitionTidList&);
    [[maybe_unused]] static int Constify(const PartitionTidList&, const PartitionTidList&);
    [[maybe_unused]] static SimpleTidList Violations(const PartitionTidList&,
                                                     const PartitionTidList&);
    [[maybe_unused]] static SimpleTidList Violations(const PartitionTidList&,
                                                     const PartitionTidList&, int& e);
    [[maybe_unused]] static bool ViolatedInCleaned(const PartitionTidList& x,
                                                   const PartitionTidList& xa, int);
};
} //namespace algos
