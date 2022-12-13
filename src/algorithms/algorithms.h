#pragma once

/* Functional dependency mining algorithms */
#include "algorithms/aidfd/aid.h"
#include "algorithms/depminer/depminer.h"
#include "algorithms/dfd/dfd.h"
#include "algorithms/fastfds.h"
#include "algorithms/fd_mine.h"
#include "algorithms/fdep/fdep.h"
#include "algorithms/fun.h"
#include "algorithms/hyfd/hyfd.h"
#include "algorithms/pyro.h"
#include "algorithms/statistics/data_stats.h"
#include "algorithms/tane.h"

/* Functional dependency validating algorithms */
#include "algorithms/gfd/gfd_validation.h"

/*Association rule mining algorithms */
#include "algorithms/association_rules/apriori.h"
#include "algorithms/association_rules/borgelt.h"

/* Conditional functional dependency mining algorithms */
#include "algorithms/cfd/cfd_discovery.h"

/* Metric FD verifier */
#include "algorithms/metric/metric_verifier.h"
