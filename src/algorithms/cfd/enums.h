#pragma once

#include <enum.h>

namespace algos {

BETTER_ENUM(CfdAlgo, char,
            fd_first_dfs_dfs = 0, /* It mines fds with acceptable support threshold
                                   * through dfs lattice traversal
                                   * then it mines patterns for this fds through dfs lattice traversal */
            fd_first_dfs_bfs /* It mines fds with acceptable support threshold
                              * through dfs lattice traversal
                              * then it mines patterns for this fds through bfs lattice traversal */
)

} // namespace algos
