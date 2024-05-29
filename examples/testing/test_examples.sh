#!/bin/bash
echo "Testing afd_multiple_error_thresholds" && python3 examples/afd_multiple_error_thresholds.py | diff - examples/testing/outputs/afd_multiple_error_thresholds_output.txt

echo "Testing algebraic_constraints" && python3 examples/algebraic_constraints.py | diff - examples/testing/outputs/algebraic_constraints_output.txt

echo "Testing anomaly_detection" && python3 examples/anomaly_detection.py | diff - examples/testing/outputs/anomaly_detection_output.txt

echo "Testing comparison_pfd_vs_afd" && python3 examples/comparison_pfd_vs_afd.py | diff - examples/testing/outputs/comparison_pfd_vs_afd_output.txt

echo "Testing data_stats" && python3 examples/data_stats.py | diff - examples/testing/outputs/data_stats_output.txt

#command sed "s,\x1B\[[0-9;]*[a-zA-Z],,g" removes ANSI color codes from output
echo "Testing mine_typos" && python3 examples/mine_typos.py | sed "s,\x1B\[[0-9;]*[a-zA-Z],,g" | diff - examples/testing/outputs/mine_typos_output.txt

echo "Testing mining_afd" && python3 examples/mining_afd.py | diff - examples/testing/outputs/mining_afd_output.txt

echo "Testing mining_fd" && python3 examples/mining_fd.py | diff - examples/testing/outputs/mining_fd_output.txt

echo "Testing mining_ind" && python3 examples/mining_ind.py | diff --color=never - examples/testing/outputs/mining_ind_output.txt

echo "Testing mining_list_od" && python3 examples/mining_list_od.py | diff - examples/testing/outputs/mining_list_od_output.txt

echo "Testing mining_pfd" && python3 examples/mining_pfd.py | diff - examples/testing/outputs/mining_pfd_output.txt

echo "Testing mining_set_od_1" && python3 examples/mining_set_od_1.py | diff - examples/testing/outputs/mining_set_od_1_output.txt

echo "Testing mining_set_od_2" && python3 examples/mining_set_od_2.py --color=never | sed "s,\x1B\[[0-9;]*[a-zA-Z],,g" | diff - examples/testing/outputs/mining_set_od_2_output.txt

echo "Testing verifying_aucc" && python3 examples/verifying_aucc.py | diff - examples/testing/outputs/verifying_aucc_output.txt

echo "Testing verifying_fd_afd" && python3 examples/verifying_fd_afd.py | sed "s,\x1B\[[0-9;]*[a-zA-Z],,g" | diff --color=never - examples/testing/outputs/verifying_fd_afd_output.txt

echo "Testing verifying_mfd" && python3 examples/verifying_mfd.py | diff - examples/testing/outputs/verifying_mfd_output.txt

echo "Testing verifying_ucc" && python3 examples/verifying_ucc.py | diff - examples/testing/outputs/verifying_ucc_output.txt

#!!!there are warnings in the current dedupe.py version
echo "Testing dedupe" && python3 -W ignore examples/dedupe.py < examples/testing/inputs/dedupe_input.txt  | diff - examples/testing/outputs/dedupe_output.txt

echo "Testing mining_cfd" && python3 examples/mining_cfd.py | sed "s,\x1B\[[0-9;]*[a-zA-Z],,g" | diff - examples/testing/outputs/mining_cfd_output.txt
