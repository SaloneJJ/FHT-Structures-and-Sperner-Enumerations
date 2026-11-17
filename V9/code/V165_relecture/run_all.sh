#!/bin/bash
export OPENBLAS_NUM_THREADS=1
cd /home/claude/fht/relecture2
timeout 1500 python3 rv_flow.py wil flat 3 -1 48 123 G2 > rv_flow_wil_flat_N3.log 2>&1
timeout 1200 python3 rv_flow.py wil tree 3 -1 32 123 G2 > rv_flow_wil_tree_N3_M2.log 2>&1
timeout 1200 python3 rv_spectres.py BC > rv_spectres_BC.log 2>&1
timeout 1200 python3 rv_jauge.py 3 p4 16 > rv_jauge_p4_N3_16.log 2>&1
timeout 1500 python3 rv_jauge.py 3 p4 64 > rv_jauge_p4_N3_64.log 2>&1
timeout 1200 python3 rv_flow.py sgO flat 3 1 48 13 G2 > rv_flow_sgO_flat_N3.log 2>&1
timeout 1200 python3 rv_jauge.py 3 flot 96 b > rv_jauge_flot_b_N3_96.log 2>&1
timeout 1200 python3 rv_flow.py wil tree 3 -1 128 13 G2 > rv_flow_wil_tree_N3_128.log 2>&1
timeout 1500 python3 rv_flow.py wil tree 4 -1 48 13 G2 > rv_flow_wil_tree_N4.log 2>&1
echo FINI > run_all.done
