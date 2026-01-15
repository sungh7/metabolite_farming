#!/bin/bash
# Run all docking jobs

vina --config /data/ethylene/results/docking/novel_candidates/config_8E83_6-O-Malonylgenistin_.txt
vina --config /data/ethylene/results/docking/novel_candidates/config_8E83_6-O-Acetyldaidzin_CI.txt
vina --config /data/ethylene/results/docking/novel_candidates/config_8E83_6-O-Malonyldaidzin_C.txt
vina --config /data/ethylene/results/docking/novel_candidates/config_1EYQ_Liquiritigenin_CID11.txt
vina --config /data/ethylene/results/docking/novel_candidates/config_1EYQ_Daidzein_CID5281708.txt
vina --config /data/ethylene/results/docking/novel_candidates/config_8EA1_Daidzin_CID107971.txt
vina --config /data/ethylene/results/docking/novel_candidates/config_8EA1_Genistin_CID5281377.txt
