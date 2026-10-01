@echo off
setlocal
python run_external_adaptive.py freeze --plan configs/protocol_external_adaptive.json --out results/external_adaptive --reason "Prospective external-comparator tranche: H1/H5, delta 0.30, all four high-precision tolerances, 15 paired repetitions with a predeclared mechanical extension rule to 31."
endlocal
