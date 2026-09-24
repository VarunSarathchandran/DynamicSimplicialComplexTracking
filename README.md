# dynamic_sc_paper

Code used for the experiments and figures in
*Tracking Dynamic Simplicial Complexes via Constrained State-Space Estimation*
(ICASSP submission). 

## Reproduce

```matlab
startup_dynamic_sc                            % adds this folder's paths
run_initial_experiment                        % Sec. 4: simulated experiment
run_primary_school_tracking_experiment        % Sec. 5: primary-school experiment
```

Both read their configuration from `configs/` and write to `results/`
(created on first run). The simulated experiment runs 10 realizations of
T = 30 steps; the real-data experiment tracks class 1B (25 students) over
09:20-13:20 on day 2 after calibrating on day 1.

The paper's side-by-side figure is produced with
`dynsc.plot.plotRecoveryMetricsPanels`; see `tools/` below for the export step.


## Requirements

- MATLAB R2024b (developed and run on this release)
- Optimization Toolbox, Statistics and Machine Learning Toolbox
- Gurobi 13 with its MATLAB interface, on the MATLAB path or discoverable via
  `GUROBI_HOME` or `/Library/gurobi*/` (see `startup_dynamic_sc.m`). Used by
  the constrained MAP and constrained RLS estimators and the static baseline.

## Data

`src/+dynsc/+data/` contains the SocioPatterns primary-school contact data
in the AHORN/ScHoLP timestamped maximal-simplex format
(`contact-primary-school.txt`) together with the official SocioPatterns
contact list and metadata used to reconstruct the node relabelling. See
Stehle et al. (2011) and Benson et al. (2018) for the data and format.


