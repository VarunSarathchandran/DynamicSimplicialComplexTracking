# dynamic_sc_paper

Exactly the code used for the experiments and figures in
*Tracking Dynamic Simplicial Complexes via Constrained State-Space Estimation*
(ICASSP submission). This is a pruned copy of the `dynamic_sc` development
repository, produced by static dependency analysis
(`matlab.codetools.requiredFilesAndProducts`) from the two experiment entry
points and the paper-figure plotting function. Nothing else was carried over.

## Reproduce

```matlab
startup_dynamic_sc                            % adds this folder's paths
run_initial_experiment                        % Sec. 4: simulated experiment
run_primary_school_tracking_experiment        % Sec. 5: primary-school experiment
```

Both read their configuration from `configs/` and write to `results/`
(created on first run). The simulated experiment runs 10 realizations of
T = 30 steps; the real-data experiment tracks class 1B (25 students) over
09:20-13:20 on day 2 after calibrating on day 1. The real run takes roughly
2.5 h on a laptop because the calibration day is filtered as well.

The paper's side-by-side figure is produced with
`dynsc.plot.plotRecoveryMetricsPanels`; see `tools/` below for the export step.

## Configs match the paper runs

`configs/default_config.m` and `configs/primary_school_tracking_config.m`
are set to the values stored in the `config.mat` saved alongside the reported
results, not to the development repository's current defaults. The fields
that were reset for the simulated experiment are: `experiment.name`,
`experiment.numRealizations = 10`, `signals.numNodeSignals = 1`,
`signals.numEdgeSignals = 1`, `baseline.identityTransition.enabled = true`,
`plot.stateSnapshots.interval = 2`, `output.saveRealizations = true`. For the
real-data experiment only `experiment.name` and
`validation.compareRawPairTimeEvents = false` were set. Note in particular
that the simulated experiment uses one node signal and one edge signal
(M0 = M1 = 1), while the real-data experiment uses fifty of each; the paper
states both values.

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

## tools/

- `exportEmbeddedPdf.m`: exports a figure as a vector PDF with embedded,
  subsetted fonts at an exact physical size, which `print -dpdf` cannot do.
  It compensates for a 0.5 device scale that `exportgraphics` applies on
  Retina Macs and refuses (with a warning) if that scale is not present, so
  check the warning on other machines.
- `checkfonts.py`: `python3 tools/checkfonts.py file.pdf` reports embedding
  and subsetting of every font, the check IEEE PDF eXpress performs.

## Deliberately excluded

Results and figures; the graph-only (edge-birth-death) experiment and its
`+graph` package; the regime-change experiment; the SC-construction and
state-snapshot sanity scripts; the test suite; the legacy `opt/` and `utils/`
solvers apart from the six generation helpers in `utils/Generation` that the
signal model calls; development notes and backup files.
