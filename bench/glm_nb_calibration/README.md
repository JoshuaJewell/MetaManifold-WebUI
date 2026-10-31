<!-- SPDX-License-Identifier: CC-BY-SA-4.0 -->
# Reproduce the small-sample glm_nb calibration result

This is the simulation published in [upstream issue #64](https://github.com/JoshuaJewell/MetaManifold-WebUI/issues/64).
The simulation and fitting logic match the issue script. The licence header and per-method summary grouping have been updated. It runs the actual
`_FIT_R` block from the supplied `src/analysis/differential.jl`, rather than
reimplementing the fitter. This branch is based on upstream main at
`a325589d3d318847c15b21cd34fb7d6e7a6d2f7b`, which includes PR #62.
The differential-analysis source is unchanged from that PR.

## Run

From the repository root, using R with MASS installed:

```sh
Rscript bench/glm_nb_calibration/calibrate-glm-nb.R
```

The script writes `calibration-replicates.csv` in the current directory and
prints the R/MASS versions, source MD5, summaries and session information.
To keep generated output outside the checkout, run from a separate directory
and pass the absolute source path as the first argument:

```sh
mkdir -p /tmp/metamanifold-calibration
cd /tmp/metamanifold-calibration
Rscript /path/to/MetaManifold-WebUI/bench/glm_nb_calibration/calibrate-glm-nb.R \
  /path/to/MetaManifold-WebUI/src/analysis/differential.jl > calibration-output.txt
```

Defaults are ten datasets of 600 taxa for each of 4 and 12 samples per group.
The fixed seed for each dataset is recorded in the CSV. `CALIBRATION_REPS`
and `CALIBRATION_TAXA` override the defaults; changing either is a different
experiment, not a reproduction of the published table.

`COMPARE_GLMGAMPOI=1` additionally requires exactly glmGamPoi 1.22.0 and runs
its optional comparison. It fails if that package/version is unavailable.

## Scope

The generator has no true group effect. The reported rate counts raw p-values
below 0.05 among valid tests. Failed fits are counted separately; boundary
fits retain their p-values, matching the production block. The standard error
is calculated across independent datasets, not by treating taxa sharing
estimated TSS offsets as independent observations.

This experiment does not test BH-adjusted false-discovery control, exercise
the Julia/HTTP wrapper, establish behaviour on every dataset, or select a new
default test. The original issue used R 4.3.3 and MASS 7.3-60.0.1, rather than
the project's pinned R 4.5 environment.

## Reproduction on 30 October 2026

The full default experiment completed successfully on Ubuntu 24.04.3,
R 4.3.3 and MASS 7.3-60.0.1. Source revision:
`a325589d3d318847c15b21cd34fb7d6e7a6d2f7b`; source-file MD5:
`fd8c875dfd50fe699e556d37d1670b6c`.

All counts match the table in issue #64:

| Samples per group | Retained | Valid | Failed | Boundary | Rejected | Pooled rejection rate |
|---|---:|---:|---:|---:|---:|---:|
| 4 | 6,000 | 5,874 | 126 | 61 | 782 | 13.3129% |
| 12 | 6,000 | 5,999 | 1 | 0 | 445 | 7.4179% |

The mean replicate rates were 13.3171% and 7.4179%, with replicate
standard errors of 0.5881 and 0.3471 percentage points, respectively.
The output contained 20 replicate rows, whose aggregate counts were checked
against the issue's table. This is a rerun of the same seeds and generator,
not an independent simulation design. Generated CSV and runtime downloads
are intentionally kept outside the source tree.

The run used R's `--vanilla --slave -f` entry point with the source path
passed after `--args`; it executes the same script as the Rscript command
above. BLAS and LAPACK were Ubuntu's system libraries (3.12.0).
The optional glmGamPoi comparison, pinned R 4.5 environment and Julia/HTTP
integration were not run. No production analysis or dependency pin changes
are included.
