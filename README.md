# CovGuided-ESPRIT

**Covariance-Guided DFT Beam Selection for Beamspace ESPRIT in Hybrid mmWave Sensor Arrays**

[![arXiv](https://img.shields.io/badge/arXiv-2512.00898-b31b1b.svg)](https://arxiv.org/abs/2512.00898)
[![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.19553631.svg)](https://doi.org/10.5281/zenodo.19553631)

**IEEE Sensors Journal (major revision submitted, April 2026)**
Manuscript ID: Sensors-100615-2025

---

## Overview

This repository contains the complete MATLAB simulation code for the paper:

> R. V. Şenyuva, "Covariance-Guided DFT Beam Selection for Beamspace ESPRIT
> in Hybrid mmWave Sensor Arrays," *IEEE Sensors Journal*, under review, 2026.
> Preprint: [arXiv:2512.00898](https://arxiv.org/abs/2512.00898)

The code implements a four-stage coarse-to-fine direction-of-arrival (DoA)
estimation pipeline for hybrid analog/digital mmWave sensor arrays:

1. **Coarse TLS-ESPRIT** on a virtual centro-symmetric element-space subarray
2. **NNLS covariance fitting** with Toeplitz-PSD projection
3. **Covariance-guided DFT beam selection** (the core contribution)
4. **Fine Unitary ESPRIT** on the selected contiguous beam subset

---

## Repository structure (v1.0)

```
CovGuided-ESPRIT/
├── src/
│   └── +covguided/             % MATLAB package — all core functions
│       ├── tlsEsprit.m
│       ├── unitaryEspritSparse.m
│       ├── centeredContiguousMask.m
│       ├── estimatePowersCoarse.m
│       ├── estimatePowersFine.m
│       ├── toeplitzProjection.m
│       ├── sectorizeDftBeams.m
│       ├── selectAdjacentPairsFromSectorizedCov.m
│       ├── selectBeamsetsFromSectorizedCov.m
│       ├── computeSubspaceAngleMetrics.m
│       ├── plotParetoA1.m
│       ├── saveVectorPdf.m
│       ├── dftDictionary.m
│       ├── crbSpatialFreq.m
│       └── drawSpatialFreqs.m
├── experiments/                % Entry-point scripts (run these)
│   ├── runMainExperiment.m          % Main paper figures (Figs. 1–3)
│   ├── runTimingBenchmark.m         % Pareto/timing ablation (Fig. S3, Table S2)
│   ├── runSectorEdgeAblation.m      % Sector-edge stress test (Fig. S4)
│   ├── runFineStageSweep.m          % Fine-stage beam budget sweep (Figs. S4–S5)
│   └── runPhaseQuantizationStudy.m  % Phase-shifter quantization robustness (Fig. 4)
├── config/
│   └── getSimParams.m          % Single source of truth for all Table I parameters
├── results/
│   └── csv/                    % Committed benchmark CSVs
│       ├── A1_timing_table.csv
│       ├── S1_aggregates.csv
│       ├── A2_sector_edge_results_M32_W6_20251214_153453.csv
│       └── phase_quant_study_20260412_192029.csv
├── validation/
│   └── validateReproducibility.m  % Automated benchmark checker
├── docs/
│   ├── CITATION.cff
│   └── CHANGELOG.md
├── LICENSE                     % MIT
├── .gitignore
└── README.md                   % This file
```

---

## System requirements

- MATLAB R2020b or later (Parallel Computing Toolbox required for `parfor`)
- No additional toolboxes required beyond PCT
- Tested on MATLAB R2023b (Linux) and R2024a (Windows)

---

## Quickstart

```matlab
% 1. Add the package and config to path (run once per session)
addpath(genpath(fullfile(pwd, 'src')));
addpath(fullfile(pwd, 'config'));

% 2. Run the main experiment (reproduces Figs. 1–3 of the paper)
cd experiments
runMainExperiment

% 3. Verify reproducibility against committed CSV benchmarks
cd ..
validation/validateReproducibility   % all 4 checks should PASS at 0.00% error
```

---

## Reproducing benchmark results

The validation script `validation/validateReproducibility.m` runs 10⁴
Monte Carlo trials at ASNR ∈ {0, 15} dB using Threefry seed 5489 with
per-trial substreams, and checks four RMSE values against
`results/csv/S1_aggregates.csv`:

| Check | Method | ASNR | RMSE [rad] | FailRate | sqrtCRB [rad] |
|---|---|---|---|---|---|
| 1 | Cov-guided (kf=2) | 0 dB | 0.21484 | 0.1571 | 0.004153 |
| 2 | Cov-guided (kf=2) | 15 dB | 0.00072382 | 0.0048 | 0.000630 |
| 3 | Sectorization (kf=2) | 0 dB | 0.74106 | 0.4445 | 0.004153 |
| 4 | Sectorization (kf=2) | 15 dB | 0.00096286 | 0.0181 | 0.000630 |

Expected output: all rows `PASS` at 0.00% error.

**Hardware timing caveat:** Timing columns in `A1_timing_table.csv` are
machine-specific. The validation script checks RMSE and FailRate only, never timing.

---

## Experiment scripts

All scripts are in `experiments/` and use repo-relative paths — they can be
run from any working directory without manual `addpath` calls.

| Script | Output | Approx. runtime (6-core parpool) |
|---|---|---|
| `runMainExperiment.m` | Figs. 1–3 PDFs in `results/figures/` | ~25 min |
| `runTimingBenchmark.m` | `A1_timing_table.csv`, Pareto PDF | ~15 min |
| `runSectorEdgeAblation.m` | `A2_sector_edge_results_*.csv`, PDFs | ~5 min |
| `runFineStageSweep.m` | `S1_aggregates.csv`, `S1_trials_subset.csv`, PDFs | ~5 min |
| `runPhaseQuantizationStudy.m` | `phase_quant_study_*.csv`, Fig. 4 PDF | ~3 min |

---

## Function reference

All functions live in the `+covguided` package under `src/`. Call them as
`covguided.functionName(...)` after adding `src/` to the MATLAB path.

| Function | Signature | Purpose |
|---|---|---|
| `tlsEsprit` | `(R, d) → mu_est` | TLS-ESPRIT coarse DoA estimator |
| `unitaryEspritSparse` | `(Y, beamCols, d, M) → mu_UE` | Fine-stage Unitary ESPRIT |
| `centeredContiguousMask` | `(M, NRF) → idx` | Virtual subarray index set |
| `estimatePowersCoarse` | `(R_hat, maskIdx, mu_hat, M) → [p, s2, mu_sort, R_est]` | Element-space NNLS fit |
| `estimatePowersFine` | `(R_hat, beamCols, mu_hat, M) → [p, s2, mu_sort, R_est]` | Beamspace NNLS fit |
| `toeplitzProjection` | `(R) → T_psd` | Toeplitz-PSD projection |
| `sectorizeDftBeams` | `(mu, M, W, zeta) → beam_groups` | DFT sectorization |
| `selectAdjacentPairsFromSectorizedCov` | `(R_model, beam_groups, opts) → pairs` | K=2 beam selector |
| `selectBeamsetsFromSectorizedCov` | `(R_model, beam_groups, K, opts) → sets` | General K-beam selector |
| `computeSubspaceAngleMetrics` | `(U_true, U_est, d) → theta` | LPA subspace metric |
| `dftDictionary` | `(M) → F` | Centered DFT dictionary (FFT phase convention) |
| `crbSpatialFreq` | `(R_signal, noiseVar, N, M) → sqrtCRB` | CRB approximation helper |
| `drawSpatialFreqs` | `(mu_true, d) → mu` | Random spatial frequency draw |

---

## Citation

If you use this code, please cite:

```bibtex
@article{senyuva2026covguided,
  author  = {Şenyuva, Rıfat Volkan},
  title   = {Covariance-Guided {DFT} Beam Selection for Beamspace {ESPRIT}
             in Hybrid mmWave Sensor Arrays},
  journal = {IEEE Sensors Journal},
  year    = {2026},
  note    = {Major revision submitted, April 2026. Preprint: arXiv:2512.00898}
}

@misc{senyuva2026covguided_code,
  author    = {Şenyuva, Rıfat Volkan},
  title     = {{CovGuided-ESPRIT}: Simulation code for covariance-guided
               {DFT} beam selection for beamspace {ESPRIT}},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.19553631},
  url       = {https://doi.org/10.5281/zenodo.19553631}
}
```

---

## License

MIT License — see [LICENSE](LICENSE).  
Copyright © 2026 Rıfat Volkan Şenyuva, Maltepe University, Istanbul, Turkey.
