# CovGuided-ESPRIT

**Covariance-Guided DFT Beam Selection for Beamspace ESPRIT in Hybrid mmWave Sensor Arrays**

> Rıfat Volkan Şenyuva  
> Dept. of Electrical-Electronics Engineering, Maltepe University, Istanbul, Turkey  
> Manuscript ID: Sensors-100615-2025 — IEEE Sensors Journal (under review)  
> Preprint: [arXiv:2512.00898](https://arxiv.org/abs/2512.00898)

---

## Overview

This repository contains the complete MATLAB simulation code for the paper above. The code implements a four-stage covariance-guided pipeline for direction-of-arrival (DoA) estimation in hybrid analog/digital mmWave sensor arrays with a strict DFT beam budget:

1. **Coarse element-space ESPRIT** on a virtual centro-symmetric subarray
2. **NNLS power/noise fitting** to reconstruct the full-aperture covariance
3. **Toeplitz-PSD projection** to obtain a denoised signal covariance
4. **Covariance-guided beam selection + sparse beamspace Unitary ESPRIT** for the fine DoA estimate

The proposed method selects contiguous DFT beam subsets using the denoised covariance rather than raw beamspace estimates, improving DoA accuracy near the detection threshold without increasing the beam budget.

---

## Repository Structure

```
CovGuided-ESPRIT/
├── setup_path.m              ← Run first: adds all subfolders to MATLAB path
├── params/
│   └── getSimParams.m        ← Single source of truth for Table I parameters
├── functions/                ← All 12 core functions (unchanged from submission)
│   ├── tls_esprit.m
│   ├── centered_contiguous_mask.m
│   ├── estimate_powers_course.m   (note: filename typo for 'coarse'; preserved for compatibility)
│   ├── estimate_powers_fine.m
│   ├── toeplitz_projection.m
│   ├── sectorize_dft_beams.m
│   ├── select_adjacent_pairs_from_sectorized_cov.m
│   ├── select_beamsets_from_sectorized_cov.m
│   ├── unitary_esprit_sparse.m
│   ├── compute_subspace_angle_metrics.m
│   ├── plot_pareto_A1.m
│   └── save_vector_pdf.m
├── scripts/                  ← Entry-point scripts (one per figure group)
│   ├── paperCodeMain.m       ← Figs. 1–5 (main paper figures)
│   ├── paperCodeA1.m         ← Fig. S3 + Table S2 (Pareto / timing ablation)
│   ├── paperCodeA2.m         ← Fig. 3 (sector-edge stress test)
│   └── paperCodeS1.m         ← Figs. S4–S5 (fixed-Kf beam budget ablation)
├── data/                     ← Benchmark CSVs committed to the repository
│   ├── A1_timing_table.csv              (8 rows — Table S2 benchmark)
│   ├── S1_aggregates.csv                (24 rows — Figs. S4–S5 benchmark)
│   ├── A2_sector_edge_results_*.csv     (68 rows — Fig. 3 benchmark)
│   └── S1_trials_subset.csv             (120,000 rows — ECDF benchmark; see §Notes)
├── validation/
│   └── validateReproducibility.m   ← Automated RMSE/FailRate check vs. benchmark CSVs
├── results/                  ← gitignored; generated .mat files go here
├── figures/                  ← gitignored; generated .pdf files go here
└── .gitignore
```

---

## Quick Start

```matlab
% 1. Clone the repository
%    git clone https://github.com/<your-username>/CovGuided-ESPRIT.git

% 2. In MATLAB, navigate to the repository root
cd /path/to/CovGuided-ESPRIT

% 3. Add all subfolders to path
setup_path

% 4. Run any entry-point script
paperCodeMain     % generates Figs. 1–5
paperCodeA1       % generates Fig. S3 + Table S2
paperCodeA2       % generates Fig. 3
paperCodeS1       % generates Figs. S4–S5

% 5. Verify numerical consistency against benchmark CSVs
validateReproducibility
```

---

## Reproducing Paper Figures

All scripts use `RandStream('Threefry','Seed',5489)` with `Substream = iter` inside `parfor`. This guarantees that trial `iter` uses the same random draw regardless of which parallel worker processes it, making results fully reproducible across runs and machines.

| Figure(s) | Script | Trials | ASNR range | Approx. runtime* | Benchmark CSV |
|-----------|--------|--------|------------|-------------------|---------------|
| Figs. 1–5 | `paperCodeMain.m` | 10,000 | −10 to +20 dB (31 pts) | ~2–3 h | — |
| Fig. S3 + Table S2 | `paperCodeA1.m` | 10,000 | −4 to +6 dB (11 pts) | ~45 min | `A1_timing_table.csv` |
| Fig. 3 | `paperCodeA2.m` | 10,000 | 3, 6 dB × 17 δ-points | ~45 min | `A2_sector_edge_results_*.csv` |
| Figs. S4–S5 | `paperCodeS1.m` | 10,000 | 0, 5, 10, 15 dB | ~30 min | `S1_aggregates.csv` |

\* Estimated on the development machine (see §Hardware Note below). Runtimes assume a parallel pool with 8 workers.

---

## Validating Reproducibility

Run `validateReproducibility` to check that the installed codebase produces RMSE and failure-rate values consistent with the stored benchmark CSVs. The script runs a focused 10,000-trial test at ASNR ∈ {0, 15} dB for both methods (Cov-guided and Sectorization) with the default beam budget (K_f = 2), then compares against `data/S1_aggregates.csv`.

**Tolerances:**
- RMSE: relative error < 3% (accounts for Monte Carlo variance at 10,000 trials)
- Failure rate: absolute difference < 0.5 percentage points

**Expected output:**
```
Method        ASNR    RMSE_new      RMSE_csv      Err_rel%    FailRate   Status
------------------------------------------------------------------------
Cov-guided    +0 dB   0.XXXX        0.XXXX        X.XX        XX|XX      PASS
Cov-guided   +15 dB   0.0XXX        0.0XXX        X.XX        X|X        PASS
Sectorization +0 dB   0.XXXX        0.XXXX        X.XX        XX|XX      PASS
Sectorization+15 dB   0.XXXX        0.XXXX        X.XX        X|X        PASS

[RESULT] ALL CHECKS PASSED.
```

---

## Hardware Note — Timing Benchmarks

The CSV files committed to this repository were generated on the **original submission machine** (Desktop A, older hardware). The development machine has since been upgraded (newer CPU, motherboard, and RAM). Consequently:

- **RMSE and FailRate columns** in all CSVs are hardware-independent and serve as exact benchmarks. The `validateReproducibility` script checks these.
- **Timing columns** (`t_cov_ms`, `t_sel_ms`, `t_es_ms`, `t_total_ms` in `A1_timing_table.csv`) reflect Desktop A timings. The current machine will produce systematically faster timings. Do not use timing CSVs for regression testing; they are retained only as the paper-submitted values for Table S2.

**Paper Table S2** reports the Desktop A timings. If you re-run `paperCodeA1.m` on different hardware, the RMSE values will match the CSV exactly (within Monte Carlo tolerance), but the timing values will differ.

---

## Dependencies

| Toolbox | Required for | Version tested |
|---------|-------------|----------------|
| MATLAB | Core language | R2022b or later |
| Parallel Computing Toolbox | `parfor` loops in all entry-point scripts | R2022b |
| Optimization Toolbox | `quadprog` in `estimate_powers_course.m`, `estimate_powers_fine.m` | R2022b |

No other toolboxes are required. The code does not use `Statistics and Machine Learning Toolbox`; the `ecdf` call in `paperCodeMain.m` uses MATLAB's built-in `ecdf` function available since R2019b.

---

## Simulation Parameters (Table I)

All shared parameters are defined in `params/getSimParams.m`. Entry-point scripts call `p = getSimParams()` and override only script-specific fields.

| Parameter | Symbol | Value | Notes |
|-----------|--------|-------|-------|
| Array size | M | 32 | ULA, λ/2 spacing |
| RF chains | N_RF | 12 | Virtual subarray size |
| Sources | d | 3 | Uncorrelated, Gaussian |
| Snapshots | N | 100 | Per training phase |
| True DoAs | µ | [−2.1, 0.5, 2.5] rad | Fixed across all scripts |
| Source powers | p | [0.95, 0.50, 0.10] | — |
| ASNR range | — | −10:1:20 dB | Main; others override |
| Sector width | W | 4 beams | A2 overrides to 6 |
| Beam budget | K_f | 2 | Per sector, default |
| Threshold factor | k_thr | 3 | For failure-rate P_fail |
| MC trials | — | 10,000 | Per (ASNR, method) |
| RNG seed | — | 5489 | Threefry; do not change |

---

## Key Implementation Notes

**Typo in `estimate_powers_course.m`:** The filename uses 'course' (a typo for 'coarse'). This is preserved across all files for compatibility with existing results and the submission ZIP. A rename to `estimate_powers_coarse.m` is planned for a future refactoring pass; all call sites will be updated simultaneously.

**`toeplitz_projection.m` vs. Supplementary App. I:** The implementation uses diagonal averaging followed by eigenvalue clipping — a valid fast approximation of the full QP described in Supplementary Appendix I (Eqs. S4–S6). Both enforce the Toeplitz+PSD constraint; the QP's PSD constraint via sampled spectral density inequalities mainly prevents negative eigenvalues, which the clipping achieves directly.

**Two beam selectors:** `select_adjacent_pairs_from_sectorized_cov` is the K=2 specialization with adaptive per-sector conditioning penalty λ. `select_beamsets_from_sectorized_cov` is the general K≥2 selector with fixed α=0.5 used in ablation scripts (A1, S1).

**`paperCodeA2.m` uses d=2 sources:** The sector-edge stress test moves one source across a DFT beam boundary, requiring a simplified 2-source scenario (pows=[0.95,0.50], mu_base=[−2.1;0.5]). All other scripts use d=3 (Table I).

---

## Citation

If you use this code in your research, please cite:

```bibtex
@article{senyuva2025covguided,
  author  = {{\c{S}}enyuva, R{\i}fat Volkan},
  title   = {Covariance-Guided {DFT} Beam Selection for Beamspace {ESPRIT}
             in Hybrid {mmWave} Sensor Arrays},
  journal = {IEEE Sensors Journal},
  year    = {2025},
  note    = {Under review. Preprint: arXiv:2512.00898}
}
```

---

## License

MIT License. See `LICENSE` for details.

---

## Contact

Rıfat Volkan Şenyuva  
Dept. of Electrical-Electronics Engineering, Maltepe University  
34857 Istanbul, Turkey  
rifatvolkansenyuva@maltepe.edu.tr
