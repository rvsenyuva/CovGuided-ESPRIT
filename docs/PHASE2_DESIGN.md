# Phase 2 Refactoring Design — `runCoarseStage`
## CovGuided-ESPRIT (Sensors-100615-2025)
**Prepared:** 2026-03-31 | **Status:** Design complete — implementation deferred pending review decision

---

## 1. Motivation

The coarse pipeline (Steps 1–6 below) is replicated nearly verbatim inside the
`parfor` body of all four entry-point scripts (`paperCodeMain`, `paperCodeA1`,
`paperCodeA2`, `paperCodeS1`). Phase 2 extracts this block into a single function
`runCoarseStage.m` and moves the shared progress-bar helper `updateProgress` to a
`utils/` subfolder. No algorithm, result, or CSV value changes.

---

## 2. Coarse pipeline steps (shared across all four scripts)

| Step | Operation | Key variables |
|------|-----------|---------------|
| 1 | Synthetic data generation | `S`, `X`, `Z` → `Y_coarse` (M×N) |
| 2 | Subarray slice | `Yb_coarse = Y_coarse(maskIdx,:)` (NRF×N) |
| 3 | Sample covariance + Hermitian symmetrization | `RY_coarse` |
| 4 | Forward-backward averaging (FBA) | `RY_fba` |
| 5 | TLS-ESPRIT | `mu_coarse = tls_esprit(RY_fba, d)` |
| 6 | NNLS power fit + Toeplitz-PSD projection | `[~,~, mu_courseSort, R_course] = estimate_powers_course(...)` |

**Step 7** (`sectorize_dft_beams`) is **outside the scope** of `runCoarseStage`.
The sectorization window width `W` differs across scripts (4 in Main/A1, 6 in A2,
deferred to the kf loop in S1) and the call is one line — keeping it in each
caller is cleaner than parameterizing it inside the function.

---

## 3. Function specification

```matlab
function out = runCoarseStage(cfg)
% =========================================================
% Function: runCoarseStage
% Purpose:  Execute coarse pipeline (Steps 1-6) for one
%           Monte Carlo trial inside parfor.
% Paper:    Sensors-100615-2025
% =========================================================
%
% INPUT  cfg — scalar struct
%
%   Required fields:
%     .Atrue      [M×d]    precomputed steering matrix at mu_true
%     .Psqrt      [d×d]    diag(sqrt(pows/2)); scales source signals
%     .d          scalar   number of sources
%     .N          scalar   snapshot count
%     .M          scalar   array size
%     .nvStd      scalar   noise std = sqrt(noiseVar/2)
%     .maskIdx    [1×NRF]  centered contiguous subarray indices
%     .mIdxFlip   [1×NRF]  NRF:-1:1  (FBA reflection index)
%
%   Optional fields (default: false when absent):
%     .do_baseline          logical  also run wide-ES baseline branch
%                                    Yb_wideES = Y_coarse(1:NRF,:)
%     .return_intermediates logical  also return RY_fba and mu_coarse
%                                    (for A1 external tic/toc timing)
%
% OUTPUT out — scalar struct
%
%   Always present:
%     .mu_courseSort  [d×1]    sorted coarse spatial-freq estimates
%     .R_course       [M×M]    Toeplitz-PSD full-aperture covariance
%
%   Present when do_baseline = true:
%     .mu_wideSort    [d×1]    wide-ES TLS-ESPRIT estimates (sorted)
%     .R_wide         [M×M]    Toeplitz-PSD covariance (wide branch)
%
%   Present when return_intermediates = true:
%     .RY_fba         [NRF×NRF]  FBA-smoothed subarray covariance
%     .mu_coarse      [d×1]      raw TLS-ESPRIT output before NNLS
%     .RY_wideES      [NRF×NRF]  wide-ES sample cov (if do_baseline)
%     .mu_wide_raw    [d×1]      raw wide-ES TLS-ESPRIT (if do_baseline)
```

---

## 4. Per-script divergences and resolutions

| Script | Divergence from shared core | Resolution |
|--------|----------------------------|------------|
| **paperCodeMain** | Runs wide-ES baseline branch | `do_baseline = true` |
| **paperCodeA1** | Times `estimate_powers_course` with `tic`/`toc` for Pareto figure | `return_intermediates = true`; A1 calls step 6 externally with its own timers after receiving `RY_fba` and `mu_coarse` from the function |
| **paperCodeA2** | No baseline; uses local `pad_truncate_to_nrf` at fine stage | `do_baseline = false`; `pad_truncate_to_nrf` stays as a local function in A2 — it operates at the fine stage, outside `runCoarseStage` scope |
| **paperCodeS1** | Runs wide-ES baseline; sectorization deferred to kf loop | `do_baseline = true`; sectorization call remains in the kf loop |

---

## 5. Call-site sketches

### paperCodeMain and paperCodeS1 (inside `parfor`)
```matlab
cfg_c = struct('Atrue', Atrue, 'Psqrt', Psqrt, 'd', d, 'N', N, 'M', M, ...
               'nvStd', nvStd, 'maskIdx', maskIdx, 'mIdxFlip', mIdxFlip, ...
               'do_baseline', true);
out_c = runCoarseStage(cfg_c);
% Outputs used: out_c.mu_courseSort, out_c.R_course,
%               out_c.mu_wideSort,   out_c.R_wide

beamGroups_fine = sectorize_dft_beams(out_c.mu_courseSort, M, W, 1);
```

### paperCodeA2 (inside `parfor`)
```matlab
cfg_c = struct('Atrue', Atrue, 'Psqrt', Psqrt, 'd', d, 'N', N, 'M', M, ...
               'nvStd', nvStd, 'maskIdx', maskIdx, 'mIdxFlip', mIdxFlip, ...
               'do_baseline', false);
out_c = runCoarseStage(cfg_c);
% Outputs used: out_c.mu_courseSort, out_c.R_course

beamGroups = sectorize_dft_beams(out_c.mu_courseSort, M, W, 1);
% pad_truncate_to_nrf remains a local helper below
```

### paperCodeA1 (inside `parfor`, with timing preserved)
```matlab
cfg_c = struct('Atrue', Atrue, 'Psqrt', Psqrt, 'd', d, 'N', N, 'M', M, ...
               'nvStd', nvStd, 'maskIdx', maskIdx, 'mIdxFlip', mIdxFlip, ...
               'do_baseline', false, 'return_intermediates', true);
out_c = runCoarseStage(cfg_c);   % covers Steps 1-5; returns RY_fba, mu_coarse

% Step 6 kept external so A1 can time it independently
t_cov_tic1 = tic;
[~,~, mu_courseSort, R_course] = estimate_powers_course( ...
    out_c.RY_fba, maskIdx, out_c.mu_coarse, M);
t_cov1 = toc(t_cov_tic1);

beamGroups_fine = sectorize_dft_beams(mu_courseSort, M, 4, 1);
```

---

## 6. Inline helper extraction

| Helper | Decision | Rationale |
|--------|----------|-----------|
| `updateProgress` | **Extract** → `utils/updateProgress.m` | Identical in Main, A1, S1 — ~25 duplicated lines; one-file change for any future edits |
| `wilsonCI` | Keep in `paperCodeMain.m` | Single user |
| `prettyAxes` | Keep in `paperCodeMain.m` | Single user |
| `pad_truncate_to_nrf` | Keep in `paperCodeA2.m` | Single user; fine-stage only |
| `pareto_front` | Keep inside `plot_pareto_A1.m` | Already encapsulated in its own function file |

---

## 7. Execution order

Phase 2 implementation proceeds in this order once the trigger condition is met:

1. Create `utils/` subfolder. Move `updateProgress` there. Update the three
   calling scripts (one-line change each). Commit as `phase2/extract-update-progress`.
2. Implement `functions/runCoarseStage.m` per the specification in §3.
3. Update call sites: `paperCodeMain`, `paperCodeA2`, `paperCodeS1` fully;
   `paperCodeA1` with `return_intermediates = true`.
4. Run `validation/validateReproducibility.m`. All four benchmark checks must
   PASS at **0.00% RMSE error** before proceeding.
5. Commit as `phase2/extract-coarse-pipeline`.

**Trigger condition:** If the review decision requires new simulation scripts
(Scripts 1–3 per the paper's anticipated-concerns table), write those scripts
first. They become the primary test cases for `runCoarseStage` and define the
final call-site interface. Implementing `runCoarseStage` before its
revision-script users are known risks a redundant second refactor.

---

## 8. Validation gate

After every commit in Phase 2, re-run:

```matlab
cd validation
validateReproducibility
```

Expected output for all four checks:

```
Check 1 PASS  (Cov-guided,    Kf=2, ASNR=0  dB): RMSE = 0.21484 rad, error = 0.00%
Check 2 PASS  (Cov-guided,    Kf=2, ASNR=15 dB): RMSE = 0.00072 rad, error = 0.00%
Check 3 PASS  (Sectorization, Kf=2, ASNR=0  dB): RMSE = 0.74106 rad, error = 0.00%
Check 4 PASS  (Sectorization, Kf=2, ASNR=15 dB): RMSE = 0.00096 rad, error = 0.00%
```

Any deviation from 0.00% on any check indicates an unintended algorithmic
change and must be resolved before committing.
