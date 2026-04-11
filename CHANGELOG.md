# Changelog

All notable changes to this project are documented in this file.  
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).  
Versioning follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [1.0.0] — 2026-04-17

**Zenodo-archived release.** This version corresponds to the DOI:
`10.5281/zenodo.XXXXXXX` *(placeholder — update after minting Thu 16 Apr)*.

This is the code snapshot accompanying the major revision submission of
"Covariance-Guided DFT Beam Selection for Beamspace ESPRIT in Hybrid mmWave
Sensor Arrays" (Manuscript ID: Sensors-100615-2025, IEEE Sensors Journal,
April 2026).

### Added
- `src/+covguided/` — MATLAB package namespace containing all 12 core and
  utility functions with camelCase names (e.g., `tlsEsprit`, `unitaryEspritSparse`,
  `selectBeamsetsFromSectorizedCov`, etc.)
- `experiments/` — four descriptively-named entry-point scripts:
  `runMainExperiment.m`, `runTimingBenchmark.m`, `runSectorEdgeAblation.m`,
  `runFineStageSweep.m`
- `config/getSimParams.m` — single source of truth for all simulation parameters
- `results/csv/` — committed CSV benchmark results (3 files: A1 timing table,
  S1 aggregates, A2 sector-edge ablation)
- `validation/validateReproducibility.m` — automated benchmark checker; all 4
  checks pass at 0.00% RMSE error against committed CSVs
- `LICENSE` — MIT License
- `CITATION.cff` — GitHub-recognized citation metadata
- `CHANGELOG.md` — this file
- `.gitignore` — excludes `*.asv`, `*.mex*`, `slprj/`, `results/figures/*.pdf`

### Changed
- Repository structure modernized from flat `functions/` + `scripts/` layout to
  MATLAB-package + standard layout (see folder structure in README)
- All function filenames updated from `snake_case_m` convention to `camelCase`
- Typo fixed: `estimate_powers_course` → `estimatePowersCoarse` (the algorithm
  implements coarse-stage estimation; "course" was a spelling error)
- README updated to reflect new structure; Zenodo badge added

### Fixed
- `estimate_powers_course.m`: filename typo corrected to `estimatePowersCoarse.m`
  (function logic unchanged)

---

## [0.1.0] — 2026-03-30

**Phase 1 bootstrap.** Initial public release establishing reproducible codebase.

### Added
- All 16 MATLAB source files organized into `functions/` and `scripts/`
- `params/getSimParams.m` — simulation parameter defaults
- `validation/validateReproducibility.m` — benchmark checker (all 4 checks pass)
- `setup_path.m` — MATLAB path setup helper
- `README.md` — reproduction guide with quickstart, dependency table,
  reproduction table, hardware timing caveat
- `LICENSE` — MIT License
- `.gitignore` — excludes `.mat` result files, figure PDFs, MATLAB auto-saves
- 3 committed CSV benchmark files (`data/` folder)
