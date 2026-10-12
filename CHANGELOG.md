# Changelog

All notable changes to this project are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project uses [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.3.0] - 2026-10-12

### Changed

- **Breaking:** default pKa values now follow Thurlkill et al. 2006
  ([doi:10.1110/ps.051840806](https://doi.org/10.1110/ps.051840806)), matching
  Faunus. Previous values were uncited free-amino-acid textbook numbers.
  Changes: ASP 3.65 → 3.67, HIS 6.00 → 6.54, CYS 8.18 → 8.55,
  TYR 10.07 → 9.84, LYS 10.54 → 10.40, ARG 12.48 → 12.0, CTR 3.1 → 3.67
  (GLU 4.25 and NTR 8.0 unchanged). Net charges shift accordingly, mostly
  through His near neutral pH; e.g. BSA (PDB 4F5S) at pH 7 in the ideal limit
  goes from −15.5 to −12.7 e. Pass the old values with `--pka` to reproduce
  earlier results.
- **Breaking (library):** `NTERM_GROUP.res_name` and `CTERM_GROUP.res_name` are
  now `"NTR"` and `"CTR"` (were `"NTERM"` and `"CTERM"`), matching the bead
  types and `--pka` names.
- PQR files get one `REMARK pKa NAME VALUE` line per group, and scan files a
  `# pKa ...` line after the column header, recording the pKa set used.

### Added

- `--pka NAME=PKA` option, repeatable, to override the pKa of any titratable
  group (ASP, GLU, HIS, CYS, TYR, LYS, ARG, NTR, CTR). Names are
  case-insensitive; unknown names, non-finite values and repeats are rejected.
- Library: `residue::PkaOverrides` (validated on every route, including serde),
  the `ChargeCalc::pka_overrides` field and the `ChargeCalc::with_pka` builder.
  Overrides apply to both Henderson-Hasselbalch and MC.

## [0.2.4] - 2026-10-12

### Added

- `coarse_grain_jointly`: coarse-grain several molecules into one shared Faunus
  topology, so pair potentials cover cross-molecule bead pairs. Each molecule is
  still titrated on its own.

[Unreleased]: https://github.com/mlund/cgkitten/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/mlund/cgkitten/compare/v0.2.4...v0.3.0
[0.2.4]: https://github.com/mlund/cgkitten/releases/tag/v0.2.4
