# Repository Guidelines

## Project Structure & Module Organization
- Core code lives in `python_project/`: `src/gravite_py/` (1D/2D/3D lattice physics utilities), `scripts/` (repro runs), `outputs/` (JSON/PNG results), `reports/` (research notes), and `requirements.txt`.
- Historical research logs are in the root (`IA1`, `IA2 v2`, `IA3`, `IA4`, `IA5`, `GEMINI.md`, `mega_rapport_IA1_IA5.md`). Treat them as read-only; append new findings instead of editing.
- Large datasets/archives (e.g., `Enzo_64.tar.gz`, `Documents drive/`) stay at the root. Do not move or overwrite; reference paths explicitly when used.
- Place new experiments under `python_project/outputs/` and summarize in `python_project/reports/` with a short Markdown note.

## Build, Test, and Development Commands
- Use a virtual env inside `python_project/`:
  - `cd python_project && python3 -m venv .venv && source .venv/bin/activate`
  - `pip install -r requirements.txt`
- Key runs (keep sizes small when iterating):
  - Forces (1D/2D/3D): `python3 scripts/run_qit_scaling.py --n 128 --onsite 0 --defect-masses 2 4 8 16`
  - Geometry (MI profile): `python3 scripts/run_qit_geometry.py --n 64 --onsite 1e-6 --defect-mass 20`
  - Higher-dim MI: `python3 scripts/run_2d_gravity.py` and `python3 scripts/run_3d_gravity.py`
  - Dynamics: `python3 scripts/run_dynamics_lensing.py` (deflection) and `python3 scripts/run_orbit.py` (capture)
  - Vector field experiment: `python3 scripts/run_vector_scaling.py`
- Artifacts land in `python_project/outputs/`; keep names descriptive (e.g., `qit_scaling_N128.json`, `gravity_well_heatmap.png`).

## Coding Style & Naming Conventions
- Follow PEP 8, 4-space indent, snake_case for functions/vars, PascalCase for classes; add type hints and succinct docstrings.
- Keep code comments brief and purposeful; default to English in code, French acceptable in reports.
- Use deterministic seeds where possible; avoid non-ASCII in code unless required by existing files.

## Testing Guidelines
- Prefer fast sanity runs before heavy jobs: small `--n` (e.g., 64) and modest defect masses.
- Verify MI positivity and stable fits: inspect generated JSON/PNG (log-log lines for power laws, heatmaps for geometry).
- Re-run a known baseline (e.g., `run_qit_geometry.py` with default params) after major changes to confirm behavior.
- Record parameters and seeds in the accompanying report note; do not delete prior outputs.

## Commit & Pull Request Guidelines
- Use conventional commit prefixes (`feat:`, `fix:`, `docs:`, `data:`, `refactor:`).
- Scope commits narrowly; include command snippets or parameter sets in the description when adding outputs.
- Never overwrite historical reports; add new files with clear names (date or phase), and reference predecessors.
- Large artifacts should be compressed and placed in `outputs/` with approximate sizes noted in the report.

## Security & Configuration Tips
- No secrets are expected; avoid embedding credentials in scripts or configs.
- Keep dependency changes confined to `requirements.txt`; document version bumps in commit messages.
- When handling external datasets (e.g., Enzo), document source and checksums in a short note under `reports/` or `outputs/`.
