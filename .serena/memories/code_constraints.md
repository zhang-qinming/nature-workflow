# Code modification constraints

- `code.backup/` — Original paper code backup. **NEVER modify** anything in this directory. It serves as a reference only.
- `src/paper_pipeline/dataprep/` — Data processing code adapted from the original. Safe to edit but **do not alter the main pipeline logic**.
- `src/paper_pipeline/plot/` — Plotting code adapted from the original. Safe to edit but **do not alter the main pipeline logic**.
- `src/paper_pipeline/figures/` — New reusable figure workflows combining original code patterns with generated data. **Can be modified for better generalization** — this is the preferred layer for improvements.

**Why:** User explicitly stated these boundaries to protect the paper reproduction fidelity while allowing the newer `figures` layer to evolve.
**How to apply:** When making changes, treat `code.backup` as read-only. For `dataprep` and `plot`, make only minimal, compatibility-preserving edits. Direct refactoring and generalization work to `figures`.
