# Done checklist
- After code changes, run Python static compile check on `src`.
- If shell generators changed, run `bash -n` on each modified generator.
- For new `figures` workflows or scripts, validate task building with a minimal `LoadedConfig` if full YAML loading is unavailable.
- Clean `__pycache__` before finishing if `compileall` was used.
- In this Windows environment, Linux/Slurm/R runtime validation is limited; call that out explicitly in final notes.