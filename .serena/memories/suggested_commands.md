# Suggested commands
- Python static compile check:
  - `@'
import compileall
print(compileall.compile_dir('src', quiet=1))
'@ | python -`
- Run the main CLI:
  - `paper-pipeline plan --config configs/paper-pipeline.default.yaml all`
  - `paper-pipeline run --config configs/paper-pipeline.default.yaml genebayes --dry-run`
- Generate Slurm scripts:
  - `bash generate_run_once_batch_scripts.sh`
  - `bash generate_trait_batch_scripts.sh`
  - `bash generate_plot_batch_scripts.sh`
  - `bash generate_figure_sbatch.sh`
- Validate shell generators:
  - `bash -n generate_run_once_batch_scripts.sh`
  - `bash -n generate_trait_batch_scripts.sh`
  - `bash -n generate_plot_batch_scripts.sh`
  - `bash -n generate_figure_sbatch.sh`
- Windows-friendly cleanup:
  - `Get-ChildItem -Recurse -Directory -Filter __pycache__ | Remove-Item -Recurse -Force`
- Search:
  - `rg pattern src`
  - `rg --files src`
- Editable install in target envs (Linux examples from project practice):
  - `conda run -n paper-pipeline-control pip install -e . --no-deps`
  - `conda run -n paper-pipeline-plot pip install -e . --no-deps`