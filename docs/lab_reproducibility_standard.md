# Lab Reproducibility Standard

## Required
- README with purpose/inputs/run/outputs
- config/config.yaml + config/samples.tsv
- metadata/ (dictionary + provenance)
- resources/manifest.tsv (reference versions)
- results/ and logs/ exist but are not committed

## Environments
Use at least one:
- conda env YAMLs (workflow/envs/)
- renv for R (analysis/renv.lock)
- containers if required

## Paths
Avoid hard-coded absolute paths in scripts. Put site-specific paths in config.

## Output policy
Do not commit results/ or logs/. Use releases/archives if sharing outputs.

## Minimal tests
Provide a smoke test in tests/ that verifies basic runability.
