# Changelog

## v1.2 (Nature Communications revision, 2026-10-23)

Changes made in response to Reviewer 4 and the editor's requests ahead of
Nature Communications publication (NCOMMS-26-022296A):

- **Hardcoded paths removed.** Every absolute, machine-specific path (in
  `rules/qc_extended.smk`, `rules/parent_of_origin.smk`, and elsewhere) is
  now a `config.yml` key with a placeholder default; real values for a given
  machine go in an optional, gitignored `config.local.yml`
  (see `workflow/config.local.yml.example` and `workflow/README.md`).
- **Conda environments pinned.** All resolvable packages in `workflow/env/`
  are now pinned to the exact version used for the published results,
  instead of resolving to whatever is latest at install time.
- **Results, logs, and benchmark files removed from the repository.**
  These are now archived separately with a DOI (see Data Availability /
  Code Availability in the manuscript); `.gitignore` updated accordingly.
- **Repository archived on Zenodo** for a permanent, versioned citation
  (see `CITATION.cff`).
