# SynVoy: Agent Context
*Strict directives, operational rules, and debugging context for AI agents working on SynVoy.*

## 1. Agent Directives / Rules
- **Environment:** Always run `conda activate synvoy_env` before executing any Python scripts or test modules locally.
- **Execution:** Always use the `-resume` flag when running Nextflow to avoid redundant compute (e.g., `nextflow run main.nf -resume ...`).
- **Paths:** Always use relative paths from the project root. Do not hardcode absolute machine-specific paths (e.g., `/home/faw/...`).
- **Code Modifications:**
  - Python scripts (`bin/*.py`) should be modular and avoid heavy third-party dependencies unless strictly necessary (e.g., use the built-in `sequence_utils.py` instead of BioPython for basic FASTA parsing).
  - When modifying `main.nf`, ensure you fully understand the channel logic. Do not break iterative search parallelization.
- **Testing:** If introducing new logic to core python scripts (like `cluster_grs.py` or sequence handling), verify if tests exist in `tests/` and update them accordingly. Ask the user for permission or guidance when making significant modifications to `main.nf`.

## 2. Project Context
**SynVoy (Synteny Voyager)** is a Nextflow-based pipeline designed to find highly divergent orthologous genes across evolutionary distances using genomic synteny (conservation of gene order) rather than relying solely on pure sequence similarity.

- **The Macro-Syntenic Approach:**
  1. Finds the Gene of Interest (GOI) in a high-quality "home" genome.
  2. Extracts the *n* flanking genes around the GOI.
  3. Searches for those flanking genes in target genomes to identify conserved syntenic blocks.
  4. Runs a rigorous localized search for the GOI *inside* the identified blocks.
- **Key Scripts (`bin/`):**
  - `iterative_search_runner.py`: The heavy-lifting mini-pipeline. Handles complex fallback logic, miniprot annotations, and classifies discoveries via `_classify_goi_evidence`.
  - `cluster_grs.py`: Clusters hits spatially, scoring regions by unique flanking coverage, strand, and order consistency.
  - `plot_synteny.py`: Generates the interactive Plotly HTML visualization.

## 3. Debugging Guide
When a pipeline run crashes, hangs, or produces unexpected results, follow this checklist sequentially:

1. **Check Process Logs:** Nextflow's stdout `.nextflow.log` is good for topology issues, but for process failures, immediately navigate to the `work/<hash>` directory indicated in the error and read `.command.err` and `.command.log`.
2. **Quality Control Prompts / Timeouts:** If `FETCH_HOME_GENOME` or `FETCH_RELATED_GENOMES` hangs or fails with `Home assembly rejected by user/timeout due to low quality`, it is because the parameter `bad_quality_policy` is set to `"ask"`. You may need to bypass this or choose a better assembly.
3. **DSL2 Warnings (Ignore):** The log often throws `WARN: The operator 'first' is useless when applied to a value channel`. Ignore this; it is a known harmless Nextflow DSL2 quirk in this codebase.
4. **Disappearing GOI Candidates:** If a GOI seems missing from the final `synvoy_report.json` or plot, check if `cluster_grs.py` dropped the region. To prevent true true loci from vanishing, it has a fallback mechanism that injects "GOI-anchor" regions.
5. **Plotting Quirks:** `plot_synteny.py` aggressively deduplicates and caps GOI entries (`MAX_GOI_PER_GENOME = 10`) to keep the UI clean. If a specific GOI copy isn't plotted, it might have been filtered out due to low identity relative to other hits at the same locus.
