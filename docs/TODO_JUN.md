# SynVoy — June 2026 TODO (melittin deep-dive follow-up)

**Created:** 2026-06-09. Source: a full empirical investigation of the melittin runs
(`local_runs/melittin_full` = 19 Hymenoptera, home *Apis mellifera*, query melittin
P01501; `local_runs/melittin_reg` = 5 bees) against NCBI / on-disk ground truth, plus an
independent subagent verification of every claim. Companion to `docs/TODO.md` (this file is
the actionable next-step list; §-numbers below map back to `docs/TODO.md` where relevant).

> **Historical working document — dated 2026-06-09, not maintained.** Per-item ✅ markers
> below record what shipped. As of 2026-07-26 the outstanding items are:
> **A4** (rank-based wave binning) — *implemented but opt-in*, `--rank_wave_binning true`,
> default off pending validation; **A5** (eggNOG OG-gating) — *not implemented*, no eggNOG
> code exists in the repo. Everything else marked ✅ is in the codebase.
>
> **Note for readers on GitHub:** `docs/TODO.md`, referenced throughout, is a local
> working file excluded by `.gitignore` and is not published. The §-numbers are internal
> cross-references; the user-facing behaviour they describe is documented in
> [USAGE.md](USAGE.md) and [PARAMETERS.md](PARAMETERS.md).

---

## HUMAN INPUT

check why and what the different evidences for gois are (tandem, miniprot, fallback) check critically if this is even sound and makes sense, there eis high doubt about that, check in wirth testing data (melettins, etc)

## Diagnosis (verified)

The synteny engine is **sound** — it recovers the conserved melittin flanking neighbourhood
(GR2) across all 19 Hymenoptera and lands the *Tetramorium* melittin candidate on the
paper's curated scaffold **OV788322**. Do **not** "fix" the synteny. The failures are all
downstream of block placement:

1. **Misses are a coverage-beats-bitscore inversion (the "fumble", `docs/TODO.md` §1.2/A2).**
   The real *Xylocopa* xylopin is at `scaffold-133:~10,190,330` (signal `MKFLFSVLVLFMIALVSLSNA`
   + acidic pro `EPEPDPEA`), in the syntenic gap between flankers LOC726827@10.189 Mb and
   LOC726866@10.197 Mb. The search **finds it** (mmseqs@9.5 50%/e≈1e-3, tblastn 48%, SW) in the
   correct flanking-seeded block **b0** — but b0's fallback was **rejected at `qcov=0.18`
   (`valid_fallback=False`)** while the spurious block **b1** (`scaffold-133:370,534`, 28%, not
   melittin) was **accepted at `qcov=0.40`**. b0 had the higher bitscore (91 vs 45) and lost on
   coverage. Cause: the 26-aa mature peptide is ~18% of the 70-aa query, so qcov is structurally
   low even on a perfect hit. (miniprot additionally cannot model the short 2-exon melittin even
   at `--outs=0.2 --outc=0.1`; the SW/tblastn fallback in the correct block is the intended path.)

2. **The wave iterative design exists but is neutered.** Genomes are binned into
   phylo-distance waves (`iterative_search_runner.py:~5498`), parallel within a wave
   (`ProcessPoolExecutor ~5522`), `latest_db` grows across waves. But
   `_SEED_QUALIFYING_CLASSES = {confident_goi, probable_goi}` (`~4040`) withholds
   `tandem_goi_copy`/fallback/ambiguous GOI from the expanding query DB, so divergent melittin
   orthologs never seed the next wave — even **Apis cerana's 100%-identical melittin seeded
   nothing** (it was classed `tandem_goi_copy`). 9/19 genomes contributed 0 GOI-derived seeds.
   The wave distance metric is also **saturated** (all 4 waves ≈0.99), so there is no real
   close→far gradient.

3. **Overcalls are fabricated local-window identities.** *Vollenhovia* "79.4% melittin, 4 copies"
   has a real melittin identity of only **~23%** (tblastn, e≈0.6), no melittin signature, and the
   copies aren't similar to each other. `tandem_copy` features inherit `best_hit['pident']`
   (a single local-window identity, `annotate_goi_exons.py:940`) and are written with
   **`query_cov=None`** (`iterative_search_runner.py:3077`), so the short-window inflation is
   invisible. *Polistes* "71%" is a 7-aa micro-window amplified by assembly fragmentation
   (calls scattered across 4 scaffolds).

4. **Reporting bug.** `.fa` is missing from `generate_report.py:KNOWN_SUFFIXES`. The hull-rescue
   module names its (often empty) GFF from the bare stem (`Colletes_gigas.hull_rescue.gff` →
   canonical `Colletes_gigas`), while the real GFF is `Colletes_gigas.fa.gff` → canonical
   `Colletes_gigas.fa`. So each `.fa` Pro-mode genome splits into a real record + an empty phantom
   that lands in `summary.goi_absent_genomes` — falsely reporting **all 5 solitary/stingless bees
   as "melittin absent"** when they were recovered. `.fna` genomes escape (`.fna` is stripped).

### Hard constraints (must hold for every item)
- **No biology-specific gates** — SynVoy is a general tool. Use only generic signals (true
  coverage, RBH / eggNOG OG consistency, the existing §1m ownership / §1j self-consistency).
- **Do not regress LY6 (single-seed) or TP53 (paralog)** with any fallback rework.
- **Do not kill the *Tetramorium* melittin ortholog** (OV788322, 40% identity × **76% coverage**).
  No gate keys on the identity×coverage product — only coverage floors that 76% clears.
- **No new `main.nf` process calls** except the eggNOG item (JVM method-size limit, §1m).
- User commits on `main`; this is a plan, not commits.

---

## Quick wins (low-risk, hours) — ✅ DONE 2026-06-09 (all 3 shipped, 323 tests pass)

### QW1 — `.fa`/`.fasta` → `KNOWN_SUFFIXES`  ✅ map: NEW (near §1h)
- **File:** `bin/generate_report.py`, `KNOWN_SUFFIXES` (~line 15).
- **Change:** add `".fasta"` and `".fa"` to the suffix list, longest-first ordering so
  `.fna`/`.faa`/`.fasta` match before `.fa`. Makes `canonical_genome_id()` collapse the empty
  hull-rescue phantom into the real genome id.
- **Risk:** very low. **Acceptance:** `canonical_genome_id("Apis_mellifera.fa.gff") ==
  "Apis_mellifera"`; re-aggregated `melittin_full` report has empty `goi_absent_genomes`;
  unit case added to `tests/`.

### QW2 — record true `query_cov` on tandem-copy GOI features  ✅ map: §1.1 / §1j
- **File:** `bin/iterative_search_runner.py`, tandem emit block (~3056–3079); copy dict from
  `bin/annotate_goi_exons.py:detect_tandem_duplications` (~933–950) already carries `qstart/qend`.
- **Change:** compute `query_cov = (qend - qstart + 1) / len(parent_query_seq)` per copy and pass
  it into the GFF attrs instead of `query_cov=None`. Purely additive metadata; does not change
  confidence by itself (that's QW3 / item A2-tandem).
- **Risk:** low. **Acceptance:** Vollenhovia's tandem "79%" feature carries a low `QueryCoverage`
  (~0.3, not None); Tetramorium OV788322 (not a tandem feature) unchanged.

### QW3 — `identity_coverage_decoupled` self-consistency flag  ✅ map: §1j / §1.1
- **File:** `bin/generate_report.py`, `build_self_consistency` (~line 914).
- **Change:** add a generic advisory flag for any GOI call with `identity ≥ high_min_identity`
  but missing/low `QueryCoverage` (`< 0.35`). Report-only, no demotion.
- **Risk:** low. **Acceptance:** Vollenhovia/Polistes appear as `identity_coverage_decoupled`;
  Tetramorium (76% coverage) does not.

---

## ⚠️ A0 — REPRODUCIBILITY (NEW, top priority, found 2026-06-09 rerun)

**The iterative search is non-deterministic run-to-run, severely.** Re-running the exact
melittin_full benchmark (`local_runs/melittin_full_qw`, only QW1–3 metadata edits vs the
Jun-4 `local_runs/melittin_full`) produced wholesale-different GOI calls:

| | Jun-4 run | Jun-9 rerun |
|---|---|---|
| Tetramorium OV788322 (paper's key locus) | **40.1% MEDIUM** | **26.3% LOW** |
| Vollenhovia | 79.4% MEDIUM | 24.1% LOW |
| Colletes | 58.8% MEDIUM | 33.3% LOW |
| headline | 2 HIGH + 20 MED | 2 HIGH + 9 MED |
| GOI-derived seed genomes | 10 | 7 (only **4** overlap) |

**Proven NOT caused by QW1–3:** `git diff bin/iterative_search_runner.py` is the QW2
metadata-only change (adds `QueryCoverage`; cannot touch pident or the tandem/fallback
decision); the **flanking gene identities are byte-identical** between runs (LOC726817=90.9
%, LOC100578330=62.5%, …) — only the divergent GOI calls moved.

**Mechanism:** mmseqs runs `--threads 2` (non-deterministic prefilter/alignment on marginal
hits) → a divergent genome's GOI flips confident/probable↔ambiguous between runs → that
flips whether it **seeds** the expanding `latest_db` (`_SEED_QUALIFYING_CLASSES`) → which
changes every later wave's query DB → a cascade. The seed sets diverged: only Apis_florea,
Tetragonula, Solenopsis, Cardiocondyla seeded in *both* runs. So the OLD run's higher GOI
identities (the Vollenhovia "79%", etc.) were local-window matches to *accumulated bee-
melittin seeds* that were present in one run's expanded DB and absent in the other's — which
also reframes the §C overcall finding: the inflation is seed-feedback + local pident, not a
fixed bug.

**Why it's blocking:** (1) publication risk — the paper's Tetramorium melittin (OV788322)
swings MEDIUM↔LOW between identical-input runs; a reviewer re-running gets different numbers.
(2) It **confounds all A1–A5 validation** — you can't measure a fix against a baseline that
moves by itself.

**CONFIRMED 2026-06-09 by a two-rep test (same code, same params, twice):** **14/19 genomes
give different GOI calls.** Headline swung **"2 HIGH + 9 MEDIUM" (rep1) vs "2 HIGH + 20 MEDIUM"
(rep2)**. Two attractor states: rep2's 10-genome seed set == the Jun-4 `melittin_full`; rep1's
7-genome set == `melittin_full_qw`. So the "MEDIUM 20→9" earlier was just attractor choice, NOT
the QW edits (fully exonerated). Flanking is stable in 15/19 (the perturbation is in marginal/
divergent hits; the wave-seed cascade amplifies it).

**Two pinned software sources:**
1. **Non-deterministically ordered query DB** — `iterative_search_runner.py:5613`
   `write_fasta([... for g in wave_results], ...)` writes each wave's new genes in
   `as_completed` order (line 5532), so `expanded_db.faa` (next wave's mmseqs query) has a
   run-dependent order → marginal-hit tie-breaks flip → seed set flips → cascade. **Cheap fix:
   sort `wave_results` deterministically (by genome then id) before `write_fasta`.**
2. **mmseqs marginal jitter** — `--threads 2` + `--mmseqs_split_memory_limit 8G` (memory-split +
   multithread) perturbs borderline hits (flanking jitters in 4/19). **Fix:** deterministic
   split (`--split 1` / fixed `--split-mode`) and/or `--threads 1` for the GOI search; measure
   the speed cost.

**A0 tasks:**
- [ ] **Decisive measurement:** run the *same* (current) code twice, diff GOI identities +
  seed sets, to quantify the noise floor. (~47 min/run.) **← the remaining step: rerun
  `local_runs/melittin_full` twice on the patched code and confirm GOI identities + seed
  sets now match. Blocked only on compute time, not code.**
- [x] **Make the GOI search deterministic.** ✅ SHIPPED 2026-06-10 (see below). (a) augmented
  GOI mmseqs forced to `--threads 1` via `--deterministic_goi_search` (default on); tblastn was
  already single-threaded. (b) `parse_hits` + the hit-merge / candidate-selection sorts now
  carry a total tie-break order. (c) the wave seed set still keys on confident/probable identity
  — the *stable-signal* seed gate is A3 (deferred); A0 instead removes the upstream jitter that
  was flipping that signal.
- [x] **`-resume` determinism guard** — the per-wave expansion DB is now assembled in the wave's
  fixed genome order (not `as_completed`), so the wave/DB byte-order is stable across runs and
  across `-resume`. Full A4 wave-order recalibration still deferred.

### ✅ A0 software fixes — SHIPPED 2026-06-10 (`bin/iterative_search_runner.py`, all 328+4 tests pass)
Three always-/default-on changes, no `main.nf` process additions:
1. **Pinned source #1 — deterministic wave assembly.** The wave loop now collects each genome's
   new genes / tree-extras into per-name dicts and, after the `ProcessPoolExecutor` drains,
   assembles `wave_results` by walking the wave in its fixed genome order and sorting each
   genome's genes by id. So `expanded_db.faa` (the next wave's mmseqs query DB) is byte-identical
   run-to-run. Live progress logging still prints in `as_completed` order. (Replaces the old
   `wave_results.extend(...)` in completion order at ~line 5573/5613.)
2. **Task 2b — deterministic hit ordering.** `parse_hits` imposes a total order
   `(query, chrom, start, end, strand, -bits)` on its output, so a run-dependent .m8 row order
   (mmseqs `--threads>1`) can no longer flip downstream *stable*-sort tie-breaks. The merged-hit
   greedy-dedup sort (`-bits, evalue, …`) and the GOI/flanking candidate-selection sorts gained
   explicit `(chrom, start, end, query)` / `(start, record.id)` tie-breakers.
3. **Pinned source #2 — `--deterministic_goi_search` (default `true`).** Forces the per-region
   augmented GOI mmseqs search to `--threads 1`. The region is tiny, so the speed cost is small,
   but it removes mmseqs' multithread marginal-hit jitter — the root of the GOI confidence flips.
   Wired through `nextflow.config` + `modules/iterative_search.nf`. Set `false` for the legacy
   fast (non-deterministic) path. **Note:** the *genome-wide* flanking easy-search still runs
   multithreaded (flanking was stable in 15/19); narrowing that residual is left for follow-up
   under the same flag if the A0 two-run diff still shows flanking drift.

Regression: `tests/test_a0_determinism.py` (parse_hits row-order invariance + flag default).
**Caveat:** the new flag changes the ITERATIVE_SEARCH task hash, so the next run will not
`-resume` from a pre-patch work dir (expected — A0 validation needs fresh runs anyway).

This re-prioritises: **A0 first**, then A3+A4 (seeding/wave determinism are the same root),
then A1 (fumble), then A2, then A5.

## Validation (cheap, before the architectural work)
- Re-run the melittin case under `-profile preset_short_peptide` (lower classify/qcov thresholds)
  and check whether Xylocopa/Colletes true loci are recovered — tests whether the qcov gate is the
  whole story.
- Check whether the **frozen paper copy** of SynVoy predates the hull-rescue module; if not, the
  paper's `goi_absent` counts are wrong for any `.fa` target.

---

## Architectural work (days)

### A1 — Fix the qcov fallback inversion = "the fumble"  ⭐ map: §1.2/A2 — ✅ CODE SHIPPED 2026-06-10
- **File:** `bin/iterative_search_runner.py`, per-block hit→model→emit path (`process_single_genome`
  block loop ~2895–3132; `valid_fallback` gate ~2962–3051; greedy non-overlapping selection ~3754).
- **Change:** (a) rank GOI fallback candidates by **flanking-supported bitscore**, not raw qcov, so
  the real-but-divergent b0 model beats a high-coverage/low-bitscore other-block artifact; (b)
  compute qcov against the **mature/aligned core** (or relax the short-query qcov floor) so a 26-aa
  hit on a 70-aa query isn't structurally rejected; (c) ensure the fallback is emitted **inside the
  block that carries the GOI hits**, not a different block. Stay in-script (no `main.nf` changes).
- **Impact:** highest science value (helps single-copy divergent genes too). **Risk:** med-high
  (candidate-selection core) — gate the ranking change on GOI-role candidates only; reuse §1m
  collinearity ranking. **Acceptance:** Xylocopa recovered at `scaffold-133:10.19 Mb` ≥ MEDIUM, the
  0.37 Mb/28% artifact no longer the sole GOI; Tetramorium OV788322 unchanged; melittin_reg 5 GT
  bees still on the paper scaffold; LY6/TP53 HIGH counts unchanged.

**✅ Implementation (2026-06-10, `bin/iterative_search_runner.py`, 9 new tests, 337 pass):**
Root cause confirmed by reading the code: `b0` and `b1` are **different blocks**. `b0`'s block
produced no GOI candidate at all because the `valid_fallback` gate vetoed the single 26-aa hit on
`qcov<0.25` (the 70-aa melittin preprotein makes a perfect mature-peptide hit structurally
low-coverage); `b1`'s block passed at qcov 0.40 → the GOI was emitted in the wrong block. Fix
unifies (a)+(b)+(c) into the gate, now a testable pure helper `is_valid_fallback(...)`:
- **(b)+(a)** For queries shorter than `--fallback_short_query_len` (150) the fallback gates on
  **absolute aligned length** (`--fallback_short_min_aln_aa`, 15) **+ bitscore**
  (`--fallback_short_min_bits`, 30) — coverage is no longer a veto. So `b0` (aln 26, bits 91)
  passes; a 7-aa micro-window or a weak (bits<30) hit is still rejected. Long queries keep the
  qcov floor (LY6/TP53 untouched); Tetramorium (76% cov) clears either path. `=0` ⇒ legacy gate.
- **(c)** Because `b0`'s block now emits, the GOI lands in its own flanking-seeded block; `b1`'s
  28% call still emits but is classified LOW by identity, so it is "no longer the sole GOI".
- **(a) ranking** The fallback candidate now carries `bits`/`flanking_support`, and the per-block
  GOI greedy selection breaks identity/length ties by flanking-supported bitscore.
- Wired `fallback_short_*` through `nextflow.config` + `modules/iterative_search.nf` (read as
  `params.X`, like the §1m bridge params — no `settings`/preset wiring needed; defaults already fix
  the general short-peptide case). Tests: `tests/test_a1_fallback_gate.py`.
- **⏳ Remaining (compute, not code):** a `local_runs/melittin_full` rerun to confirm the live
  acceptance — Xylocopa recovered at `scaffold-133:10.19 Mb` ≥ MEDIUM, Tetramorium OV788322
  unchanged, melittin_reg 5 GT bees still on the paper scaffold, LY6/TP53 HIGH counts unchanged.
  Run **on top of A0** so the before/after diff is against a deterministic baseline.

### A2 — Coverage-aware tandem classification  map: §1.1 — ✅ CODE SHIPPED 2026-06-10
- **File:** `bin/iterative_search_runner.py`, tandem branch of `_classify_goi_evidence` (~894–897);
  `CLASSIFY_THRESHOLDS` (~163–164).
- **Change:** MEDIUM only when `identity ≥ tandem_min_identity AND query_cov ≥ tandem_min_qcov`
  (new ~0.35); else LOW. **Depends on QW2** supplying real coverage.
- **Risk:** med — demoted tandems must still flow to `tree_extra` for the phylogeny. **Acceptance:**
  Vollenhovia/Polistes demoted/flagged; Tetramorium stays; LY6/TP53 unchanged.

**✅ Implementation (2026-06-10, 5 new tests):** the `tandem_copy` branch now requires BOTH
`identity ≥ tandem_min_identity` AND `qcov ≥ tandem_min_qcov` (new `CLASSIFY_THRESHOLDS` key,
default 0.35, CLI `--classify_tandem_min_qcov`, wired via `params.X` in
`modules/iterative_search.nf` — NOT `settings`, to avoid editing the JVM-method-size-limited
`workflow{}` body in `main.nf`). A high-identity short local window (Vollenhovia "79%"/Polistes
"71%" — now carrying their true low qcov from QW2) is demoted to LOW with a distinct reason
`goi_tandem_copy_low_coverage` (vs `_low_identity`), so QW3's `identity_coverage_decoupled` flag
and A2's demotion are complementary. Tetramorium is an exon_annotation/fallback, not a tandem, so
it's untouched; LY6/TP53 aren't tandems either. Tests: `tests/test_a2_a3_seed_tandem.py`.
- **Deliberate deviation from the risk note:** an A2-demoted tandem becomes LOW and therefore does
  NOT flow to `tree_extra` (existing design: `test_low_confidence_dropped_from_both` — "LOW feeds
  neither seed nor tree, too noisy"). I kept that: the A2-demoted cases are precisely the fabricated
  overcalls (Vollenhovia is ~23% real identity, not melittin), so keeping them in the published
  phylogeny would be wrong. If the rerun shows a *genuine* low-coverage ortholog dropped from the
  tree, add a carve-out routing `goi_tandem_copy_low_coverage` (identity-passed) tandems to
  `tree_extra` — but that would also re-admit the overcalls, so validate first.
- **⏳ Remaining (compute):** confirm on the melittin_full rerun that Vollenhovia/Polistes drop to
  LOW and Tetramorium/LY6/TP53 are unchanged.

### A3 — Un-neuter the wave seed  ⭐ map: §1.2/A2 (wave machinery) — ✅ CODE SHIPPED 2026-06-10 (opt-in)
- **File:** `bin/iterative_search_runner.py`, `_SEED_QUALIFYING_CLASSES` (~4040),
  `_classify_goi_for_seed_and_tree` (~4044–4078).
- **Change:** admit a GOI-role feature to `expanded_db.faa` on a **generic** signal — strong
  flanking support (`flanking_support ≥ high_min_flanking`) AND real coverage (`query_cov ≥` floor)
  — regardless of `goi_class`, so divergent/tandem orthologs bridge to the next wave. Behind a kill
  switch (`--seed_on_flanking_support`); default decided after the rerun.
- **Risk:** med (paralog drift) — coverage+flanking floor, not identity; §1m collinearity limits
  drift; A4 (eggNOG) is the backstop. **Acceptance:** Apis cerana melittin seeds `expanded_db`; a
  previously-missed divergent bee gains ≥ MEDIUM; Tetramorium unchanged; TP53 keeps the intended
  paralog.

**✅ Implementation (2026-06-10, 7 new tests, default OFF):** `_classify_goi_for_seed_and_tree`
gained keyword args `seed_on_flanking_support` / `seed_min_flanking` / `seed_min_qcov` (defaults
preserve the old 3-arg behaviour, so the existing tandem-tree tests are untouched). When the flag is
ON, a **HIGH/MEDIUM** GOI feature ALSO seeds — regardless of `goi_class` — when
`block_flanking_support ≥ seed_flanking_min_count` (2) AND `query_coverage ≥ seed_flanking_min_qcov`
(0.5), parsed from the GFF attrs via a safe `_safe_float` (empty `QueryCoverage` ⇒ can't seed).
Keys on coverage+flanking, not identity. **Confidence floor is HIGH/MEDIUM by design:** a strong-
flanking *fallback* is already MEDIUM `probable_goi` and seeds via the normal path, so the only thing
A3 newly admits is the MEDIUM `tandem_goi_copy` (Apis cerana) — LOW noise is still excluded. Wired
`seed_on_flanking_support`/`seed_flanking_min_count`/`seed_flanking_min_qcov` via `params.X` in
`modules/iterative_search.nf` (no `workflow{}`/settings edits). CLI `--seed_on_flanking_support`.
Tests: `tests/test_a2_a3_seed_tandem.py`.
- **Default OFF** per "default decided after the rerun" — flipping it on without measuring risks
  paralog drift (TP53). **A0 is the precondition:** A3's whole mechanism is the wave-seed cascade,
  which was non-deterministic before A0 — only now can the rerun attribute a change to A3.
- **⏳ Remaining (compute):** rerun melittin_full with `--seed_on_flanking_support true` (vs the
  default-off baseline) and confirm Apis cerana seeds, a divergent bee gains ≥ MEDIUM, and
  TP53/Tetramorium are unchanged; if clean, flip the default to `true`.

### A4 — Recalibrate the saturated wave-distance metric  map: §1.3 / §1.2 (NEW sub-item)
- **File:** `bin/iterative_search_runner.py`, distance normalization + wave binning (~5385–5434);
  upstream `compute_tree.py` / phylo_sort.
- **Change:** diagnose the ≈0.99 saturation (likely `max_dist` normalization ~5388–5394) and switch
  to a rank/quantile wave assignment so close→far actually grades (precondition for A3's closest-first
  seeding).
- **Risk:** med — changes execution order → `-resume` hashing; behind a flag, verify determinism.
  **Acceptance:** wave log shows graded `dist≈` (≥2 tiers); combined with A3 a missed bee is recovered
  via a closer relative's seed.

### A5 — eggNOG OG-gating (generic family-consistency safety net)  map: §1.4/§2a
- **File:** NEW `bin/` step on the reconstructed candidate GOI protein; demotion in
  `generate_report.py` (alongside §1m `paralog_not_goi`, ~1101); behind `--eggnog_classify` (default
  off). The one item that adds a `main.nf` module — reuse the lean wiring pattern (§1m caveat).
- **Change:** compute candidate OG vs query OG; mismatch ⇒ relabel as the real paralog + demote.
  Biology-agnostic backstop for A1/A3 drift on **any** family.
- **Risk:** med (DB/runtime footprint) — opt-in until benchmarked. **Acceptance:** with
  `--eggnog_classify true`, true melittin orthologs (incl. OV788322) keep the query's OG and stay
  HIGH/MEDIUM; a wrong-paralog overcall (or TP53→TP63/TP73 mis-seed) gets demoted; default-off path
  byte-identical to today.

---

## Suggested sequencing
QW1 → QW2 → QW3 (independently shippable, make every later rerun legible) → cheap validation →
**A1 (the fumble)** → A2 → A3 + A4 together → A5 as the generic safety net. Items A1–A5 each need a
`local_runs/melittin_full` rerun; LY6/TP53 reruns are blocked on refetching their target genomes
(`docs/TODO.md` §1m).
