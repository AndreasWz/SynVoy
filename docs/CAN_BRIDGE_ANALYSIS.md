# `_can_bridge()` — code analysis, figure alignment, and whether to adapt it

Prompted by the observation that fig3's two "refused" cases both look bridgeable.
Everything below was produced by **executing** `identify_synteny_blocks` /
`_can_bridge` / `_longest_collinear_run`, on synthetic cases and on two real runs.
Nothing is quoted from a summary.

Date: 2026-09-09. Code state: working tree at `dev` + uncommitted F9/AMBIGUOUS work.
(That work was committed on 2026-10-05, `ab31c7a` and `96f1d53`. Two-sided bridging is in the
code as an opt-in, `--synteny_bridge_two_sided`, default off.)

---

## 1. What the function actually does

`_can_bridge(block, candidate, home_rank, max_rank_gap, min_anchors)`
(`bin/iterative_search_runner.py:1819`) is called from the **sequential left-to-right
walk** in `identify_synteny_blocks` Step 3. Three properties matter and none of them
are visible on the figure:

1. **It is evaluated per candidate LOCUS, not per group.** At a gap it sees the block
   built so far and exactly **one** locus — the first one after the gap. Whatever
   follows that locus has not been read yet.
2. **It is one-sided.** The `min_anchors` test is `len(set(_block_anchor_ranks(block))) >= 3`
   — the anchors on the **left** of the gap only. The right side contributes one gene.
3. **Direction is the block's own**, from `_longest_collinear_run` over its anchors
   (LIS in both directions, ties → `'+'`). `'+'` ⇒ frontier `max(anchors)`, candidate
   must satisfy `frontier < cand <= frontier + 5`; `'-'` ⇒ mirror image.

Also worth knowing, because the figure implies bridging is the only thing standing
between a gap and a dead zone: `cluster_distance = 150 kb` cuts, then
`merge_synteny_blocks` re-joins any two windows within `region_padding = 150 kb` of
each other. **So a cut only becomes a real dead zone for gaps between ~300 kb and
6 Mb.** Below 300 kb the merge fixes it; above 6 Mb bridging was never going to fire.

---

## 2. Do figures 1–4 match the code?

**Yes — every claim I checked is accurate.** Spot-checked against source:

| Figure | Claim | Verdict |
|---|---|---|
| 3 | the four bridge/refuse cases | **executed all four — all reproduce exactly** |
| 3 | cut at 150 kb; bridge < 6 Mb; ≥3 anchors; rank jump ≤5 | correct (`nextflow.config:77,130-132`) |
| 3 | direction is the block's own LIS, inversions bridge | correct — `F E D …gap… C B A` bridges |
| 3 | HIGH needs collinear run ≥3 once flank ≥3 | correct (`CLASSIFY_THRESHOLDS`, `collinear_ok`) |
| 4 | full classifier tree + every threshold | correct, value for value |
| 4 | `id ≥ 30` floor on the strong-flanking clause | correct (`fallback_strong_min_identity_floor = 30.0`) |
| 2 | weights 0.4/0.3/0.3, `+0.15` bonus, `0.45` floor | correct |
| 2 | `p = (#{null≥obs}+1)/(n+1)`, n=200, seed=42 | correct (`cluster_grs.py:894,1181`) |
| 2 | sort by `(−B, −final_score, +p)` | correct (`cluster_grs.py:1207`) |
| 1 | 10 flanking **per side**, 20 total | correct as a *default* — see gap C.3 |

I verified the rank-gap boundary by brute force: with block `A B C` (frontier rank 2),
candidates D…H (jumps 1–5) bridge and I, J (jumps 6–7) do not. Exactly as documented.

### Three gaps between the figures and reality

**C.1 — fig3 and fig1 overstate the dead zone.** Neither mentions that
`region_padding` re-merges sub-300 kb gaps. A reader concludes every refused bridge
is a silent false negative; in fact only the 300 kb – 6 Mb band is at risk.

**C.2 — fig4 predates the AMBIGUOUS tier** (uncommitted work from this week). Its
classifier tree now has a fourth leaf.

**C.3 — fig1's "20 anchors" is the default, not the flagship run.** The melittin
benchmark's `initial_db_locus_1.faa` carries **10 distinct flanking genes**, not 20
(the GOI-similarity filter removes the melittin-family neighbours). That number is
`N` in fig2's `C = max(U,B)/N`, so it is load-bearing for every coverage score in the
paper's headline run.

### One code defect found while checking (not a figure error)

`_longest_collinear_run` uses a **non-decreasing** LIS (`seq[j] <= seq[i]`), so
repeated identical ranks inflate the run: `[5,5,5] → 3`. A block that is one gene hit
four times therefore reports `collinear_chain_len = 4`. Confirmed on real data —
oskar `NC_045757.1`, `genes=['gene-Dmel_CG31100']`, `loci=4`, `chain=4`.

**Scope: block ranking and cap survival only.** I checked whether it reaches the HIGH
gate and it does not — `process_region_block` (line ~3095) builds `ordered_ranks` from
`flank_first_pos`, **one representative position per gene**, so the classifier's
`collinear_support` is correctly deduplicated. The fig4 "biglycan guard" is sound. But
`identify_synteny_blocks` sorts blocks by the inflated number, so a single repeated
gene can outrank three genuinely ordered ones.

---

## 3. Is the objection right?

**Yes, and the real problem is sharper than "min_anchors is too strict".**

Executed comparisons, same 5 genes (A B C D E), perfectly collinear, gap fixed at
~2 Mb, only the **position of the gap** moving:

```
A B   …gap… C D E     ->  CUT        <-- 5 collinear anchors, refused
A B C …gap… D E       ->  BRIDGED    <-- same 5 genes, admitted
A B C D …gap… E       ->  BRIDGED
A     …gap… B C D E   ->  CUT
```

The evidence is identical in all four; only the accident of where the rearrangement
fell decides. That is not a conservatism trade-off, it is an artefact of testing the
left side alone.

The same asymmetry breaks the case you *want* refused:

```
A B C …gap… A A D     ->  CUT       (duplicates first — refused, as intended)
A B C …gap… D A A     ->  BRIDGED   (same genes, duplicates second — admitted)
A B C …gap… D J E     ->  BRIDGED   (rank 9 spliced mid-run — never inspected)
```

Only the **first** locus after the gap is ever tested. `D` passes, then `A A` and `J`
join by ordinary 150 kb proximity. So the current rule refuses `A B …gap… D E F` on
thin evidence and admits `A B C …gap… D A A` on none.

Turning the knob does not fix this. Measured on the melittin benchmark (19 genomes):
`min_anchors 3→2` changes exactly one genome; `max_rank_gap 5→10` and `5→20` change
nothing at all. The shape of the rule is wrong, not its constants.

---

## 4. Would a two-sided rule be better?

Prototyped: at a gap, assemble the **right cluster** (loci after the gap joined by
normal proximity), and bridge iff the two sides *together* carry ≥ `min_anchors`
distinct anchors, their combined target-order rank sequence has a collinear run ≥
`min_anchors`, and the rank step across the gap continues that direction by ≤
`max_rank_gap`.

**Synthetic — it reproduces the stated intuition on every case:**

| case | current | two-sided |
|---|---|---|
| `A B C …gap… D E F` | BRIDGE | BRIDGE |
| `F E D …gap… C B A` | BRIDGE | BRIDGE |
| `D E F …gap… C B A` | cut | **BRIDGE** |
| `A B   …gap… D E F` | cut | **BRIDGE** |
| `A B   …gap… D`     | cut | **BRIDGE** |
| `A B C …gap… A A D` | cut | cut |
| `A B C …gap… D A A` | **BRIDGE** | **cut** ← fixes a current false admit |
| `A B C …gap… J`     | cut | cut |

**Real data — the change is small and auditable:**

- melittin, 19 genomes: **59 → 57 blocks**; bridges 1 → 3 (adds *Nomia melanderi*,
  *Osmia bicornis*; *Colletes* already bridged).
- oskar ladder, 2 divergent genomes: **no change at all** (22/21 blocks either way).
  Their flanking land as isolated singletons — most blocks hold one distinct anchor,
  so no rule keyed on collinear continuation can fire. Worth stating plainly: at
  oskar-level divergence, bridging is not the mechanism doing the work.

*Caveat:* the melittin reconstruction feeds `identify_synteny_blocks` the flanking
models that **survived** filtering and were emitted to the region GFFs, not the raw
MMseqs loci. It is a lower bound on how often bridging could apply.

### The cost, which is the part that decides this

- searched span across the 19 genomes: **2.47 Mb → 4.77 Mb (1.93×)**, from two
  genomes alone (Nomia 0.08→0.84 Mb, Osmia 0.09→1.64 Mb).
- `BlockFlankingSupport` on those blocks: **6 → 8**.

That second number is the problem. `flanking_support` is a property of the **block**,
and `fallback_strong_min_flanking = 5` gates a MEDIUM call on it — the exact mechanism
fig4 names for the yeast STE2 blow-up (32 MEDIUM at 21–27 % identity). **Bridging more
aggressively directly feeds the over-calling problem**: every junk hit anywhere in the
enlarged block inherits the stronger flanking context.

So this is not a free recall win. It trades a silent false negative for louder false
positives, in a pipeline whose measured weakness is already precision
(`locus_precision = 0.0755` on the current benchmark).

---

## 5. Recommendation

**Adapt it — the two-sided rule, not a constant change — but sequence it after the
over-calling work.**

1. **Now, independent of bridging:** fix the duplicate-rank inflation in
   `longest_collinear_run` for block *ranking* (count distinct ranks, or dedupe before
   the LIS). Low risk, no effect on the classifier, removes a case where one repeated
   gene outranks three ordered ones.
2. **Now:** correct fig3/fig1 for the `region_padding` re-merge (C.1), and fig1's
   anchor count for the melittin run (C.3).
3. **After AMBIGUOUS + phylo placement are landed and measured:** implement the
   two-sided rule behind `--synteny_bridge_two_sided`, **default off**, and measure
   HIGH/MEDIUM/AMBIGUOUS counts and §1y coordinate precision with it on and off. Flip
   the default only if precision holds. The reason to wait is that the AMBIGUOUS tier
   is the thing that will absorb the extra flanking-inherited calls this change
   produces; shipping the bridge first means shipping the false positives naked.
4. **One guard to add when implementing:** the prototype bridges `A …gap… D E F` (one
   anchor on the left). The union is collinear so it is defensible, but it lets a
   single hit extend a window by up to 6 Mb. Cap total bridged block span, or require
   ≥2 anchors on the side that grows.

If the rule changes, fig3's four-case panel must be redrawn — three of its four rows
would flip or change reasoning.

---

## 6. What was implemented (2026-09-09)

All of items 1, 3 and 4 above, plus the F9 blocker found while checking. **Uncommitted.**

| Change | File | Default |
|---|---|---|
| Two-sided bridge rule (`_right_cluster`, `_can_bridge_two_sided`) | `bin/iterative_search_runner.py` | **off** (`--synteny_bridge_two_sided`) |
| `bridge_max_per_block` runaway guard | same | **2** |
| Duplicate-rank inflation fix in block ranking | same | always on |
| Params + module wiring | `nextflow.config`, `modules/iterative_search.nf` | — |
| 13 tests | `tests/test_synteny_collinearity.py` | 17 → 30 pass |

The block walk became index-based so a bridge can consume the whole right cluster at
once. The legacy path is untouched and still the default; `test_two_sided_is_off_by_default`
pins that.

**Verified on the production code path** (not the prototype): melittin 59 → 57 blocks,
2 of 19 genomes changed, span 2.47 → 4.77 Mb, max flanking 6 → 8; oskar unchanged.
Identical to the prototype numbers in §4.

### The F9 blocker found on the way

The first F9 cluster run (job 5777904) reported
`evaluated 8, concordant 0, discordant 0, insufficient_data 8` — no verdicts at all.
Not a refutation: `_ID_SUFFIX` required a token *after* the locus id, so the tandem-copy
form `GOI_copy_1|Apis_cerana_fna_b0_l1` never matched and resolved to a species absent
from `sorted_genomes.txt`. **7 of 8 calls silently lost their species distance**, leaving
1 fitted against `min_calls=4`.

Three fixes in `bin/phylo_placement_check.py`:

- **resolve by longest known genome stem** instead of parsing arbitrary suffixes — the
  stems are the keys of `sorted_genomes.txt`. This also covers the hull-rescue
  (`…_<chrom>_hull_rescue`) and strong-synteny-rescue (`GOI_rescue_<genome>_locus_…`)
  forms, which the regex never handled at all — and fig4 attributes **6 of 7 HIGH calls**
  to the hull rescue.
- **fail loud** when a call cannot be matched, listing the unmatched ids and the known
  stems. The whole failure mode was that this degraded silently while reporting success.
- **`--min_target_len 20`**: the one call that *did* resolve was a **7 aa** protein at
  40 % identity. Over 7 residues that is noise, and one such point can set the slope.

Also added `--min_distance_levels 3`. `PHYLO_SORT` emits near-constant sentinel
distances when the taxonomy walk cannot separate targets — `local_runs/melittin_full_qw`
has 14 of 19 genomes at exactly 999, and `melittin_reg` runs 998–1000. The existing
guard only required ≥2 distinct values, which such a set passes while carrying no signal.
The 19-genome benchmark has a real spread (2, 6, 7, 9, 10, 12, 13), so F9 *is* measurable
there — it was only ever the parsing that failed.

## 7. The A/B run

LRZ job **5777925** (`bridge_ab`), submitted 2026-09-09, queued behind the Ant_Venoms
`newtaxa_synteny` run. Both arms of the 19-target melittin benchmark in one job — arm A
current rule, arm B two-sided, everything else identical — and both carry the F9 fix, so
the phylo measurement is re-taken at the same time.

What decides whether B ships:

- HIGH up **and** §1y coordinate recall up → B found something real.
- MEDIUM/AMBIGUOUS up with coordinate precision flat or down → B is manufacturing calls
  out of inherited block flanking support. Do not ship.
- everything flat → no-op on this clade; keep the default off and re-test on a rearranged
  set (decorin / vertebrates), which is where bridging is supposed to matter.

The third outcome is the likely one on melittin: only 2 of 19 genomes change offline, and
neither is a species where melittin is missing. Melittin is a weak test of a rearrangement
fix — worth stating before the numbers arrive rather than after.

### Result (job 5777925, 2026-09-09)

| metric | A (current) | B (two-sided) | delta |
|---|---:|---:|---:|
| high_confidence_goi | 2 | 2 | 0 |
| medium_confidence_goi | 5 | 4 | −1 |
| ambiguous_goi | 17 | 17 | 0 |
| low_confidence_goi | 71 | 71 | 0 |
| §1y locus_recall | 0.04 | 0.04 | 0 |
| §1y locus_precision | 0.0566 | 0.0566 | 0 |

Block construction behaved exactly as designed: arm A 415 blocks / 1 bridge event, arm B
414 blocks / 4 bridge events. The coordinate score is **identical to the digit** in both
arms, and the only headline change is one Nomia MEDIUM disappearing — a `tandem_goi_copy`,
not an ortholog. **On melittin the change is a no-op**, as predicted.

**But the result is not interpretable, and the reason matters more than the result.**

1. **Arm A did not reproduce the previous run of the same benchmark** (job 5777904,
   `bench_f9`). 17 of 20 per-genome GFFs differ, and the MEDIUM set is almost entirely
   different — Polistes, Bombus terrestris, Tetragonula, Formica out; Euglossa, Megachile
   in. §1y `locus_recall` fell 0.12 → 0.04 and `false_positive_calls_on_lost_loci` rose
   1 → 2. Yet **block construction was byte-identical between them** (415 blocks, 1 bridge
   event each), and the two search-core changes are provable no-ops there: the rank fix
   cannot bite because **no melittin block contains a duplicated flanking gene**
   (`loci_count > genes_count` in 0 of 59 blocks), and `bridge_max_per_block` never fired
   (1 bridge event, cap 2).
2. **Arm A vs arm B differ in 8 of 20 GFFs**, where the offline replay predicted 2 and the
   logs show only 3 extra bridge events. Nomia and Osmia are both in the differing set —
   along with six genomes the bridge should not have touched.

Either the pipeline is not reproducible run to run, or something in the sync changed
behaviour without changing block construction. Both fit the evidence; nothing here
separates them. Job **5778368** (`det_rep`) runs the same arm **twice with byte-identical
flags** to settle it.

If R1 ≠ R2, then no melittin benchmark delta of this magnitude is interpretable — which
would retroactively apply to the §1y recall/precision numbers already reported, and would
mean the AMBIGUOUS-tier and phylo measurements need replicates before they can be read as
effects. That is a bigger finding than the bridge rule.

## 8. The pipeline is not reproducible (job 5778368, 2026-09-09)

**R1 ≠ R2.** Two runs, byte-identical flags, same job, back to back:

| | R1 | R2 |
|---|---:|---:|
| per-genome GFFs identical | **3 of 20** | **3 of 20** |
| HIGH / MEDIUM / AMBIGUOUS / LOW | 2 / 6 / 15 / 71 | 2 / 5 / 17 / 71 |
| §1y locus_recall | 0.12 | 0.04 |
| §1y locus_precision | 0.08 | 0.0566 |
| curated_models_hit | 1 | 0 |
| median_overlap_frac | 0.1544 | 0.0 |
| false_positive_calls_on_lost_loci | 1 | 2 |

R1 reproduces the old `bench_f9` numbers exactly; R2 reproduces `bridge_A` exactly. **The
two "different" results were two draws from one distribution.** Nothing in §6 changed
behaviour — the code is exonerated, and so is the A/B in §7.

### Where it enters

Everything up to and including block construction is deterministic. On *Tetramorium*,
both replicates logged: `Parsed 200 hits`, `Found 35 discrete syntenic blocks`,
`Block filter … 2/35 retained`, `Smith-Waterman found 6 additional hits` (b0 **and** b1),
`Flanking gene summary: 9/10 parents found`. Run-wide: 3779 parsed hits and 415 blocks in
both. The **first** divergence is the line after:

```
R1: Expansion payload: 0 GOI-derived / 12 total annotations (3 GOI-like withheld)
R2: Expansion payload: 0 GOI-derived / 13 total annotations (4 GOI-like withheld)
```

Same hits in, different models out. ~~The consequence is not cosmetic: the true
*Tetramorium* melittin (`OV788322.1:15,634,953`) is modelled in R1 and **absent in R2**.~~
**Corrected 2026-09-10 (see "Resolution" below):** it is modelled in **both**: R1 as a fallback mRNA
`15,634,953-15,635,084` (AMBIGUOUS), R2 as a tandem copy `15,634,947-15,635,084` (LOW). The
tandem copy is a `gene` feature, and `score_coordinates.py` read only mRNA, so R2 scored as a
miss. The difference is real (a model-type and confidence flip), but it is not a gain or a
loss of the gene.

### Prime suspect

`bin/iterative_search_runner.py:3225` iterates `unique_queries` — a **set of strings** —
and appends to `found_queries` in that order; the downstream parent-collapse and spatial
filtering are order-dependent. CPython randomizes the string hash seed **per process**, and
`PYTHONHASHSEED` is set nowhere in the repo, `nextflow.config`, or the `Dockerfile`. One
`python3` invocation per ITERATIVE_SEARCH task means one seed per run, inherited by the
forked `ProcessPoolExecutor` workers — which predicts exactly "stable within a run, varies
between runs".

This is a mechanism consistent with every observation, **not yet a proven cause**. Job
`hashseed_probe` pins `PYTHONHASHSEED=0` and re-runs the same two replicates with **no code
change**, so the answer is attributable: R1 == R2 confirms it, R1 != R2 sends the hunt to
tool threading instead. (Submission blocked 2026-09-09 — Slurm controller unreachable;
scripts staged on LRZ, resubmit with `sbatch hashseed_probe.sbatch`.)
**Superseded 2026-09-10:** the cause was reproduced and fixed locally, so do not submit the
probe. See "Resolution" below.

### What this invalidates

Any melittin-benchmark delta of the size reported in this session is noise-dominated:
`locus_recall` alone swings 0.04–0.12 (3×) between identical runs. That covers the §1y
numbers, the AMBIGUOUS-tier counts, and the bridge A/B's −1 MEDIUM. The A/B's *structural*
result (415→414 blocks, 1→4 bridge events, coordinate score identical to the digit) is
still sound, because block construction is the deterministic part.

The F4 wavefront conclusion was drawn from a "byte-identical output" comparison and should
be re-checked against this: either those runs were genuinely identical (and the variance is
config-dependent) or that comparison needs a replicate too.

### Resolution (2026-09-10): the GOI representative collapse

**Cause.** `process_region_block` keeps one GOI query per parent ID (`extract_base_gene_id`
splits on `|`). After wave 2, the expanding DB holds the *A. florea* model
`GOI_Melt|Apis_florea_fna_b0_l1_exon_ann`, which shares parent `GOI_Melt` with the home
query. **Both are exactly 70 aa.** The collapse kept the longer sequence, first-seen on a tie, and
"first" came from iterating `unique_queries`, a set of strings, i.e. from the hash seed.
`unique_queries` is rebuilt for **every search block**, and its contents, and so its
iteration order, depend on that block's hits. So the tie is resolved **per block**: each block of
each post-wave-2 genome searched with either the *A. mellifera* or the *A. florea* melittin,
and a genome with *k* affected blocks has up to 2^*k* outcomes. The prime suspect named above
was the right line; the precise mechanism is the length tie in the collapse that follows it.

**Evidence, strongest first.**

1. *A prediction the mechanism makes, checked against LRZ job 5778368.* Only genomes searched
   before the expansion can be immune. The 3 GFFs identical across R1/R2 are exactly
   *A. cerana* (wave 1), *A. florea* (wave 2) and the home. All 17 post-expansion genomes differ.
2. *The tie resolves by hash seed.* The old rule, run on these exact IDs under
   `PYTHONHASHSEED=0..11`, keeps the home melittin in 9 seeds and the *A. florea* model in 3.
3. *No expansion, no flip.* Local 2-target replays (*Bombus terrestris*, *Tetramorium*), where
   nothing is added to the DB, gave functionally identical output under seeds 0, 0, 1, 2, 3 and
   the original random seed. The only byte difference was a cosmetic `Rearranged_from=` join
   order, also a set.
4. *Local 19-target run (random seed, pre-fix):* 2 HIGH / 5 MEDIUM / 17 AMBIGUOUS / 71 LOW,
   identical to LRZ R2, so this machine drew the same coin.
5. *Causation, by forcing the choice.* Reduced-wavefront replays of ITERATIVE_SEARCH used the
   19-target run's exact inputs, with *A. cerana*, then *A. florea* + *Bombus*, then
   *Colletes* + *Tetramorium*. Two otherwise identical copies of the fixed code differed only in
   the tie-break. **Home melittin wins:** *Tetramorium* gets `tandem_copy gene
   15,634,947-15,635,084` LOW, 13 annotations (4 withheld). That is LRZ R2. **A. florea
   forced:** `fallback_hits mRNA 15,634,953-15,635,084` AMBIGUOUS. That is LRZ R1's melittin
   model exactly. *Colletes* and *Tetramorium* (post-expansion) differ functionally in GFF,
   proteins and homology; *A. cerana*, *A. florea* and *Bombus* (pre-expansion) are
   byte-identical. The forced arm gives 15 annotations (6 withheld), not R1's 12 (3), because
   R1 was a *mixed* draw (item 6).
6. *The old code under seeds 0–5 reaches all three outcomes, and they decompose by block.*
   Seeds 0 and 4 give home everywhere: 13/4, byte-identical to the home arm, which is R2.
   Seed 3 gives *A. florea* everywhere: byte-identical to the forced arm (so no other
   hash-ordered path mattered). Seeds 2 and 5 give **12/3 with the fallback mRNA, exactly R1**.
   There, *Tetramorium*'s block-0 GOI models are identical to the forced arm's block 0 and its
   block-1 models to the home arm's block 1. That is the per-block resolution, observed directly.
7. *The fix holds, including under the seeds that flipped the old code.* Post-fix seeds 0, 1,
   2 and 3 give all 19 output files (GFFs, proteins, homology, expanded DB) byte-identical, and
   identical to pre-fix seeds 0 and 4. Seeds 2 and 3 are the ones that drove the old code to the
   R1-mixed and all-*florea* outcomes. So the fix locks in LRZ R2's outcome, is invariant to the
   hash seed on this path, and changes nothing else.

**Summary of the replays** (*Tetramorium*, reduced wavefront; H = home melittin, F = *A. florea*):

| code | seed | block 0 | block 1 | annotations (withheld) | melittin locus model |
|---|---|---|---|---|---|
| old | 0, 4 | H | H | 13 (4) = LRZ R2 | `tandem_copy` gene, LOW |
| old | 2, 5 | F | H | 12 (3) = LRZ R1 | `fallback` mRNA, AMBIGUOUS |
| old | 3 | F | F | 15 (6) | `fallback` mRNA, AMBIGUOUS |
| fixed | 0, 1, 2, 3 | H | H | 13 (4) | `tandem_copy` gene, LOW |

Reproduce: copy the ITERATIVE_SEARCH task's `.command.sh` from a run's work dir, relink its
five inputs into a fresh directory, replace `filtered_sorted_genomes.txt` with the reduced
list (`Apis_cerana.fna 2`, `Apis_florea.fna 2`, `Bombus_terrestris.fna 7`,
`Colletes_gigas.fa 9`, `Tetramorium_bicarinatum.fna 13`), and run it with
`PYTHONHASHSEED=<n>` and `synvoy_env` first on PATH (~6 min each).

**Fix.** `collapse_goi_queries_by_parent()` breaks ties explicitly: the longest sequence wins;
on equal length the user's own query (id == parent) beats an expansion model; then the smallest
id. `unique_queries` and the GOI-proxy sweep (a real second order-dependent site: each
sweep marks its locus covered) now iterate `sorted()`. The cosmetic sites (`Rearranged_from`,
block gene lists, family-token match reason, a `max()` tie in `plot_synteny.py`) are sorted too.
An AST audit of every `bin/*.py` for iteration, `list()`, `join`, `max`/`min`/`sorted` or
`pop` over sets found no other functional site. `PYTHONHASHSEED=0` is set in `nextflow.config`
`env {}` as defence in depth. `tests/test_a0_determinism.py` pins the tie-break and runs the
collapse under 12 seeds, with a control that asserts the raw set order *does* vary across them.

**Side effect for F4.** The same collapse means an expansion model reaches the *region-level*
search only when it is strictly longer than the home query (before this fix: or when it won
the hash tie). Its genome-wide hits still enter through `relevant_hits`. That is a plausible
structural reason why "the wavefront changes nothing" (STATE_OF_THE_PROJECT F4, job 5768082).
It is a code-reading hypothesis, not a measurement.

**Scorer bug found on the way.** `score_coordinates.py` read only `mRNA` GOI features, and
tandem copies are written as `gene`. On the local R2-equivalent run that hid 42 of 95 GOI
calls (6 MEDIUM) and moved locus_recall 0.04 → 0.28. Fixed, with a regression test. Every §1y
number computed before 2026-09-10 undercounts and should be rescored before it is cited.

## Reproducing

The probes live in the session scratchpad; the substantive ones are:
`identify_synteny_blocks` replayed on
`local_runs/oskar_ladder_5702852/hits/*.m8` + `initial_db_locus_1.faa` (exact inputs,
`parse_hits(min_identity=40, min_length=50, evalue=1e-5)`), and on flanking models
reconstructed from
`local_runs/melittin_val_sw_20260630/plot_inputs_synteny_block_locus_1/*.gff` with
`home_rank` from that run's `initial_db_locus_1.faa`.
