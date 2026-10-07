# SynVoy — State of the Project

**Compiled 2026-07-26 by direct code and output inspection; fix logs added through
2026-10-05.** Every number below was re-derived from the repository or from run outputs
on disk; nothing is quoted from an earlier summary without re-checking. Where a claim
could not be verified it says so.

**Read the fix logs first.** Parts A–F are the 2026-07/08 audit and are kept as written;
where a later log entry supersedes a statement, the entry says so. The newest state is:
[`ALGORITHM_AUDIT_2026-09.md`](ALGORITHM_AUDIT_2026-09.md) (search core),
[`NEXT_SESSION_FAMILY_BENCHMARK.md`](NEXT_SESSION_FAMILY_BENCHMARK.md) (the 10-family
benchmark), [`VIZ_REVIEW_2026-09.md`](VIZ_REVIEW_2026-09.md) (figures), and **Part G**
below (the 2026-10-05 robustness pass).

Companion to `docs/TODO.md` (the live task list). This file answers three questions:
**what is the contribution**, **what is actually broken**, and **what has to happen to
wrap up**.

> **Fix log — 2026-07-26.** Part F items 2, 4, 5 (half), 6, 7 and 8 are done: the
> anchor-grid identity labels, the invalid p-value (**F1**), the un-floored fallback
> clause (**F6**, first half), the unauditable ranking key (**§1u**), the silent discard
> (**§1x**), strand conservation (**F2** — which turned out to be scoring the *opposite*
> of what it claimed), the two rival collinearity measures (**F5**), and the wavefront
> binning (**F4**). **F6**'s remaining half was measured and closed won't-do. Suite:
> **550 passed, 1 skipped, 0 failed** (was 515 with 2 failed). Remaining open items are
> marked ⬜ below. **F4**'s code half is fixed too (rank wave binning is now the default);
> what remains there and in **§1y** are benchmark *runs*, not code.
>
> ⚠️ **F2 and F5 both change region scores.** Old `regions/*.scores.tsv` and region names
> (`Reg1_G7_CMEDIUM_S0.47`) are not comparable to new ones, and the melittin GT fixture
> needs regenerating — do it once, now that both landed in the same pass.

> **Fix log — 2026-08-21/22 (first cluster validation of the above).** The 2026-07-26 pass
> was committed (`49df417`, `8e44d55`) and then run on LRZ for the first time. Three runs,
> melittin `P01501`, easy mode, 4 auto-picked targets:
>
> - **F1 ✅, F6 ✅, §1u ✅ confirmed on real output.** 2 of 8 GOI-overlapping region rows sit
>   at the permutation floor (pre-fix: all of them); MEDIUM identities bottom out at 30.4 %
>   against the 30.0 floor; `synteny_score` is emitted.
> - **F4 ⚠️ measured, inconclusive** — both binnings gave identical waves on this target set.
>   The experiment that answers the question is *waves vs no waves* on a uniform-depth set.
>   See **F4**.
> - **§1x ⚪ null result** — 0 rejected candidates; every GOI hit was inside the flanking
>   envelope.
> - **F8 🔴 NEW, found and fixed.** Easy mode was searching for the **flanking proteins**
>   instead of the GOI in all three rescue/reciprocal passes. Every "high-confidence
>   ortholog" in the first run was a 10 aa flanking fragment at "100 %". Fixed
>   (`main.nf:1431`); the re-run recovers the true *Apis cerana* melittin at 97.1 % as the
>   single HIGH call. **This invalidated the headline of every easy-mode run ever produced.**
>   See **F8**.
>
> Part A gained **A.8** (the rescue layer), which had never been documented algorithmically
> despite being the source of the recovered ortholog in both validated cases.

> **Fix log — 2026-09-09 → 2026-09-17 (committed 2026-10-05 as `00e4853`…`1e70fa6`).**
> Three weeks of work that Parts A–F do not describe. In brief; details in the three
> documents named above.
>
> - **Runs are reproducible.** Two identical cluster runs had differed in 17 of 20 genome
>   GFFs. Cause: two equal-length GOI queries tied when one representative per parent was
>   chosen, and set order broke the tie. Fixed, `PYTHONHASHSEED` pinned, verified on the
>   cluster under different hash seeds (20/20 identical).
> - **A 10-family benchmark** on curated ant venom families (honeybee seed, 33 targets,
>   reachable set fixed and hashed in advance). Near-single-copy families: 30/31, 28/28,
>   32/32, 31/32 of reachable genes. Held-out families: no recall lost with the fixed
>   code. 4 of 10 curated seeds sit at the wrong genomic address (a curation finding).
> - **Four precision defects** (§19): tandem copies promoted without synteny context,
>   overlapping "exons" chained, coverage measured as span, fragments lifted by one
>   flanking gene. Every final call with query coverage < 0.5 had been off-truth.
> - **Two defects that reached every output** (§20): multi-exon model proteins were
>   translated out of frame behind phase-1/2 introns (coordinates right, proteins wrong —
>   so seeds, trees and the paralog check ran on corrupted sequences), and Smith-Waterman
>   hits carried a fake E-value (six chance hits per window). Both fixed.
> - **Model ends** (§21): gene ends exactly on the curated coordinate 28/182 → 97/182.
> - **Figures** (§22): flanking genes had been drawn with invented exon positions; fixed
>   at the source. ≤ 3 lanes, fragments hidden by default, GOI column from real CDS.
> - **AMBIGUOUS tier and F9** (2026-08-31/09-09): see F9 below. Measured on the benchmark:
>   the tier holds 0 of 137 records on a curated gene; the phylo check carries no signal.
>
> **What this changes in Parts A–F:** A.4's `tandem_copy` row now also needs context
> (block flanking ≥ 5 or query coverage ≥ 0.65); A.3's Smith-Waterman input now has real
> statistics; §1p.1 (invented flanking exons) is fixed; §1y's scorer and tests are in the
> repository and CI runs on `dev` (Part C said "frozen in CI" on 2026-08-30, but those
> files were not committed until 2026-10-05).

> **Fix log — 2026-10-05 (robustness pass, Part G).** A systematic check of the wiring
> between `nextflow.config`, the `.nf` files and the scripts, plus end-to-end runs. Found
> and fixed: the plot step did not parse on Python 3.10/3.11; `--auto_params` silently
> applied nothing on current Nextflow; an auto-selected preset overwrote values the user
> set and carried TP53 family names into unrelated queries; `--x false` was read as true
> in nine places; region scores depended on the line order of the hit table; rescue-pass
> gene models were never written to the output folder; 78 % of a cluster console log was
> one repeated warning. Suite: **867 passed, 0 failed**. See Part G.

---

## 0. One-paragraph status

SynVoy is a working synteny-guided search engine wrapped in a scoring, labelling and
reporting layer that is not trustworthy. The discovery layer is validated across four
clades and is robust to large parameter changes. The adjudication layer — region ranking,
confidence labels, the p-value, the report headline — contained at least four defects
serious enough that its numbers should not be published as-is; as of 2026-07-26 those four
are fixed, and the 2026-08-21 cluster run confirmed three of them on real output.

*(Status as of 2026-10-05: the paragraphs below are the 2026-08 picture. Since then the
search core was audited and the two output-wide defects in it fixed, the result became
reproducible, and there is a ten-family benchmark behind the recall claim. What has not
changed: on melittin SynVoy still places a call in most species including those that lost
the gene — those calls are now labelled AMBIGUOUS and not counted, which is honest but is
not the same as solving it — and the wavefront still has no measured effect.)*

**The single most important thing in this document is now F8**, found on 2026-08-21: in
**easy mode** — the documented path, the one the BA students run — the rescue and
reciprocal passes were searching for the **flanking genes** rather than the gene of
interest, and reporting what they found as high-confidence orthologs. On the melittin run
all three "high-confidence orthologs" were 10 aa fragments of one flanking gene at
"100 % identity", two of them in genomes (*Bombus*) that do not have the gene. It is fixed,
and the corrected run returns the true *Apis cerana* ortholog at 97.1 % as the single HIGH
call — but **every easy-mode headline produced before 2026-08-22 is invalid**, and any
figure or count derived from one must be regenerated.

Second in importance, and still open: the signature *mechanism* of the method (the
closest-first expanding wavefront) was **inoperative in 37 of 44 runs on disk, including
the flagship melittin benchmark** — so the published melittin result was produced *without*
the mechanism the paper would attribute it to. The binning is fixed; the measurement
proving the wavefront earns its place has still not been made (see **F4**).

---

# Part A — The contribution, stated as mathematics

This is the material for the paper and the talk. It is extracted from the code, not from
the docs, and each formula carries its source location.

## A.1 What the method actually is

> Locate the query's home locus → take its *n* flanking genes → find where those flanking
> genes land in each target genome → keep the collinear neighbourhoods → search for the
> query only inside them.

The claim is that **gene order is conserved longer than gene sequence**, so for a gene too
divergent for direct homology search, the neighbourhood is a stronger locator than the
sequence. Everything below is machinery in service of that one idea.

One caveat belongs in the summary rather than a later section: the line above says "keep the
collinear neighbourhoods", but the neighbourhood that finally yields the ortholog is often
**not** a collinear block. When the gene's immediate neighbours are rearranged it falls into
a *gap between* blocks, and it is the hull rescue (A.8) — the span of the whole neighbourhood,
gaps included — that recovers it. On both validated cases where the answer was known in
advance, cow decorin and *Apis cerana* melittin, the recovered ortholog came from the hull,
not from a seeded block. The honest one-line statement of the method is therefore
"**the flanking neighbourhood, not the collinear block, is the search unit**".

## A.2 Synteny block construction — `bin/iterative_search_runner.py`

Flanking-gene hits on one target chromosome are clustered by proximity: a new block starts
whenever the genomic gap to the next hit exceeds `cluster_distance` (150 kb).

**Collinearity measure (LIS).** For a block, let `r = (r₁ … r_k)` be the home-genome ranks
of its anchored flanking genes, taken in *target* order. `_longest_collinear_run`
(`:1749`) returns

```
L(r) = max( LIS≤(r), LIS≤(reverse r) )          direction '+' if the first wins, else '-'
```

i.e. the longest non-decreasing subsequence in either direction. Inversion-tolerant by
construction: a uniformly inverted neighbourhood scores as highly as a forward one.

**Gap bridging** (`_can_bridge`, `:1788`). Two clusters on one chromosome separated by a
gap up to `synteny_bridge_max_gap` (6 Mb) are kept as **one** block iff

```
|distinct anchors in block| ≥ synteny_bridge_min_anchors      (3)
direction d = direction component of L(anchors)
d = '+' :   frontier = max(anchors),  bridge iff  frontier < r_cand ≤ frontier + Δ
d = '-' :   frontier = min(anchors),  bridge iff  frontier − Δ ≤ r_cand < frontier
Δ = synteny_bridge_max_rank_gap                                (5)
```

What is refused is a **direction reversal across the gap**, not an inversion. This is the
fix that recovered cow decorin at chr5:21.0 Mb (90 %), which sat in a 19.7–22.3 Mb dead
zone between two seeded blocks.

**Block flanking support** (`:2988`):

```
S_block = |{ base_gene_id(h.query) : h ∈ block, h.query is not a GOI proxy }|
```

Note this is a property of the **block**, not of an individual call. See flaw **F6**.

## A.3 Region scoring — `bin/cluster_grs.py`

For a candidate cluster with hits mapped to home-gene ranks, `score_flexible_synteny`
(`:761`) returns three quantities:

```
U  = |{ distinct home ranks hit }|                            unique genes
C  = U / N                                                    coverage, N = n flanking expected
K  = L(r) / k                                                 order consistency
S  = max(#agree, #disagree) / #scored                         strand conservation
```

where `r` are the home ranks of the `k` ranked hits taken in **target** order, `L` is the
same longest-monotonic-run function the search core uses (A.2), and `agree` counts hits
whose target strand matches their home gene's. Both `K` and `S` take a max over both
directions, so a **uniformly inverted neighbourhood scores 1.0** — an inversion preserves
relative order and orientation and is conserved synteny; only scrambling is penalised.

*(Both definitions were corrected on 2026-07-26 — see flaws **F2** and **F5**. `K` was a
fraction of adjacent monotonic pairs, a different measure from the search core's; `S`
never consulted the home genome at all and scored scrambled neighbourhoods highest.)*

The composite (`:1055`):

```
Q             = w_b·C + w_k·K + w_s·S          w_b,w_k,w_s = 0.4, 0.3, 0.3
synteny_score = Q · C                          ← the measured statistic
score         = synteny_score
              + 0.15                           if the region overlaps a GOI call
score        := max(score, 0.45)               if the distant-macrosynteny rescue fires
```

Expanding: **synteny_score = 0.4·C² + 0.3·K·C + 0.3·S·C**. See flaw **F3** — coverage
dominates; this is not the balanced three-component score it presents itself as.

**Significance** (`estimate_pvalue`, `:838`) is a label-shuffling permutation test:
resample `k` query labels with replacement from all hits genome-wide, re-score with the
same `Q·C` formula, and report

```
p = (#{ synteny_null ≥ synteny_obs } + 1) / (n + 1),      n = 200, seed = 42
```

so `p ∈ [1/201, 1] = [0.00498, 1]`. The tested statistic is `synteny_score`, **not**
`score`: the GOI bonus and the rescue floor are prior knowledge we inject, and a null
that cannot generate them must not be asked to compete with them (flaw **F1**, fixed).
Both numbers are emitted in `regions/*.scores.tsv` so the test is auditable.

**Ranking** (`:1119`): `sort by (−S_block, −score, +p)`. Block synteny support outranks
score. Since §1u (2026-07-26) both `synteny_score` and `goi_block_flanking` (= `S_block`)
are emitted in `regions/*.scores.tsv`, so the ordering a user sees is now reconstructible
from the output files.

## A.4 Confidence classification — `_classify_goi_evidence` (`:935`)

A decision tree over `(evidence_type, identity, exon_count, query_cov, S_block, L_block)`.
Defaults from `CLASSIFY_THRESHOLDS` (`:199`):

| evidence_type | → HIGH | → MEDIUM |
|---|---|---|
| `exon_annotation` | `exons ≥ 2` ∧ `id ≥ 50` ∧ `S_block ≥ 2` ∧ collinear_ok | `id ≥ 35` ∧ (`S_block ≥ 1` ∨ `qcov ≥ 0.65`) |
| `fallback_hit_span` | — | (`S_block ≥ 2` ∧ `qcov ≥ 0.75` ∧ `id ≥ 60`) **∨** (`S_block ≥ 5` ∧ `id ≥ 30` ∧ (`qcov ≥ 0.25` ∨ `id ≥ 35`)) |
| `tandem_copy` | — | `id ≥ 40` ∧ `qcov ≥ 0.35` ∧ (`S_block ≥ 5` ∨ `qcov ≥ 0.65`) *(context term added 2026-09-11)* |
| `rescued_exon`, `raw_hit` | — | only via PLM/structural rescue |

where the **order gate** is

```
collinear_ok  ⇔  L_block is None  ∨  S_block < 3  ∨  L_block ≥ 3
```

A HIGH-quality model in a *scrambled* neighbourhood (enough flanking, wrong order) is
demoted to MEDIUM — this is the biglycan-as-decorin guard.

*(The second `fallback_hit_span` clause carried **no identity floor** until 2026-07-26:
with `S_block ≥ 5` and `qcov ≥ 0.25`, a hit at any identity reached MEDIUM — the mechanism
behind the yeast STE2 overcalling, 32 MEDIUM calls at 21–27 % identity. The `id ≥ 30` term
is the **F6** fix, calibrated so the 34.7 % ant-melittin rescue survives.)*

⚠️ **HIGH still has no query-coverage term in the classifier itself.** The
`exon_annotation` → HIGH branch tests exon count, identity and flanking support only, so a
short high-identity window can still be *labelled* HIGH where the `tandem_copy` branch would
reject it. Since 2026-08-22 this is caught one layer later instead: `apply_coverage_demotion`
(`generate_report.py`) relabels such a call `identity_coverage_decoupled` and drops it from
the headline counts, exactly as §1m does for `paralog_not_goi`.

The demotion fires **only when coverage is known and low** (`identity ≥ 50` ∧
`qcov < 0.35`, both tunable). A call whose coverage was never recorded stays put and remains
advisory — absence of evidence is not evidence of low coverage, and before
`rescue_goi_hull.py` learned to emit `QueryCoverage` that described *every* rescue model,
including the true *Apis cerana* melittin. Disable with `--disable_coverage_demotion`.

## A.5 Distance auto-tuning — `_apply_distance_adaptive_thresholds` (`:282`)

Let `m` = median best per-gene flanking identity for the target. With `c = 70`, `f = 40`,
`R = 10`:

```
        ⎧ 0                      m ≥ c        tier close  (exact no-op)
relax = ⎨ R·(c − m)/(c − f)      f < m < c    tier mid
        ⎩ R                      m ≤ f        tier far    (flag manual_review)

high_min_identity  := max(50 − relax, 1)
medium_min_identity := max(35 − relax, 1)
```

Mutates module-global state per genome; safe because genomes run in separate worker
processes and the baseline is re-derived first.

## A.6 The wavefront — `bin/iterative_search_runner.py:6549`

The design: process targets closest-first; each genome's recovered models are appended to
the MMseqs2 query DB, so a divergent species is searched with a *nearer relative's* actual
ortholog rather than only the original query. This is the "iterative" in the method.

Distances arrive from `PHYLO_SORT` as **integer taxonomic rank-steps**, then:

```
d_i := d_i / max_j(d_j)          applied whenever max > 1
```

**Binning — since 2026-07-26 by rank-quantile** (`assign_waves_by_rank`). With `n`
targets sorted closest-first and `q = i/(n−1)`:

```
q < 0.10   → wave of 1     (closest tier: strictly serial, maximum seed propagation)
q < 0.35   → wave of 2
q < 0.70   → wave of 3
else       → wave of 5     (farthest tier: parallelise)
```

This is a pure function of the sorted order, so it is deterministic and grades close→far
however saturated the distance metric is.

The legacy binning it replaced applied *absolute* cut-offs to the divide-by-max
normalised distance:

```
d < 0.05            → wave of 1
0.05 ≤ d < 0.15     → wave of ≤ 3, members within 0.01
d ≥ 0.15            → wave of ≤ 5, members within 0.02
```

which is incoherent by construction — see flaw **F4**. `--rank_wave_binning false`
restores it for reproducing older runs.

**Measured 2026-08-21 (LRZ A/B, jobs 5757307/5757308):** the two binnings are *not*
universally different. On a target set spanning three distinct taxonomic ranks — *Apis
cerana* (genus), two *Bombus* (family), *Agapostemon* (order) → normalised distances
`0.200 / 0.700 / 1.000` — both produce identical waves `[1, 2, 1]` and identical results.
Legacy binning fails specifically on a **uniform-depth** target set, where divide-by-max
collapses everything to 1.0. Rank binning is the safer default because it cannot collapse;
it is not a change of behaviour on well-spread sets.

## A.7 Post-hoc paralog adjudication

`bin/build_home_paralog_panel.py` aligns the query against the whole home proteome and
keeps homologs above `paralog_panel_min_identity` (28 %) with a coverage guard.
`bin/reciprocal_best_paralog_check.py` then Smith-Waterman-scores each recovered target
against every panel member; a call whose best home match is a family paralog rather than
the GOI is relabelled `paralog_not_goi` and excluded from the headline counts. This is
what demotes the cow/mouse "decorin" that is really biglycan (BGN 950 vs DCN 507).

*(Until 2026-08-22 there was a hole here: `ASSIGN_LOCUS_OWNERSHIP` consumes the
ITERATIVE_SEARCH region files, while rescue models (A.8) are mixed in later at
`STAGE_REGION_GFF` — so the one guard designed to catch a mislabelled GOI never saw the
models most likely to carry one. `rescue_goi_hull.py` now also emits the rescued model's
protein (`--output_faa`, translated from the CDS with strand/phase handling), and the
sub-workflow feeds it in as its own `(locus, genome.hull_rescue)` ownership task.
`canonical_genome_id()` already strips that suffix, so the rows fold back onto the real
genome.)*

## A.8 The rescue layer — the part that is not block-based

A.2–A.4 describe the main path: seed blocks from flanking hits, search inside them, classify.
Two passes run *outside* that path, and on the melittin run they were the only source of a
HIGH call — so they are part of the method, not a footnote.

**Strong-synteny rescue** (`bin/rescue_strong_synteny.py`, §1e Phase B). Fires on a block
`cluster_grs` classed `goi_missing_but_strong_synteny`: ≥ `strong_synteny_min_flanking` (5)
HIGH flanking genes and no GOI model at all. Runs a relaxed `miniprot` (`--outc 0.05`) over
the block window and emits a **LOW**-confidence model, `EvidenceType=relaxed_miniprot_rescue`.
Deliberately conservative: it asserts "something GOI-like is here", not orthology.

**GOI synteny-hull rescue** (`bin/rescue_goi_hull.py`, §1m/§17). Addresses the failure the
block model *cannot* see: a GOI sitting in a **gap between** flanking-seeded blocks, because
its immediate neighbours were rearranged. Cow decorin at chr5:21.0 Mb fell in the dead zone
between blocks at 17.8–19.5 and 22.5–24 Mb and was simply never searched.

Rather than a block, it takes the **hull** of the neighbourhood (`compute_hull`): find the
dominant flanking cluster (single-linkage at `goi_hull_cluster_max_gap`), take its centre,
and keep every flanking gene whose midpoint lies within `goi_hull_max_window/2` of it —

```
centre = midpoint(dominant cluster)
hull   = span{ f ∈ flanking : |midpoint(f) − centre| ≤ goi_hull_max_window / 2 }
window = hull ± goi_hull_window_pad          (skipped if wider than max_window)
```

so the window spans the gap while still excluding a far rearranged outlier (cow chr5's
105 Mb singleton, ~80 Mb away). It fires only when the hull holds ≥ `goi_hull_min_flanking`
(4) HIGH flanking genes **and no HIGH GOI model already overlaps it** — so it is a no-op
wherever the main path already succeeded. Relaxed miniprot again, then the best model by
identity that clears both gates:

```
accept  ⇔  identity ≥ goi_hull_min_identity (40)  ∧  cov ≥ goi_hull_min_coverage (0.5)
cov     =  Σ CDS length / 3 / |query|
conf    =  HIGH   if identity ≥ classify_high_min_identity   (50)
           MEDIUM if identity ≥ classify_medium_min_identity (35)
           LOW    otherwise
```

Note what confidence is **not** a function of here: flanking support, collinearity, or exon
structure. Inside the hull, identity alone decides the label — the coverage term is a gate,
not a grade. That is defensible (the hull already established the neighbourhood) but it means
the rescue path skips every order-based guard A.4 applies, including the collinearity gate
that exists to stop a paralog being called the GOI.

**Empirical weight.** On the melittin LRZ run (job 5757853, post-F8) the single HIGH call —
*Apis cerana* melittin, `NC_083855.1:12,201,521`, 135 aa, 97.1 % — came from the hull rescue,
not the main block path. The mechanism that recovers the flagship result is therefore the
hull, and this is worth stating plainly in any write-up.

---

# Part B — Flaws found in this audit (new, verified, not in TODO.md)

These are additional to the §1x–§1z items already tracked. Each was reproduced.

## F1 ⭐⭐ The p-value is invalid for every GOI-overlapping region — ✅ FIXED 2026-07-26

**`bin/cluster_grs.py:1061-1090`.** The observed score has `+0.15` (GOI overlap) and/or a
`max(·, 0.45)` floor (distant rescue) applied **before** it is handed to `estimate_pvalue`.
The null model recomputes `Q·C` and never applies either term. An inflated observed value
is therefore tested against an un-inflated null.

Controlled reproduction — same cluster, only the bonus toggled:

| | observed score | p |
|---|---|---|
| no GOI overlap | 0.2160 | **0.194** |
| GOI overlap (+0.15) | 0.3660 | **0.005** ← the floor |
| distant-rescue floor | 0.4500 | **0.005** ← the floor |

In the real outputs (130 region rows across the `results/demo_*` runs):

- **55 rows sit exactly at the p-value floor** `1/201 = 0.004975`.
- **All 55 have `goi_overlap=True`. Zero of the 21 `goi_overlap=False` rows reach the
  floor.** Perfect separation — the signature of the bug, not of biology.

**Fix.** `cluster_grs.py` now separates the two quantities. `synteny_score = Q·C` is the
pure neighbourhood statistic and is the only thing `estimate_pvalue` sees; `final_score`
keeps the bonus and the rescue floor and is used for ranking and display. The bonus is
prior knowledge we inject, not evidence the neighbourhood provides, so a null that cannot
generate it must not be asked to compete with it.

Verified after the fix — same cluster, all three cases now agree, as they must:

| case | ranked score | p |
|---|---|---|
| no GOI overlap | 0.2160 | 0.194 |
| GOI overlap (+0.15) | 0.3660 | **0.194** |
| distant-rescue floor | 0.4500 | **0.194** |

`synteny_score` is now emitted in `regions/*.scores.tsv` next to `p_value`, so the test
statistic is auditable from the output. Regressions: `tests/test_adjudication_fixes.py`
(including a test that pins the *old* behaviour, so a reintroduction is unmistakable).

## F2 ⭐⭐ `strand_consistency` measured the *opposite* of strand conservation — ✅ FIXED 2026-07-26

**`bin/cluster_grs.py:831-834`.** `S = max(#'+', #'−') / k` over the *target* hits. The
home strand never enters the function. Verified: a cluster with every hit on `+` and one
with every hit on `−` both score `S = 1.000`. An entirely inverted neighbourhood is scored
as perfectly strand-conserved.

Second consequence: `S ∈ [0.5, 1.0]` by construction (confirmed empirically: minimum over
20 000 random clusters = 0.500). It contributes a near-constant `0.15–0.30` to `Q`. It is
documented as "Weight for strand conservation" and it is neither.

**It was worse than useless — it was inverted.** Measured on a 6-gene home neighbourhood
with alternating strands:

| cluster | legacy `S` | fixed `S` |
|---|---|---|
| conserved (matches home exactly) | **0.500** | **1.000** |
| uniformly inverted (every strand flipped) | 0.500 | 1.000 |
| scrambled (all forced to `+`) | **1.000** | **0.500** |
| half scrambled | 0.833 | 0.667 |

The legacy measure gave the *scrambled* neighbourhood the maximum and the *conserved*
one the minimum, because a scrambled-but-strand-uniform cluster is maximally uniform. A
term carrying weight 0.3 of the quality score was actively rewarding the thing it was
meant to penalise.

**Fix.** `load_home_strands` reads column 6 of the synteny BED (already present — it was
simply never loaded), and `score_flexible_synteny` now scores
`max(#agree, #disagree) / #scored` against the home strand. Taking the max of both keeps
a **uniform inversion at 1.0** — an inversion preserves relative orientation and is
conserved synteny, exactly as `_longest_collinear_run` treats gene order in both
directions — while a scrambled neighbourhood drops. The permutation null uses the same
measure, so observed and null remain commensurable. When the BED carries no strand
column the scorer falls back to the legacy value and says so loudly rather than
fabricating a conserved score. `--legacy_strand_score` reproduces old numbers.
Regressions: `tests/test_adjudication_fixes.py` (including one pinning the inversion).

⚠️ **This changes every region score**, so `regions/*.scores.tsv` and region names
(`Reg1_G7_CMEDIUM_S0.47`) from earlier runs are not comparable to new ones. The melittin
GT fixture will need regenerating — that is the intended "deliberate behaviour change"
trigger its README documents.

## F3 ⭐ The region score is effectively coverage-squared — ⬜ description fixed, formula unchanged

Because `Q` contains `w_b·C` and is then multiplied by `C`, and because `S` is pinned to
`[0.5, 1]`:

| coverage C | achievable score range | spread from K and S |
|---|---|---|
| 0.2 | 0.076 – 0.136 | 0.060 |
| 0.6 | 0.324 – 0.504 | 0.180 |
| 1.0 | 0.700 – 1.000 | 0.300 |

Coverage determines the score; order and strand modulate it by at most `0.3·C`. Presenting
this as a weighted three-component synteny score overstates it. If the paper describes the
scoring function, describe it as **coverage-dominated with order/strand as a modifier** —
that is defensible; "0.4/0.3/0.3 weights" is not.

**Deliberately not "fixed".** Rescaling `K` and `S` from `[0.5, 1]` onto `[0, 1]` would
make the weights mean what they say, but it changes the magnitude of every score ever
produced for a presentational gain, on top of the F2 change that already shifts them. F3
is a *description* problem and is now handled by describing it correctly. Revisit only if
you decide to renumber the scale anyway — in which case do it in the same pass as
regenerating the benchmark fixtures.

## F4 ⭐⭐⭐ The wavefront is inoperative in 37 of 44 runs — ✅ FIXED, and ✅ MEASURED 2026-08-31: it changes nothing

Recomputed the binning from every `sorted_genomes.txt` in `results/` and `local_runs/`
(44 runs):

- **37/44** never place their closest genome in the serial (`d < 0.05`) tier.
- The 7 that do all contain a literal `dist = 0` entry.
- The **melittin ground-truth runs** (`melettin_gt_v2 … v7`) have raw distances
  `[998, 999, 999, 999, 1000]` → normalised to `[1.0, 1.0, 1.0, 1.0, 1.0]` → **one single
  wave of five genomes, fully parallel**.

Confirmed in the run log itself:

```
results/melettin_gt_v7_smoke/logs/...ITERATIVE...  →  "Defined 1 waves of execution."
```

Two separate causes compound: distances of 998–1000 mean the taxonomy lookup fell back to
its unknown-sentinel, and divide-by-max then collapses any uniform target set to 1.0.

**Why this matters more than any other item here:** the expanding database is the
mechanism the method is *named* for. On the flagship benchmark it never iterated — every
target was searched with the initial DB only. Melittin still scored 11/12 vs 2/12 for the
baselines, which is good news for the result and bad news for the explanation: **the win
came from the flanking-neighbourhood restriction alone, not from the wavefront.**

**Root cause.** The binning tests *absolute* cut-offs (0.05 / 0.15) against a distance
that has just been divided by its own maximum. That is incoherent by construction: the
farthest target is always exactly 1.0, so the closest reaches the serial tier only if it
is ≥ 20× nearer. On a uniform target set — the common case, since you pick relatives at
similar depth — every distance collapses to 1.0 and the graded wavefront degenerates to
one parallel wave.

**✅ Code fixed 2026-07-26.** `rank_wave_binning` (§A4, already implemented and tested)
bins by phylo-distance *rank* instead, which grades regardless of how saturated the metric
is. It is now the **default**; `--rank_wave_binning false` restores the legacy path.
Measured over all 44 runs on disk:

| | closest genome searched serially |
|---|---|
| legacy absolute binning | **7 / 44** |
| rank binning | **44 / 44** |

The melittin ground-truth run goes from a single undifferentiated wave of 5 to `[1, 2, 2]`
— a real closest-first gradient. Regressions: `tests/test_adjudication_fixes.py`.

**⚠️ Measurement attempted 2026-08-21 — inconclusive, and instructive about why.**
A/B on LRZ (`lrz-cpu`, jobs 5757307 / 5757308): melittin `P01501`, easy mode, identical
inputs (arm B resumed arm A's work dir, so the fetch was the *same cached task*), one
variable — `--rank_wave_binning true` vs `false`.

**Both arms produced identical waves** `[1, 2, 1]` at normalised distances
`0.200 / 0.700 / 1.000`, and byte-identical results (3 HIGH / 13 MEDIUM / 145 LOW).

The reason is the target set, not the code. Easy mode auto-picked four genomes sitting at
three *distinct* taxonomic ranks — *Apis cerana* (genus), *Bombus turneri* + *Bombus
ignitus* (family), *Agapostemon virescens* (order) — so the distances are genuinely spread
and **legacy binning separates them correctly too**. The degenerate single-wave collapse
needs a *uniform-depth* target set, which is what the pinned GT genome set happens to be.

So: this run shows the fix does not regress, and nothing more. It does **not** validate the
wavefront. Note the expanding DB did do work within the run — wave 1 added 0 new genes,
wave 2 added 6, wave 3 added 1 — but whether those seeds changed any final call is still
untested, because both arms had the same wave structure.

**✅ MEASURED 2026-08-31 (job 5768082) — clean negative result. The wavefront changes nothing.**

Three arms on the pinned pro-mode melittin GT set, one variable, separate work dirs so no
cached task could mask a difference. The third arm is the one the 2026-08-21 attempt lacked:
an identical **replicate** of arm A, giving a noise floor, because a difference between two
arms is uninterpretable while the run's determinism is unknown.

| arm | setting | waves | DB growth between waves |
|---|---|---|---|
| A | wavefront ON | **3 graded** `[1, 2, 2]` | 1 → 4 → 6 genes |
| A′ | ON, replicate | 3 graded `[1, 2, 2]` | 1 → 4 → 6 genes |
| B | `--disable_wavefront true` | **1 wave of 5** | none — initial DB only |

The control genuinely controlled this time — that was the specific failure of the 2026-08-21
run, where both arms produced identical waves and therefore tested nothing. Here arm B searched
every genome against the initial database only, while arm A searched genomes 2–5 against a
database augmented by nearer relatives' recovered orthologs.

**All three arms produced byte-identical results**: 0 HIGH / 7 MEDIUM / 13 LOW, coordinate
distance 0 for every pairing. `diff(A,B) = diff(A,A′) = 0`.

So: **the expanding wavefront made no difference to any call on this target set**, and the run
is deterministic here (worth knowing in itself — the suspected non-determinism did not appear).

**What this does and does not license.** It licenses removing the wavefront from the *melittin
result's* explanation: that result is produced by neighbourhood restriction, and the mechanism
section must say so. It does **not** license "the wavefront is useless", and the negative must
not be over-read — **the experiment did not stress the mechanism.**

The wavefront can only pay off through a **stepping-stone**: a target whose gene is too divergent
to hit from the original query, but reachable from an ortholog recovered in a *nearer* species and
added to the database. That needs a genuine divergence ladder. This set has no such ladder — the
normalised distances were

```
wave 1: dist 0.667   wave 2: dist 0.667   wave 3: dist 1.000
```

i.e. **two distinct values across five genomes**, all bees comfortably within reach of the *Apis*
query (every target was recovered in wave 1's database anyway). There was no target that needed a
stepping stone, so no stepping stone could be observed. A flat gradient cannot test a gradient
mechanism.

**⬜ The experiment that would actually test it** — a divergence ladder, targets at increasing
distance where the farthest is beyond direct reach of the query but within reach of an
intermediate's recovered ortholog. The oskar ladder set (mosquito → Nasonia → Gryllus) is the
natural candidate: chosen for exactly that spread, and already validated on the cluster. Run
`--disable_wavefront` A/B there before concluding anything general about the mechanism.

Honest summary for the paper: *the wavefront was measured on the melittin benchmark and changed
nothing there, so the melittin result is attributable to neighbourhood restriction alone. Whether
it helps across a real divergence gradient is untested.*

## F5 Two incompatible definitions of collinearity — ✅ FIXED 2026-07-26

`cluster_grs.score_flexible_synteny` used an **adjacent-pair monotone fraction**;
`iterative_search_runner._longest_collinear_run` used a proper **LIS**. The same concept
scored differently in ranking than in seeding, so a block could be judged collinear when
it was seeded and non-collinear when it was ranked. The code comment in `cluster_grs`
even said "LIS could be used".

**Fix.** The LIS moved to `bin/sequence_utils.py:longest_collinear_run` — the repo's
single-source-of-truth module — and both callers now use it (`iterative_search_runner`
keeps a thin alias). Adjacent-pair counting was also fragile: one transposed gene inside
an otherwise perfect run of 5 breaks two pairs (0.75) but costs an LIS one element
(0.80). Regressions: `tests/test_adjudication_fixes.py`.

## F6 ⭐⭐ `flanking_support` is per-block, not per-call — the overcalling mechanism — ✅ RESOLVED

Every GOI feature emitted inside a block is classified with
`flanking_support = block_flanking_support` (`:3373, 3435, 3457, 3540, 3605`). A junk hit
anywhere in a 5-flanking block inherits `S_block = 5` and, via the un-floored
`fallback_hit_span` clause (A.4), reaches **MEDIUM at any identity** provided
`qcov ≥ 0.25`.

This is the mechanism behind TODO §1v (yeast STE2: 32 MEDIUM at 21–27 %; Drosophila
Defensin: 170 GOI calls for 3 real orthologs).

**✅ Done — the identity floor.** New
`--classify_fallback_strong_min_identity_floor` (default **30.0**, `0` = legacy) is AND-ed
into the strong-flanking clause. Calibrated against the two real populations rather than
picked: the genuine *Tetramorium* melittin rescue sits at **34.7 %**
(`BlockFlankingSupport=7`, `QueryCoverage=0.686`); the yeast false positives run
**21–27 %**. 30 separates them with margin on both sides.

Measured impact over **671 GFFs** from every run on disk — 1 436 MEDIUM calls rest on this
clause:

| | count | identity range |
|---|---|---|
| kept | 953 (66 %) | ≥ 30 % |
| **demoted to LOW** | **483 (34 %)** | 17–29 % |

The melittin ground-truth runs (`melettin_gt_v4 … v7`) lose **zero** calls, and the §1y
coordinate-exact ant melittin hit at `15,634,953` / 34.7 % is **kept** in both
`mel_det_rep2` and `melittin_full`. The yeast runs lose their junk. Regressions:
`tests/test_adjudication_fixes.py`.

**🛑 Per-call flanking support — investigated and DELIBERATELY NOT IMPLEMENTED.**

`flanking_support` is still the block's count, inherited identically by every call inside
it, and the obvious next step was to compute support in a window around each call. Before
writing it I measured what it would actually buy, over the 1 436 affected calls in every
`plot_inputs_*` GFF on disk. It buys nothing, and it would cost real orthologs:

1. **The motivating case is already solved by the identity floor.** Of the 129 yeast STE2
   false positives on this clause, **0 survive** it. Of 115 melittin calls, **97 (84 %)
   do**, at median identity 50 % and median 6 local flanking genes. There is no residual
   yeast population for a local gate to catch.
2. **Local flanking does not separate the two populations anyway.** Within 100 kb the
   yeast false positives have a median of 2 flanking genes and only 6 % have none — they
   are not isolated hits, so no threshold cleanly divides them from melittin's 6.
3. **It would demote real orthologs.** Of the 953 calls that survive the identity floor,
   36 have zero flanking within 100 kb — and they are the *highest-identity* calls in the
   set (86.5 %, 85.9 %, 72.9 %, 71.3 %; median 38.4 %), concentrated in `Y-e3_evolution_v2`
   (MRJP/Yellow tandem family) and `luciferase_rerun_auto`. A gene recovered at 86 %
   identity in a rearranged neighbourhood is precisely what SynVoy exists to find. A
   local-flanking gate would demote it.

Adding an unvalidated threshold that fixes nothing measurable and demotes the tool's best
results is the exact pattern flagged in Part 0 — *"parameters set by an agent when coded
and never changed"*. **Closing this as won't-do**, with the numbers, rather than leaving it
open as a plausible-sounding task. Reopen only if a case appears that the identity floor
misses.

Corroborating evidence for the knock-on: **109 of 130 region rows have
`goi_overlap=True`** — junk GOI calls are scattered so widely that the GOI-overlap
prioritisation in `cluster_grs` carries almost no information.

## F8 ⭐⭐⭐ Easy mode searched for the *flanking* genes and reported them as the GOI — ✅ FIXED 2026-08-21

Found by the LRZ run above, not by any test. [`main.nf:1431`](../main.nf#L1431) read:

```groovy
rescue_query_ch = ( params.query
    ? channel.value(file(params.query))
    : EXTRACT_FLANKING.out.faa.map { rec -> rec[1] }.first() )
```

In **easy mode `params.query` is null** — the user passes `--query_id` — so the fallback
fires and the channel carries `flanking_proteins_locus_<n>.faa`. Three processes then run
against the flanking proteins instead of the gene of interest: `RESCUE_GOI_HULL` (§17),
`RESCUE_STRONG_SYNTENY` (§1e Phase B), and `RECIPROCAL_BEST_PARALOG` (§1j).

**What it produced.** All three HIGH calls in job 5757307 were 32 bp (~10 aa) fragments of
`gene-LOC107964339_exon_1` — a *flanking* gene — at "100 % identity", class
`synteny_hull_rescue`, one per genome. The report's headline, *"3 high-confidence GOI
ortholog annotation(s)"*, was entirely this artefact. Melittin at 100 % in two *Bombus*
genomes is the tell; that is not biology. The plausible real ortholog (*Apis cerana*
`NC_083855.1:12,201,524`, ~134 aa, 100 %) was ranked only **MEDIUM**.

The rescue log states it plainly: `13 HIGH flanking, no HIGH GOI -> rescued model
id=100.0% conf=HIGH`.

**Why nothing caught it.** `ASSIGN_LOCUS_OWNERSHIP` (§1m), which exists to demote exactly
this class of mislabel, consumes `paralog_inputs_ch` — the ITERATIVE_SEARCH region files.
Rescue GFFs are mixed in later, at `STAGE_REGION_GFF`, so **rescue models bypass ownership
entirely**. `self_consistency` *did* flag 6 HIGH rows `identity_coverage_decoupled`
(query coverage 0.21–0.29) — the detector works; the headline ignores it. This is the
two-layer thesis in miniature: the search layer found the neighbourhood correctly, and the
adjudication layer named the wrong gene in it.

**Why it survived to now.** Every §1m/§17 validation was done in **pro mode** (decorin/DCN),
where `params.query` is set and the branch is correct. Easy mode is the path README and
QUICKSTART document, and the one the BA students run.

**✅ Fixed** — `rescue_query_ch = normalized_gene_ready_ch.first()`. `NORMALIZE_QUERY.out.fasta`
is the real GOI protein in both modes. This also repairs pro mode, where the old branch
passed the *raw* `--query` file — DNA if the user supplied DNA, which miniprot cannot use.
The edit is 3 lines → 1, so the `workflow {}` body shrinks (safe for the JVM method-size
limit). 550 tests pass; re-run verification in progress.

**Both follow-ups closed 2026-08-22:**

1. **A 10 aa 100 % match can no longer reach the headline**, whichever query produced it.
   `apply_coverage_demotion` relabels a known-low-coverage, high-identity call
   `identity_coverage_decoupled` and excludes it from the HIGH/MEDIUM counts — turning the
   existing QW3 detector from an annotation into a demotion. Only fires on *known* low
   coverage, so `rescue_goi_hull.py` now emits `QueryCoverage` (it computed it as a gate all
   along and threw it away). See A.4.
2. **Rescue models now reach the ownership check.** The rescue emits the model's protein and
   the sub-workflow feeds it into `ASSIGN_LOCUS_OWNERSHIP` as its own task. See A.7.

Both landed behind the `main.nf` sub-workflow extraction described below — the wiring in (2)
was impossible before it.

**Follow-up (2) verified against the real v2_decorin run — 2026-08-30, and it is a no-op there.**
The wiring alone was not enough: `build_locus_ownership` indexed ownership rows by the raw
`genome` column, and the rescue task writes that column as `<genome>.hull_rescue` (the tag
exists so its output filename cannot collide with the seeded task's). The dedup records on
the other side of the lookup are already canonicalized, so the key never matched and every
rescue-derived call still skipped the RBH check. Fixed by canonicalizing the row's genome
(`generate_report.py`, committed in `b3fcfe6`); regression
`tests/test_locus_ownership.py::test_rescue_tagged_genome_still_matches`.

Replaying **both** code versions over the staged ownership TSVs of job 5760122 (25 files,
792 rows, 6 of them rescue-tagged) gives **byte-identical verdicts** — 15 records evaluated,
11 `paralog_misassignment`, same classes, same owners. The reason is worth recording, because
it is the difference between "the guard was broken" and "the guard was redundant here":

> All 6 rescue rows are **shadowed** by a seeded row at the same coordinates carrying the
> same verdict. For the cow decorin at `NC_037332.1:21,015,210` the rescue row scores
> `locus_1|gene-DCN` at bit 1713 and the seeded row scores the same gene at bit 886 — the
> rescue row is the better alignment, but it is not the *deciding* one. §18's collinearity
> bridging already recovers that gene inside the seeded block, so the seeded ownership task
> models the same span independently.

**Consequence: the 2026-08-24 decorin headline stands** — the 4 HIGH decorins were RBH-checked
after all, via their seeded twins, and the 2 biglycans were correctly demoted. The fix is
defensive rather than corrective: it matters for a call the rescue path finds *alone*, which
on this dataset does not occur (0 of 6 orphans). It should not be described as having changed
any published number. Verify on new data with `scripts/lrz/verify_ownership_fix.py`.

## F7 Documentation-layer defects (fixed in the 2026-07-26 docs pass)

`--keep_intermediate` is declared and read by nothing; `--no_anchor_grid` is a
`plot_synteny.py` switch that the Nextflow module never forwards; `--help` does not exist
though `main.nf:714` still advertises it; 59 parameters were documented nowhere. All now
corrected in the docs — but `main.nf:714` and the two dead parameters are still live in
the code.

*Update 2026-10-05:* the validation message no longer advertises `--help`; every plot
switch is reachable through `--plot_extra_args`; and the docs now mark **five** declared
parameters as inert (`keep_intermediate`, `max_retries`, `multi_profile`,
`multi_profile_max_jobs`, `locus_ownership_synteny_window`). A static test
(`tests/test_nextflow_param_wiring.py`) now fails when a module passes a flag its script
lacks or reads a parameter no config declares.

## F9 ⭐⭐⭐ Over-calling is a MISSING-EVIDENCE problem, not a threshold problem — ⚠️ measured 2026-08-31

The benchmark re-run (D.1b) called PRESENT in **16 of 19 species and ABSENT in none**, with the
same 5 false positives as May. Every fix since then has adjusted *thresholds*. This section
records why that could never have worked, using the per-call evidence dumped from the run.

**The evidence does not separate the two populations.** HIGH/MEDIUM calls, grouped by whether
the species actually has the gene:

| | identity | query cov | flanking support |
|---|---|---|---|
| gene **PRESENT** | 84.3 · 73.5 · 63.5 · 59.8 · 54.5 · **30.2** | 1.00–0.27 | 6–10 |
| gene **LOST/ABSENT** | **93.8** · 57.7 · 37.0 · 31.9 · 30.4 · 30.2 | 0.56–0.26 | 7–8 |

The **highest-identity call in the entire benchmark (93.8 %, *Bombus terrestris*) is a false
positive**, and the lowest true call (*Formica*, 30.2 %) ties the lowest false one. Nearest
neighbours across the boundary — *Colletes* 54.5 %/0.427/4 exons (real) versus *Bombus impatiens*
57.7 %/0.443/4 exons (lost) — are indistinguishable on every feature the classifier has.

**Why flanking support carries no signal here, structurally.** It is 6–10 on *both* sides,
because a gene lost from a conserved locus leaves the neighbourhood intact — that is precisely
what gene loss looks like. Synteny tells you *where to look*; using it as evidence that *what you
found is the gene* is a category error, and it is the mechanism behind the over-calling.

**The RBH check cannot rescue it either, and the reason is fixable.** Best-home-match is
`gene-Melt` for 14/22 PRESENT rows and 15/23 LOST/ABSENT rows — no discrimination. The cause is
in the last column of every ownership TSV: **`n_paralogs_compared = 2`**. The home panel holds
only `gene-Melt` and `gene-LOC726866`. *Bombus terrestris* RBHs to melittin at bit 89 / 53.3 % —
**stronger than real *Colletes* (78)** — because **bombolitin, the melittin-family gene that
replaced melittin in *Bombus*, is not in the panel for it to match instead.** §1m's
`find_family_paralogs` recovered BGN/ASPN for decorin but only one extra member here.

### The AMBIGUOUS tier (shipped 2026-08-31, default ON)

A fourth verdict for a candidate promoted **only** by flanking support:
`AMBIGUOUS` / `syntenic_candidate_unconfirmed`. It stays in the JSON and is drawn in the plots
(palest, most broken outline in the legend) but is **excluded from the ortholog counts** and
named separately in the headline. Kill switch `--disable_ambiguous_tier` *(until 2026-10-05
this existed in the script only; the pipeline did not pass it, so it had no effect on a run)*.

**Measured cost, replaying all 53 real benchmark calls through both classifiers:**

| | TP | FP | FN | TN | precision | recall | F1 |
|---|---|---|---|---|---|---|---|
| MEDIUM (pre-2026-08-31) | 6 | 4 | 2 | 3 | 0.600 | 0.750 | 0.667 |
| AMBIGUOUS tier | 2 | 1 | 6 | 6 | 0.667 | **0.250** | 0.364 |

11 of 13 HIGH/MEDIUM calls become AMBIGUOUS. It removes 3 false positives (*Bombus impatiens*,
*Melipona*, *Solenopsis*) and costs 4 true positives (*Cardiocondyla*, *Colletes*, *Euglossa*,
*Formica*). **That is a near 1:1 trade, and it is the expected result**: the tier adds no
information, it only stops the tool asserting what it cannot support. Precision barely moves
because the errors were never separable.

So the tier is an **honesty fix, not an accuracy fix**, and it should be argued as one. Part of
the F1 drop is an artefact of a metric that rewards guessing PRESENT — there are more PRESENT
species than LOST ones in the truth set.

### The "populate the paralog panel" fix — PROPOSED, TESTED, REFUTED (2026-08-31)

The obvious next move was to fill the 2-member panel with the melittin gene family so the RBH
check could reject the *Bombus* call. **It was tested against the real *Apis mellifera* proteome
(9 935 proteins) and it does not work.** Recorded because the negative is more useful than the
proposal:

| home gene | len | SW score | score / self |
|---|---|---|---|
| `gene-Melt` (melittin) | 70 | 290 | **0.797** |
| `gene-LOC726866` | 282 | 26 | 0.071 |
| `gene-Apamin` | 46 | 21 | 0.058 |
| `gene-Mcdp` | 50 | 16 | 0.044 |

Apamin and MCD peptide *are* annotated in the proteome and *are* bee-venom relatives, but against
a self-score of 364 they score 21 and 16. They share only the venom-gland signal peptide; the
mature peptides are unrelated (apamin is a disulfide neurotoxin, melittin an amphipathic lytic
peptide). **They are a functional family, not a sequence family** — melittin is effectively
single-copy in *Apis*, so there is nothing to put in the panel.

**The real reason RBH cannot fix this is architectural.** The *Bombus terrestris* call matches
`gene-Melt` at bit 89 / 53.3 % — **stronger than real *Colletes* (78 / 37.1 %)**. SynVoy is
finding a genuine melittin-like sequence in the *Bombus* syntenic locus; the literature says
melittin is absent from *Bombus* venom and **bombolitin** replaced it, and bombolitin is not
annotated in the *Bombus* proteome at all. So the recovered sequence is a **target-lineage-specific
paralog**, and a home-proteome RBH can only reject a call that matches a *home* paralog better.
There is no *Apis* gene for bombolitin to match. **No amount of home-panel curation reaches this
class of error.**

**This strengthens the case for AMBIGUOUS rather than weakening it.** "Melittin ortholog vs
bombolitin paralog" is not decidable from the home genome and the target sequence alone, so
AMBIGUOUS is the *correct* verdict here, not a stopgap for a missing threshold.

### Architecture profiling — PROTOTYPED, test INCONCLUSIVE (2026-08-31)

Third hypothesis: the melittin precursor is signal peptide -> acidic propeptide -> amphipathic
cationic mature, and the 2026-07-19 six-frame audit used exactly that to sort 11 real family
members from ~13 overcalls by hand. Prototyped as a **query-derived** profile (6-bin Kyte-Doolittle
hydropathy + 6-bin net charge, compared by normalised L1) rather than a hardcoded melittin rule,
so it would generalise to any query. The query's own signature comes out clearly:

```
hydropathy  +0.41 +0.99 -1.68 -1.05 +1.31 +0.06      (hydrophobic signal, then ...)
charge      +0.00 -0.08 -0.42 -0.42 +0.08 +0.25      (... acidic propeptide, cationic mature)
```

Scored against the best melittin-like protein in each species' proteome it did **not** separate:
`min(PRESENT) = 0.828 < max(LOST) = 0.898`.

**But that test is inconclusive, not a refutation, and the reason matters.** The test bed is
broken: *Apis cerana*'s best proteome match scores SW 49, while its recovered ortholog scores
**342**. The truth table already says why — *"NCBI annotation gap (RefSeq does not currently
annotate melittin)"*. For most of these species the proteome does not contain melittin at all, so
"best melittin-like protein" is noise, and a discriminator cannot be evaluated against noise.

**⬜ To actually test it**, score the profile against SynVoy's *recovered models* (the 12
HIGH/MEDIUM calls of D.1b), not against proteome proxies. That needs the models translated from
the staged genomes — a small batch job, not login-node work. Two design notes for whoever does it:
the coarse 6-bin profile is probably too permissive (many proteins are hydrophobic-N /
acidic-middle / basic-C), and melittin's defining property is the **amphipathicity** of the mature
helix, which a linear profile cannot see — a helical-wheel moment is the feature with real
discriminating power.

**⬜ Other routes**, in increasing order of cost:
1. **Phylogenetic placement.** A true ortholog groups with the query in a species-tree-consistent
   way; a lineage-specific paralog does not. SynVoy already builds trees (`compute_tree.py`) but
   the tree never feeds back into the call (§1p — verdicts never reach the tree or plot). This is
   the route that could reach the *Bombus* class of error, because it does not depend on the home
   genome containing a paralog to match against.
2. **Target-side annotation / synteny of the paralog itself** — expensive, and out of scope.

### Summary: three hypotheses, none cheap

| # | hypothesis | outcome |
|---|---|---|
| 1 | tune identity / coverage / flanking thresholds | **refuted** — distributions overlap completely; the highest-identity call in the benchmark is a false positive |
| 2 | populate the home paralog panel | **refuted** — melittin has no sequence family in *Apis*; the confounder is a *target*-lineage paralog a home RBH can never see |
| 3 | architecture profiling | **inconclusive** — plausible, but untestable against proteomes that do not annotate the gene |

The over-calling is **not reachable by reweighting evidence the pipeline already has.** That is
the finding, and it is what justifies the AMBIGUOUS tier as the shipped response rather than a
placeholder: where the evidence cannot decide, the honest output is a named candidate with
orthology explicitly unasserted.

Until one exists, the defensible output for this class is exactly what the tier now emits: a
named candidate in a conserved neighbourhood, with orthology explicitly not asserted.

---

# Part C — Known open items (from `docs/TODO.md`, re-verified)

| § | Item | Severity | State |
|---|---|---|---|
| **§1x** | **Silent discard** — a confident hit on the correct gene is rejected by the synteny gate and never mentioned. `demo_Def`: *Anopheles* defensin found at 54–72 % inside the real RefSeq gene, reported as "no ortholog". | ⭐⭐ | ✅ **Fixed 2026-07-26** — report block, per-rejection logging, and the plot no longer claims absence. Off-block marker ⬜ (needs a `main.nf` channel) |
| **§1y** | **Coordinate regression** — the ant melittin was hit in 8 runs; the newest run is 2 543 bp off and nobody noticed, because scoring is species-presence not coordinate overlap | ⭐⭐ | ✅ **Scorer + 28-model truth set + CI freeze built 2026-08-30**, and on its first real use (D.1b) it caught a live overclaim: the benchmark's only recall gain is a model **63,936 bp** off the gene. ⬜ Cardiocondyla accession mapping (13 models unscoreable) |
| §1v | Classifier ignores the synteny evidence the run computed | ⭐⭐ | ✅ **Fixed 2026-07-26** — identity floor; per-call support measured and closed won't-do (see F6) |
| §1u | `goi_block_flanking` drives ranking but is not emitted | ⭐ | ✅ **Fixed 2026-07-26** — emitted in the scores TSV alongside `synteny_score` |
| §1z | `QueryCoverage` read as target recovery; target coverage never reported | 🐛 | Open |
| §1p | Post-processing verdicts never reach the plot or the tree | ⭐ | Open. Partly eased 2026-10-05: rescue-pass gene models are now written to `rescue/locus_<N>/` (before, a hull-rescued ortholog existed in the output only as a span in the report). The plots and `homology.tsv` still do not show rescue models or paralog verdicts. |
| §1p.1 | Target exon composition not passed to plot → synthetic uniform exons drawn | ⭐ | ✅ **Fixed 2026-09-17** at the source: the search had been collapsing every in-block flanking model to one CDS row. Verified on the cluster (`fix7`): 100 % of multi-exon flanking models now carry real CDS rows (was 1–3 %). Old run outputs still hold collapsed rows. |
| §1t | Every distance parameter is metazoan-scale (found by pointing it at an 11 Mb yeast) | ⭐ | Open |
| §1w | Duplicate regions with byte-identical spans | 🐛 | Open |
| §1q, §1s | `strong_synteny` read before assignment; numeric CLI params crash the run | ✅ fixed | Committed 2026-08-21 (`49df417`) |
| — | **Ivan's annotation agent deletes coding sequence** on validation failure, silently, marked only `class=modified` | ⚠️ external | Blocking the next annotation round |

Counted: **38 unchecked action items** across Part 0 + Part 1 of `docs/TODO.md`.

---

# Part D — Benchmark

## D.1 Current numbers (`benchmark_results/`, last scored **2026-05-03**)

| tool | TP | FP | FN | TN | precision | recall | F1 |
|---|---|---|---|---|---|---|---|
| **synvoy** | 7 | 5 | 1 | 2 | 0.583 | **0.875** | **0.700** |
| tblastn | 2 | 0 | 6 | 7 | 1.000 | 0.250 | 0.400 |
| miniprot | 2 | 0 | 6 | 7 | 1.000 | 0.250 | 0.400 |
| mmseqs_genome | 2 | 0 | 6 | 7 | 1.000 | 0.250 | 0.400 |

The shape of the result is exactly the story: **SynVoy trades precision for recall.** The
5 FP are the overcalling problem (F6) and the 0.875 recall is the contribution. Per tier,
SynVoy is the only tool with any tier-2/tier-3 recall at all (baselines: 0.000).

### D.1b Re-run on current code — 2026-08-31 (job 5768083, 19 targets)

| scoring | TP | FP | FN | TN | precision | recall | F1 |
|---|---|---|---|---|---|---|---|
| **2026-05-03 (D.1)** | 7 | 5 | 1 | 2 | 0.583 | 0.875 | 0.700 |
| **2026-08-31, current code** | **8** | 5 | **0** | 2 | 0.615 | **1.000** | **0.762** |

On species presence the tool improved: recall 0.875 → **1.000**, F1 0.700 → 0.762. The entire
gain is one species — *Tetramorium bicarinatum* moving from FN to TP.

**That single improvement is an artefact, and the coordinate scorer catches it.** Tetramorium is
called PRESENT at HIGH confidence over `OV788322.1:15,483,596–15,743,000` — a **259 kb** region
that does contain the true melittin at 15,634,953. But the actual GOI *model* is placed
**63,936 bp** away from the gene. In `mel_det_rep2` the same gene was modelled at offset **0**,
132 bp overlapping. So between those runs the **gene model got worse while the benchmark score
got better** — the §1y thesis demonstrated end-to-end on the flagship benchmark in a single run.

> **Erratum 2026-09-10: rescore before citing anything in this section.** Until 2026-09-10
> `score_coordinates.py` read only `mRNA` GOI features and silently skipped tandem copies,
> which SynVoy writes as `gene` features. On a local rerun of this benchmark, that hid **42 of
> 95** GOI calls (6 of them MEDIUM). The fix moves locus_recall **0.04 → 0.28** on the same
> output, and it turned LRZ det_rep R2's "*Tetramorium* melittin absent" into a hit at the
> correct locus. The "model placed 63,936 bp away" claim below was scored the old way, so a
> tandem copy on the true gene would have been invisible to it. Its run outputs are on LRZ;
> rescore them with the fixed scorer.

Coordinate scoring of the same run:

```
locus recall      4/12  searched truth models overlapped by a call
substantial       1/9   curated models recovered over >=50% of their span
false positives   1     call on a locus where the gene is LOST
unscoreable      13     curated on a different assembly (Cardiocondyla)
```

Two further findings from the same run:

- **SynVoy called PRESENT in 16 of 19 species and ABSENT in none.** All four species whose truth
  is `LOST` (*Bombus terrestris*, *B. impatiens*, *Melipona*, *Tetragonula*) got a MEDIUM call,
  as did *Solenopsis* (truth `ABSENT`). The 5 FP are unchanged from May — **none of the
  2026-07/08 adjudication work moved precision at all.** The identity floor and coverage demotion
  fixed the *labels*; they did not stop the tool asserting presence nearly everywhere.
- **13 of 28 truth models could not be scored**, because the curated Cardiocondyla coordinates
  are on `GCA_019399895.2` (`CM079762.1`) while the run used `GCF_019399895.1` (`NC_091864.1`).
  These are the same linkage groups in different accession namespaces; an accession mapping would
  make them scoreable. Worth doing — Cardiocondyla has the highest GR2 copy count in the paper,
  its call sits ~70 kb from the curated array, and it is currently counted as a presence TP that
  has never been checked by coordinate.

⚠️ `median_overlap_frac = 1.0` in that output is a median over **n = 1** curated hit and should
not be read as a typical recovery. The statistic needs a minimum-n guard.

Three caveats that must be stated if this is published:

1. **It is stale.** Scored 2026-05-03, before §1e, §1m, §1n, §18 and every 2026-07 fix.
2. **n is tiny** — 15 scored species, 8 positives.
3. **It scores species presence, not coordinates.** §1y shows this is precisely what
   allowed the melittin overclaim. Re-score by coordinate overlap before publishing.

## D.2 TOGA — a complete non-result

`benchmark_results/toga/calls.tsv`: **7 rows, all `ABSENT`, all `extra=toga_no_output`.**
TOGA does not appear in the confusion table at all. It has produced zero data.

Failure point, from `benchmark_results/toga/work/chains/_run_Cardiocondyla_obscurior/run.log`:
the pipeline reaches `### Lastz Alignment Step ###`, dispatches 18 jobs through
make_lastz_chains' own Nextflow layer, and the log ends there. Two known blockers:

- **RAM.** The machine has 15.3 GB; TOGA's chain step needs ≥ 22 GB. Logged as §3a.
- **Soft-masking.** make_lastz_chains 2.0.8's py2bit path strips soft-masking, which
  overflows LASTZ's 32-bit scoring at sensitive K. Workaround already identified:
  `twoBitToFa` + `--seq1_chunk`.

**Recommendation: stop trying to run TOGA locally.** Two viable routes, in order of
effort:

1. **Run it on the LRZ Linux-Cluster (CoolMUC).** It is a CPU job, it is exactly what the
   cluster is for, and the RAM blocker disappears. This is the cheapest path to a real
   number.
2. **Replace it.** TOGA is a whole-genome annotation projector and is arguably the wrong
   comparator for a single-gene locator anyway — a reviewer could reasonably call it an
   unfair matchup in both directions. **GENESPACE** or **MCScanX** (already wired,
   `run_mcscanx.sh`) is the honest synteny-based comparator; the archived
   OrthoFinder/SonicParanoid wrappers were dropped for annotation-coverage reasons that
   are worth stating explicitly in the methods rather than leaving as a silent omission.

If TOGA has to be in the table, run it on CoolMUC for the three tier-3 species only. If it
cannot be run, **say in the paper that it was attempted and why it was excluded** — an
explained omission is fine; a silent one is not.

---

# Part E — Infrastructure

## E.1 Cluster (LRZ)

Reconstructed from `lrz_backup_2026-07-21/CLUSTER_RUN_ANALYSIS.md`. **Nothing has run on
the cluster since 2026-07-08.** All run output directories were deleted; only logs
survive. Ten jobs, of which five were bring-up failures (ESMFold CUDA-OOM on the 16 GB
V100) and four succeeded after the fix.

Strongest surviving evidence — job `5702387`, the two-phase GPU rework:

```
plm:        {requested: true, embedded: 247, failed: 0}
structural: {orf_candidates: 215, orf_folded: 40, folded_items: 72, failed_items: 0}
```

`failed_items: 0` against a previous `1501/1501` failure rate, with TM-scores spread
0.153 → 0.948. The GPU layer is genuinely working now, and Foldseek's earlier silent
no-op (CUDA-in-fork inside `ProcessPoolExecutor`) is fixed.

**Assessment: the cluster is a solved problem parked in a good state, not an open risk.**
Leave it. The one thing worth doing is running TOGA on the CPU side (D.2), because that
unblocks the benchmark rather than adding new capability.

## E.2 Website

**Nothing exists.** No `mkdocs.yml`, no `_config.yml`, no `docs/index.md`, no Pages
workflow — the only CI is `.github/workflows/test.yml`.

The good news is that the docs pass on 2026-07-26 left `README` / `QUICKSTART` / `USAGE` /
`PARAMETERS` / `OUTPUT` accurate and cross-linked, which is 90 % of a documentation site's
content. Turning that into a site is a few hours of MkDocs Material configuration, not a
project. **Do it last, and do it from the existing markdown** — it needs no new writing.

---

# Part F — What to do, in order

The organising principle: **nothing in the adjudication layer should be published before
it is fixed, and the mechanism claim must be corrected before the paper is written.**

### Immediately (hours) — stop the bleeding

1. ✅ **Committed 2026-08-21** (`49df417`). The adjudication fixes are on `origin/dev`.
2. ✅ **Anchor-grid identity labels — fixed 2026-07-26.** The GOI column was routed through `draw_goi_array_grid`, which bypasses `draw_arrow` — the only code that draws identity labels. So every flanking column carried a number and the one column a reader cares about most carried none, while the legend still said "number = % identity". Labels restored on each enrolled copy (both the grid and the track-plot variant). The obsolete `×2` badge test was replaced: the copy count is now conveyed structurally by drawing one arrow per copy, so the test asserts both copies appear with their own identity instead.
3. ✅ **Committed 2026-08-21** (`8e44d55`).

### Before any figure or number is published (days) — the adjudication layer

4. ✅ **F1 — p-value fixed 2026-07-26.** The tested statistic is now `synteny_score` (bonus- and rescue-free) and is emitted for auditing. Verified invariant to GOI overlap.
5. ✅ **F6 + §1v — resolved 2026-07-26.** Identity floor added and calibrated (483/1 436 MEDIUM calls demoted; 0/129 yeast false positives survive, 97/115 melittin calls do). Per-call flanking support was measured and **closed as won't-do** — it fixes nothing the floor misses and would demote 36 high-identity orthologs in rearranged neighbourhoods. See F6 for the numbers.
6. ✅ **§1u — `goi_block_flanking` now emitted** in `regions/*.scores.tsv`.
7. ✅ **§1x — fixed 2026-07-26.** `select_dispersed_goi_seeds` now records every hit that clears the quality bar and is then refused by a synteny gate, with the gate that refused it and the distance to the nearest flanking anchor, and logs one INFO line per rejection. Written to `hits/<genome>.nonsyntenic.tsv` — which `generate_report.nf` already stages wholesale, so this needed **no new channel and no `main.nf` edit** (that method body is at the JVM 64 kB limit). The report gains a `rejected_candidates` block, and `goi_absent_genomes` is split from the new `goi_found_but_not_syntenic` — opposite claims that used to share a field. Opt out with `--report_nonsyntenic_candidates false`. ⚠️ **Plot, partially:** the anchor-grid legend said **"no ortholog"** for an empty cell — a claim the figure is not entitled to make, and precisely the misreading §1x exists to prevent. Now reads **"not placed here"**, with a caption line stating that an empty cell means no ortholog was *placed* in this neighbourhood, not that the gene is absent, and pointing at `synvoy_report.json → rejected_candidates`. The track plot's per-row "✗ GOI not found" is now "✗ GOI not placed". ⬜ Drawing the specific rejected candidate as an off-block marker still needs the TSV to reach `PLOT_SYNTENY`, which needs a new channel in `main.nf` — blocked on the sub-workflow refactor (CLAUDE.md §17: that method body is at the JVM 64 kB limit, where a 4-line comment has already re-triggered `UTF8 string too large`).
8. ✅ **F2 — fixed 2026-07-26**, and it was worse than reported: the legacy measure scored a *scrambled* neighbourhood 1.000 and a *conserved* one 0.500. Now scored against the home strand, inversion-tolerant. ⬜ **F3 left as a description matter** — see the section above for why rescaling is not worth doing on its own.

8b. ✅ **F8 — easy mode was searching for the flanking genes and reporting them as the GOI. Fixed 2026-08-21** (`main.nf:1431`). This invalidated the headline of *every easy-mode run*: the three "high-confidence orthologs" on the melittin LRZ run were 10 aa fragments of a flanking gene at 100 %. ✅ Both follow-ups closed 2026-08-22: the QW3 `identity_coverage_decoupled` detector now **demotes** instead of merely annotating (only on *known*-low coverage, so `rescue_goi_hull.py` now emits `QueryCoverage`), and rescue models are routed into `ASSIGN_LOCUS_OWNERSHIP` via a new `--output_faa`. See F8.

8c. ✅ **`main.nf` sub-workflow extraction — the JVM method-size blocker is gone (2026-08-22).** The adjudication + staging + report block (two rescue passes, reciprocal-best paralog, locus ownership, per-locus staging, `GENERATE_REPORT`) moved to `subworkflows/adjudicate_and_report.nf`. The `workflow {}` body went 56 693 → 49 975 bytes, and — the actual point — the extracted block now has its own method budget, so new adjudication steps no longer risk `UTF8 string too large`. **New steps of that kind belong in the sub-workflow, not in `main.nf`.**

### Before the paper's mechanism section (days)

9. ⚠️ **F4 — the measurement is prepared and ready to submit (2026-08-30).** `scripts/lrz/f4_wavefront.sbatch` + `f4_wavefront_payload.sh` run **three** pro-mode arms on the pinned melittin GT set: wavefront ON, an identical ON **replicate**, and `--disable_wavefront true`. The replicate is not padding — SynVoy has known run-to-run non-determinism, so `diff(A,B)` is uninterpretable without the noise floor `diff(A,A')`; the rule is that a wavefront effect exists only if `diff(A,B) > diff(A,A')`. `scripts/lrz/f4_compare.py` scores the arms by **coordinate**, not by count (§1y is precisely the lesson that counting species presence hides a wrong call), and prints the verdict. Binning fixed 2026-07-26; measurement attempted 2026-08-21 and inconclusive.** The LRZ A/B (jobs 5757307/5757308, identical inputs, one variable) returned *identical* wave structure and byte-identical results, because easy mode auto-picked targets at three distinct taxonomic ranks — a set legacy binning already grades correctly. **The right experiment is waves vs no waves on the uniform-depth GT genome set**, not rank-vs-legacy binning. Until it exists, attributing the melittin result to the wavefront is still not permitted. See F4.
10. ✅ **§1y — coordinate scoring built and frozen in CI (2026-08-30).** `tests/benchmark_truth/melittin_loci.tsv` is a new **28-model coordinate truth set** — one row per curated model, not per species — generated from the curated UGENE/manual GFFs in the Ant_Venoms tree, with every scaffold's assembly membership verified against the downloaded FASTAs. `scripts/benchmark/score_coordinates.py` scores a run by overlap and reports the **fraction** of each truth model recovered, because the fraction is what separates "recovered the gene" from "landed on its first exon". Wired into `run_all_sequential.sh` so a benchmark run now emits **both** numbers. Frozen by `tests/test_coordinate_benchmark.py` (18 cases), and CI now runs on `dev` as well as `main`.

    **What it shows on the two real runs** (re-verified 2026-08-30 against `local_runs/`):

    | run | `A0A6M3Z554_1` | overlap | frac of gene | region recall |
    |---|---|---|---|---|
    | `mel_det_rep2` | **hit**, start exact to the base | 132 bp | **15.4 %** | 2/12 |
    | `melittin_val_sw_20260630` | **missed**, nearest model 2 543 bp away | 0 | 0 % | **3/12** |

    Note the last column: the run that **lost** the gene scores *higher* on region recall. That is the §1y trap in one line — presence-style scoring can rise while the actual result goes to zero, which is exactly how the loss went unnoticed.

    Two things the truth set forced into the open: a curated model must not be scored against a run that used a **different assembly** (Cardiocondyla is curated on `GCA_019399895.2` while the truth table names `GCF_019399895.1` — 13 models, now reported `assembly_mismatch` and excluded rather than counted as fabricated misses), and a 110 kb *syntenic region* is not a gene model, so region rows are excluded from the fraction statistics.

### Benchmark completion (days, mostly waiting)

11. **TOGA on CoolMUC**, or a written exclusion + MCScanX/GENESPACE as the synteny comparator.
12. ✅ **Re-run done 2026-08-31** (job 5768083, 19 targets) — see **D.1b**. Presence F1 0.700 → 0.762 (recall 1.000), but the single gain is an artefact the coordinate scorer catches, and **precision did not move at all**. `scripts/lrz/benchmark_rerun.sbatch`. Must run on the **cluster**, not the laptop: the laptop holds 14 of the 19 target genomes, and the five Koludarov assemblies missing there include *Colletes*, *Euglossa*, *Xylocopa* (PRESENT) **and both deliberate false-positive tests**, *Melipona* and *Tetragonula* (LOST). Scoring without them changes the confusion matrix materially. All five are already staged on LRZ in `gt_melittin/targets`; the other 14 need uploading (5.6 GB).
13. Write `docs/BENCHMARK_RESULTS.md`. *(2026-10: the family benchmark is written up in
    `NEXT_SESSION_FAMILY_BENCHMARK.md` §0; a consolidated results page is still owed.)*

### Last

14. **Calibration pass** — one line of derivation for every parameter in `nextflow.config`, or delete it. Start with the distance family (§1t) and the classify thresholds.
15. **Website** from the existing markdown.
16. Talk to Ivan about the annotation agent (blocking, but not on the critical path for the tool).

---

# Part G — Robustness pass, 2026-10-05

A systematic check rather than a bug hunt: every flag a module passes against the
arguments its script defines, every parameter read against the ones declared, every
documented default against the config, static analysis for the read-before-assignment
class, edge-case inputs, and three end-to-end runs (the *Drosophila* Adh demo on the
committed and on the fixed code, and the yeast STE2 demo through the launcher). The
trigger was a crash report from a student that turned out to be a bug fixed on `dev` in
August and never merged to `main`.

## G.1 Defects found and fixed

| # | Defect | Effect before the fix | Evidence |
|---|---|---|---|
| G1 | `bin/plot_synteny.py` used a backslash inside an f-string expression | `SyntaxError` on Python 3.10 and 3.11, both allowed by `environment.yml`: PLOT_SYNTENY could not start. CI on `dev` was red. | CI run on `1e70fa6`; reproduced with a local 3.10. Now 782 tests pass on 3.10, 3.11 and 3.12. |
| G2 | `--auto_params true` applied nothing | Nextflow ≥ 25 ignores `params.put` without an error. The run logged every estimate as applied and every step kept its launch value. | Minimal Nextflow test (`put ok: n 5 -> 5`); yeast run: 5 estimated, 0 applied, search ran with `--max_intron 20000`. Now a warning with the flags to re-run with. |
| G3 | A preset overwrote values the user set | The resolver was promised "defaults < preset < user value" but was never told which values were the user's. | Code + new tests. Now any value that differs from the shipped default is kept (`user_set`, `user_set_kept` in `effective_params.json`). |
| G4 | `preset_paralog_discrimination` carried the TP53 family names | The preset is auto-applied to any large family, so an alcohol-dehydrogenase query ran with strict tokens `TP53, TP63, TP73, …`. | Adh run on committed code: 13 of 16 GOI models tagged `strict_family_downgrade`. No confidence changed on that dataset (the 13 were LOW anyway); on annotated targets a weak true call could have been demoted. |
| G5 | Nine boolean parameters tested by Groovy truthiness | Nextflow 26 hands every command-line value over as text, and `"false"` is true. `--require_paralog_panel false` (the documented opt-out) did nothing; `--disable_goi_hull_rescue false` switched the rescue **off**; `--enable_plm_search false` would request a GPU. | Minimal Nextflow test of parameter types; static test now forbids the pattern. |
| G6 | Region scores depended on the line order of the hit table | Hits were sorted by (chrom, start); ties kept input order, and the search writes its table in a thread-dependent order. | Old-vs-new Adh runs: identical hits as a set, different order, consistency 0.896 vs 0.917. Five orders of one file gave five scores; now one. Gene calls were never affected. |
| G7 | Rescue-pass gene models were never published | A hull-rescued ortholog (27 of the 51 confident calls of the SP benchmark run) existed in the output only as a span in the report. | Output trees of all runs on disk. Now `rescue/locus_<N>/`. |
| G8 | Three switches existed in scripts only | `--disable_ambiguous_tier` (documented as the kill switch), `--disable_distant_synteny_rescue`, `--legacy_strand_score` had no effect on a pipeline run. | Wiring audit. Now parameters; all plot options reachable through `--plot_extra_args`. |
| G9 | The query was passed through verbatim | Lowercase, alignment gaps, digits and spaces of numbered lines, and a trailing `*` reached the search inside the "normalised" query. | 15 awkward inputs. `normalize_query.py` now cleans or rejects; it had no tests. |
| G10 | 78 % of a cluster console log was one warning | `env { LLM_API_KEY = … ?: '' }` exported an empty variable, and Nextflow warned for every task (12,929 of 16,528 lines). | Console logs of the `fix7` runs. |
| G11 | The end-of-run summary could contradict itself | "GOI found in 4 genome(s)" (genomes with *any* model) next to "absent in 1". A default easy-mode run also warned `max_genomes=0 … <3 target genomes` although 0 means automatic. | Adh run. The summary now prints the report headline. |
| G12 | Assembly quality columns were misread in two more NCBI queries (found in G.4) | `xtract` drops absent fields unless `-def NA` is given, and an assembly record has no top-level scaffold or contig count. The scaffold and contig N50 were therefore read as scaffold and contig *counts* in `fetch_related_genomes.get_related_species` (every easy-mode target search) and `fetch_home_genome.find_any_genome`. A scaffold-level assembly with an N50 above 500 kb failed the quality gate and was dropped, and within one species the assembly with the smallest contig N50 ranked first. The same defect was fixed in May for one call site only (§1l). | The easy-mode run logged *S. paradoxus* as `scaf=903028, contigs=835812, N50=NA/NA`; confirmed on the live NCBI record. `tests/test_assembly_docsum_columns.py` (9). |

Smaller: a stale local copy of `species_from_leaf` shadowed the shared one in the matrix
plot; `rescue_goi_hull.py --help` crashed on an unescaped `%`; six `[DEBUG …]` prints
made up 31 % of a search log; NCBI query pipelines had no timeout; genome ordering and
report file lists depended on directory-listing order; the validation message advertised
a `--help` that does not exist; the manifest said version 1.0.0; the long-failing
`test_rbh_batch` (a short-peptide fallback that only triggered on an error current
MMseqs2 no longer raises).

New guards: file names with whitespace or special characters, and targets that differ
only by extension, stop the run at launch; ignored task failures are counted in the end
summary.

## G.2 What the end-to-end runs showed

| Run | Result |
|---|---|
| Adh, 4 *Drosophila*/*Anopheles* targets, committed code | 3 HIGH (98.8 / 97.3 / 90.2 %), 13 LOW; exit 0, 62 tasks |
| Same, fixed code | identical calls and identical gene models byte for byte (family-gate attributes aside); one region score differs through G6 |
| Yeast STE2 through `./run_synvoy.sh`, fixed code, `--auto_params true` | 1 HIGH (66.0 %) + 2 MEDIUM (47.8 / 41.9 %), all three matching the home STE2 gene; 55 tasks, 3 m 56 s. In July the same demo read "1 high + 32 medium". |

Home-locus selection was tested for order dependence too (original, reversed and shuffled
hit tables, a single-locus and a three-locus query): identical output.

## G.3 Found, not changed

- **`RECIPROCAL_BEST_PARALOG` (§1j Phase B) cannot fire.** It compares calls against
  several home query sequences, but since the F8 fix the query that reaches it is always
  the one normalised protein. It still spawns a task per locus and genome. Locus
  ownership does the job it was meant to do. Remove it or feed it the paralog panel.
- **A multi-sequence `--query` uses only its first record** (warned in the task log
  only). The paper's GR1 query file holds ten sequences.
- **The LLM estimate never reaches the six preset-covered parameters**, and (G2) none of
  the others either. Making it effective means routing it through the `settings`
  channel, as presets are.
- **The plots and `homology.tsv` do not show rescue models or adjudication verdicts**
  (§1p). In the yeast demo the report's HIGH call is a hull-rescue model the figure does
  not draw as such.
- **Five declared parameters are inert**: `keep_intermediate`, `max_retries`,
  `multi_profile`, `multi_profile_max_jobs`, `locus_ownership_synteny_window`. The docs
  now say so.
- **Four packages in `environment.yml` are not imported by any pipeline code**
  (`plotly`, `taxopy`, `psutil`; `biopython` only by one test).
- **Each new work directory rebuilds the conda environment** (3.4 GB, about two minutes),
  because no `conda.cacheDir` is set.
- **`main` is behind `dev`** by everything since 2026-07-21, including two crash fixes and
  the easy-mode rescue bug. Anyone who cloned `main` is running that.
- The melittin ground-truth fixture has not been regenerated since 2026-03-28;
  `scripts/reproduce_annotation.py` is dead.

## G.4 Pre-merge check on the students' use

Two more runs on the pushed code (`baccf29`), both through `./run_synvoy.sh` with its
default profile, on a 12-core, 15 GB laptop.

| Run | Result |
|---|---|
| Firefly luciferase, pro mode — *Photinus pyralis* home, query `XP_031329057.1`, targets *D. melanogaster* and *A. gambiae*. This is one of the two cases a student reported as crashing on `main`. | exit 0, 130 tasks, none failed, 39 min. Four home loci; `preset_paralog_discrimination` auto-applied, now without family tokens. No HIGH call, 3 MEDIUM (40–42 % identity, each reached from more than one home locus), 1 AMBIGUOUS, 81 LOW. |
| Easy mode — yeast STE2 by accession (`--query_id D6VTK4 --max_genomes 3`). The first easy-mode run on this code. | exit 0, 44 tasks, none failed, 2 HIGH. It took 35 min, 30 of them in three NCBI queries that ran into the 600 s limit, and it fetched 2 of the 3 requested targets. |

**The crash itself.** `tests/test_cluster_grs_strong_synteny.py`, run against `main`'s
`cluster_grs.py`, reproduces the student's traceback line for line (`line 1059 …
UnboundLocalError … 'strong_synteny'`); on `dev` it passes. The luciferase run did not
pass through the crashing state: `main`'s script also completes on that run's eight
clustering tasks. The crash needs the *first* scored cluster to overlap a GOI hit, and
that depends on what the search wrote, which differs between the two branches.

**What the luciferase calls mean.** Flies are not bioluminescent. The three MEDIUM calls
are the closest members of the same enzyme family, found from several home loci at once.
A MEDIUM call on a distant target says "best family member in a neighbourhood with some
synteny", not "luciferase".

**Easy mode, found and fixed.** G12 above, and one more:

| # | Defect | Effect before the fix | Evidence |
|---|---|---|---|
| G13 | An NCBI query could block on its own error output | `run_piped_command` reads only the last stage of `esearch \| efetch \| xtract`; the earlier stages wrote stderr into a pipe nobody drained, so a stage stopped for good once it had written about 64 kB, and whatever it had said was lost. entrez-direct prints the whole request for every call it has to retry. Each stage now writes stderr to a temporary file, and its tail is printed when the query fails or times out. | 23.5 kB of stderr for one 1,235-assembly query on a connection where each request failed once. `tests/test_ncbi_pipe_chain.py` (4): a stage that floods stderr hangs the pushed code and passes now. |

What G12 changes in practice, for the firefly case in easy mode (*Photinus pyralis* home,
six targets): the old code took the four Lampyridae references with the *smallest* contig
N50 (0.5, 2.0, 3.1 and 7.7 Mb) and left out *Aquatica leii* (10.8 Mb); for the two
Coleoptera slots it took assemblies with 0.14 and 0.16 Mb. The fixed code takes
*A. leii* first and fills the Coleoptera slots with 125 and 102 Mb assemblies.

**Easy mode, found and not fixed: the stall.** G13 is not what cost the 30 minutes. The
same 3,407-assembly query was run through the pushed and through the fixed helper with a
300 s limit: both ran into it, with nothing on stderr. The query is simply too large. Each
taxonomy level asks NCBI for the record of *every* assembly under the taxon, and how long
that takes varies widely: 2,357 records came back in under 90 s inside the pipeline run,
and not within 8 minutes an hour later. A level that hits the limit is skipped, the run
continues with fewer genomes, and the console says only "downloaded 2 target genome(s)".

| Taxon | Assemblies at NCBI (2026-10-05) |
|---|---:|
| *Apis* / Apidae / Hymenoptera | 39 / 228 / 1,753 |
| Lampyridae / Coleoptera | 12 / 1,235 |
| Hydrozoa / Cnidaria | 72 / 605 |
| *Saccharomyces* / Saccharomycetaceae | 2,357 / 3,407 |
| *Homo* / Primates / Mammalia | 3,211 / 3,792 / 6,969 |
| Insecta | 12,093 |

So a human or yeast gene in easy mode can stall at the first level, and any insect run
that has to climb to class level will. The fix is to narrow the search on the server
(reference and representative assemblies first, the rest only if the budget is not
filled). That changes which genomes are chosen and needs its own test, so it was not done
here. An `NCBI_API_KEY` raises the rate limit. Pro mode does not use these queries.

Suite after G12 and G13: 880 passed, 0 failed.

## G.5 Figure review of the latest runs (2026-10-06)

All figures of the four latest cluster runs (APYR, DPP4, SP, melittin; 17–18 Sep) and of
the four local runs of G.2 and G.4 were laid out on one page and looked at. Before any
change, replaying each cluster run's own plot command with the current plotting code
gave 78 of 78 byte-identical files: the September figures were what `dev` drew.

| # | Defect | Effect before the fix | Evidence |
|---|---|---|---|
| G14 | The grid's `×N` badge counted hits anywhere in the genome | Above a flanking arrow, `×N` was the number of models named after that home gene in the whole track, on any scaffold. A paralog-rich gene looked like a local tandem array. Now it counts the models on the drawn gene's scaffold that lie no farther from it than the neighbourhood is long; the legend says so. The GOI column is unchanged. | Yeast, easy mode: `×9` above a 100 % ortholog, the other eight being 26–48 % hits on other chromosomes; 23 such badges in that grid, none after the fix. SP, three loci: 46 badges → 23. DPP4: 14 → 11. APYR: 3 → 3 (real local copies). |
| G15 | The home GOI label of the synteny plot ran through the subtitle | The label climbs at 45° into a header of fixed height. A long product name crossed the subtitle and was cut off at the top edge. The header now grows as far as the first track's labels need, a name longer than 34 characters is cut with an ellipsis, and the small "GOI" tag is dropped when the label runs through it. | Firefly luciferase, "4-coumarate--CoA ligase 1-like ×5". Short labels in a busy track (melittin) keep the old header. |
| G16 | Home row labels | The home row printed its scaffold twice ("NC_037641.1 • NC_037641.1"), and in easy mode its species label carried UniProt's strain qualifier ("Saccharomyces cerevisiae (strain ATCC 204508 / S288c)"). | Every synteny plot; the yeast easy-mode figures. |
| G17 | Easy mode could offer the home species as a related genome | Found through G16. The related-genome search compared the home species for equality with NCBI's species name. With a strain in the home name (UniProt's parenthetical form, or NCBI's "Saccharomyces cerevisiae S288C") it never matched, and the taxonomy lookup of the long name failed over to the genus. | On the 2,357 real assembly records of the genus, the home assembly `GCF_000146045.2` itself was the fourth candidate. `tests/test_fetch_related_home_species.py` (10). |

Tests for G14–G16: `tests/test_plot_grid_copies_and_labels.py` (14). All figures were
re-drawn afterwards; the README figure is the re-drawn melittin grid.

**Exon counts against NCBI.** NCBI lists 14 exons for the apyrase home gene; the figure
shows 11. Nothing coding is missing. The gene has three annotated transcripts with 11
coding segments each, and NCBI's number also counts three non-coding 5′ UTR exons of
the two predicted variants. The home gene model is taken from the home annotation, not
re-predicted:
all eight runs report `method: gff_annotation` in `goi_info.json`, with the model
protein equal to the query (APYR: 576 aa, start and stop present). The alignment-based
model is the fallback for a home genome without a matching annotated gene.

Suite after G14–G17: 904 passed, 0 failed.

## G.6 Algorithm description, parameter docs, block merging (2026-10-06/07)

[ALGORITHM.md](ALGORITHM.md) now describes what the pipeline computes, step by step,
written from the code. Reading the code for it found two more defects.

| # | Defect | Effect before the fix | Evidence |
|---|---|---|---|
| G18 | Eight parameters were documented as doing things the code does not do | `max_consecutive_empty_blocks`, `gap_min_size`, `gap_evalue` and `gap_min_alnlen` have no effect; `gap_min_identity` and `gap_max_hits` only bound last-resort `raw_hit` rows; `gff_search_window` is the home-gene lookup distance; `gap_search_window` is the margin for GOI hits and for miniprot. There is no gap-filling search. | Code reading (`process_single_genome`, `annotate_exons_from_hit_list`). PARAMETERS.md, USAGE.md and the config comments corrected; no behaviour changed. |
| G19 | A merged block added its parts' gene counts and kept only the first part's gene list | Two copies of one flanking gene less than 300 kb apart read as a 2-gene block and passed `min_block_genes = 2`, so the window was searched although the neighbourhood had a single anchor. The genes of the later parts were invisible to the "seeded by the GOI's own gene" exception and to the 80-block cap. Now the merged block holds the union of the genes. | `tests/test_merge_block_gene_count.py` (6; 4 fail on the old code). Measurements below. |

**G19, measured.** Replaying block building on the saved genome-wide hits of two local
runs (firefly luciferase against two flies, yeast STE2) reproduces the block counts in
the run logs: 89 of the 544 blocks that passed the filter held one distinct flanking
gene. An A/B on the cluster with two code snapshots that differ only by this fix
(control `pre8`, fix `fix8`; jobs 5831701–08), plus the luciferase case locally:

| Run | Blocks searched, control → fix | HIGH / MEDIUM / AMBIGUOUS records that differ | Gene-model rows that differ |
|---|---|---|---|
| APYR, 33 genomes | 224 → 209 | 0 of 32 | 0 of 331 |
| DPP4, 33 genomes | 352 → 351 | 0 of 32 | 0 of 149 |
| SP, 33 genomes, 3 home loci | 5,675 → 5,028 | 0 of 122 | 19 LOW gone, 2 LOW new, of 2,129 |
| Melittin, 19 genomes | 46 → 46 | 0 of 10 | 0 of 8 |
| Luciferase, 2 genomes, 4 home loci (local) | 525 → 436 | 0 of 6 | 2 LOW gone, 1 MEDIUM new, 1 LOW with another identity, of 89 |

Recall, precision and the set of final calls are identical in both arms for all three
families. The one new MEDIUM row (luciferase) is the second half of the defect: a block
that was seeded by a rescued GOI hit lost that gene from its list when it was merged,
was not protected, and fell outside the 80-block cap. With the fix it is searched; the
model it yields was already in the report from another home locus, and as a new member
of the query set it changes the identity of one LOW row in the next genome. Run times of the two
cluster arms are not comparable: they ran on different nodes, and the melittin run, which
did identical work in both arms, differed by 35 %.

**The control arm is also the first cluster run of the code since G.1.** Against the
last validated run (17 Sep): APYR, DPP4, SP and melittin have the same final calls,
the same recall and precision, and identical gene-model rows in the search output. One
label differs: 14 APYR records carry the class `confident_goi` where they carried
`synteny_hull_rescue`. Both models exist in both runs, at the same coordinates and with
the same identity; the report keeps the first of two tied models, and since the October
robustness pass (`565bf81`) the result files are read in sorted order, main-path model
first, and no longer in directory order.

**Found, not changed.** The hull rescue looks for an existing HIGH GOI model inside the
span of the HIGH flanking models but searches that span plus 100 kb on each side. A gene
the main path found just outside the span is modelled again: in the APYR control run the
rescue produced a model in 28 of 33 genomes, all of which the main path already had as
HIGH. The duplicates are merged in the report, so the counts are right; the work is
wasted, and before the sorted read the reported class of such a gene was arbitrary.
Changing the check can change which model is reported where the two differ, so it needs
its own benchmark run.

The rescue also labels its models with the shipped identity thresholds (HIGH from 50 %),
not with the thresholds of an auto-applied preset. In the DPP4 run the preset raised the
main path's HIGH threshold to 55 %; 16 genes at 52.1–54.8 % identity are MEDIUM there,
and the same 16 are HIGH in the report through the rescue. All 16 are curated orthologs,
so recall and precision are unaffected; the confidence label of such a gene depends on
which path reports it.

Suite after G18–G19: 910 passed, 0 failed.

---

## What to say about SynVoy right now

> A synteny-guided locator that finds the *neighbourhood* of a divergent gene when sequence
> search cannot — validated on melittin, yeast STE2 (3/3), Drosophila Defensin (3/3), and
> me31B across ~250 My. It recovers the right locus far more often than BLAST/MMseqs2
> baselines (species-presence recall 1.000 vs 0.250). Its weakness is the complement of its
> strength: it asserts presence almost everywhere (16 of 19 benchmark species, ABSENT in
> none, precision 0.615), and locating a neighbourhood is not the same as modelling the gene
> in it — by coordinate overlap only 1 of 9 curated models is recovered over half its span.

Add, since 2026-09: *on ten curated venom-gene families across 33 ant genomes it recovered
30/31, 28/28, 32/32 and 31/32 of the reachable genes in the four near-single-copy families,
with the reachable set fixed in advance; tandem arrays and loci absent from the home genome
are out of reach by construction.* And replace the coordinate sentence above with the
September measurement: internal splice sites of recovered models match the curated ones,
the errors are at the gene ends, where 97 of 182 now land exactly on the curated coordinate.

That is defensible under review. **Three things must NOT be said** (2026-08-31):

1. **Do not attribute the melittin result to the expanding wavefront.** Measured with a proper
   control and a replicate: three waves vs one wave, byte-identical output. Neighbourhood
   restriction is what produces that result (F4). Equally, do not claim the wavefront is
   useless — that set had only two distinct distances and no target needing a stepping stone,
   so the mechanism was never stressed. Both overstatements are wrong.
2. **Do not quote a species-presence recall without the coordinate number beside it.** On the
   current benchmark the only recall gain over May is a gene model **63,936 bp** from the gene
   it is credited with finding (D.1b).
3. **Do not say the over-calling is fixed.** The 2026-07/08 adjudication pass fixed *labels*;
   benchmark precision is unchanged since May (5 FP, same species).
