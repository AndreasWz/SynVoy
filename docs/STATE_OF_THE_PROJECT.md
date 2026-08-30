# SynVoy — State of the Project

**Compiled 2026-07-26 by direct code and output inspection.** Every number below was
re-derived from the repository or from run outputs on disk; nothing is quoted from an
earlier summary without re-checking. Where a claim could not be verified it says so.

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

---

## 0. One-paragraph status

SynVoy is a working synteny-guided search engine wrapped in a scoring, labelling and
reporting layer that is not trustworthy. The discovery layer is validated across four
clades and is robust to large parameter changes. The adjudication layer — region ranking,
confidence labels, the p-value, the report headline — contained at least four defects
serious enough that its numbers should not be published as-is; as of 2026-07-26 those four
are fixed, and the 2026-08-21 cluster run confirmed three of them on real output.

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
| `tandem_copy` | — | `id ≥ 40` ∧ `qcov ≥ 0.35` |
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

## F4 ⭐⭐⭐ The wavefront is inoperative in 37 of 44 runs — ⚠️ CODE FIXED, measurement still owed

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

**⬜ Still owed — but the mechanism now exists.** The comparison that answers the question
is *waves vs no waves*, not rank-vs-legacy binning. `--disable_wavefront` (added 2026-08-22)
puts every genome in one parallel wave, so the query DB is never augmented between genomes;
run it against a normal run on the **uniform-depth GT genome set**, where the gradient
actually differs. It is a measurement flag, not a performance option.
**Until that run exists, do not attribute the melittin result to the wavefront.**

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

## F7 Documentation-layer defects (fixed in the 2026-07-26 docs pass)

`--keep_intermediate` is declared and read by nothing; `--no_anchor_grid` is a
`plot_synteny.py` switch that the Nextflow module never forwards; `--help` does not exist
though `main.nf:714` still advertises it; 59 parameters were documented nowhere. All now
corrected in the docs — but `main.nf:714` and the two dead parameters are still live in
the code.

---

# Part C — Known open items (from `docs/TODO.md`, re-verified)

| § | Item | Severity | State |
|---|---|---|---|
| **§1x** | **Silent discard** — a confident hit on the correct gene is rejected by the synteny gate and never mentioned. `demo_Def`: *Anopheles* defensin found at 54–72 % inside the real RefSeq gene, reported as "no ortholog". | ⭐⭐ | ✅ **Fixed 2026-07-26** — report block, per-rejection logging, and the plot no longer claims absence. Off-block marker ⬜ (needs a `main.nf` channel) |
| **§1y** | **Coordinate regression** — the ant melittin was hit in 8 runs; the newest run is 2 543 bp off and nobody noticed, because scoring is species-presence not coordinate overlap | ⭐⭐ | Open; best-ever hit was 132 bp of an 855 bp gene |
| §1v | Classifier ignores the synteny evidence the run computed | ⭐⭐ | ✅ **Fixed 2026-07-26** — identity floor; per-call support measured and closed won't-do (see F6) |
| §1u | `goi_block_flanking` drives ranking but is not emitted | ⭐ | ✅ **Fixed 2026-07-26** — emitted in the scores TSV alongside `synteny_score` |
| §1z | `QueryCoverage` read as target recovery; target coverage never reported | 🐛 | Open |
| §1p | Post-processing verdicts never reach the plot or the tree | ⭐ | Open |
| §1p.1 | Target exon composition not passed to plot → synthetic uniform exons drawn | ⭐ | Drawing side in progress (uncommitted); data side untouched |
| §1t | Every distance parameter is metazoan-scale (found by pointing it at an 11 Mb yeast) | ⭐ | Open |
| §1w | Duplicate regions with byte-identical spans | 🐛 | Open |
| §1q, §1s | `strong_synteny` read before assignment; numeric CLI params crash the run | ✅ fixed | **Fixed but uncommitted** — see Part F |
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

9. ⚠️ **F4 — binning fixed 2026-07-26; measurement attempted 2026-08-21 and inconclusive.** The LRZ A/B (jobs 5757307/5757308, identical inputs, one variable) returned *identical* wave structure and byte-identical results, because easy mode auto-picked targets at three distinct taxonomic ranks — a set legacy binning already grades correctly. **The right experiment is waves vs no waves on the uniform-depth GT genome set**, not rank-vs-legacy binning. Until it exists, attributing the melittin result to the wavefront is still not permitted. See F4.
10. **§1y — re-score the melittin benchmark by coordinate overlap** and freeze it in CI. Publish the coordinate number next to the species-presence number.

### Benchmark completion (days, mostly waiting)

11. **TOGA on CoolMUC**, or a written exclusion + MCScanX/GENESPACE as the synteny comparator.
12. **Re-run the whole benchmark on current code** — the numbers in D.1 predate three months of fixes.
13. Write `docs/BENCHMARK_RESULTS.md`.

### Last

14. **Calibration pass** — one line of derivation for every parameter in `nextflow.config`, or delete it. Start with the distance family (§1t) and the classify thresholds.
15. **Website** from the existing markdown.
16. Talk to Ivan about the annotation agent (blocking, but not on the critical path for the tool).

---

## What to say about SynVoy right now

> A synteny-guided locator that recovers divergent orthologs sequence search misses —
> validated on melittin (11/12 vs 2/12 for BLAST/MMseqs2 baselines), yeast STE2 (3/3),
> Drosophila Defensin (3/3), and me31B across ~250 My — with a confidence-labelling layer
> that currently over-calls and is being rebuilt to consume the synteny evidence the
> search already computes.

That is true, it is defensible under review, and it does not depend on any of the numbers
this document just invalidated.
