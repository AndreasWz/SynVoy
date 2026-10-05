# Proposal — iterative family-recovery analysis on Ivan's venom families

**For the next session.** Written 2026-09-10, self-contained: everything needed to start is
below, including where the data lives and what was already checked.

**The idea (Frank's):** take one sequence from one species, ask SynVoy to find the rest of
that family across the other species, and see whether it can. Then work out *why* it fails
where it fails — parameters? a real bug? or is the task genuinely too hard, so that "fixing"
it would just be overfitting? Change something, retry. The point of using **many** families
is that you get independent retests instead of one anecdote.

---

## 0. Status after the 2026-09-11 session — read this first

**Where it stands.** Phase 1 has run (all five easy families, default flags, *Apis* home, 33
targets), plus the two tune families. The determinism fix is confirmed on LRZ: two runs under
different hash seeds gave 20/20 identical GFFs and identical adjudicated records (job
5780371). Numbers below are family-level only; everything naming a locus or gene is in the
gitignored `tests/benchmark_truth/family_benchmark_private.md`.

**1. Where the seed sits at its locus's curated address, SynVoy finds essentially everything.**
Locus-level recall on the reachable genes, baseline code (b3fcfe6 + the determinism fix):

| family | reachable | recovered | the misses | locus precision | the off-truth calls |
|---|---:|---:|---|---:|---|
| APYR | 31 | 30 | 1 curated model not homologous to the seed (E ≈ 0.7) | 0.97 | 1: the real outgroup ortholog, at the address |
| DPP4 | 28 | 28 | — | 0.90 | 3: DPP4 at the address in 2 species the curation lists as lacking it, + a second copy |
| HYAL | 32 | 32 | — (4 of them are shared-address, `unresolvable_by_synteny`) | 0.91 | 3: SynVoy defects, fixed below |
| KAZA | 32 | 31 | 1 gene on a 3.5-kb contig (no room for a flanking gene) | 0.97 | 1 chance tandem hit, fixed below |

Every miss is layer A, and every one is explained without a search bug. Every off-truth call
was checked against the genome. Four are real orthologs at the address the curation lacks
(callcheck: the call's window is best explained by the seed's own home gene; its flanking
models sit in home order), and the rest trace to the four classifier defects below.

**2. Four of the ten curated *Apis* seeds do not sit at their own locus's address**
(SPIN, SECP, PLA2 — all three seeds — and the VA CAP seed). The outgroup labels were
assigned by sequence, not by position. The runs from those seeds searched a neighbourhood the
ant members do not share: SPIN 0/32 and SECP 0/30 locus recall, all layer A. That measures
the seed, not SynVoy. **Amendment 2** (below) audits every seed's address by hand before any
further run and adds *address seeds*: the *Apis* gene that does sit at the address, found by
reciprocal best hit and confirmed by its neighbours. For SECP it is an unannotated ORF inside
a RefSeq lncRNA; the curated SECP seed turns out to sit at the family's *other* locus address,
so it is scored under that label instead. `seedcheck`'s RBH flag alone is not a verdict: for
APH the RBH alternative is a sequence-closer paralog elsewhere and the curated seed is right.

**3. SynVoy defects found and fixed** (all layer C, all in `bin/iterative_search_runner.py`,
CLAUDE.md known-bug §19). Across the seven runs, **every final call with query coverage < 0.5
was off-truth (78) and all 157 on-truth calls had ≥ 0.5** — the misses of precision were all
cheap alignments dressed up as calls:
* tandem copies got MEDIUM on identity + coverage alone, no synteny context: a short query's
  best chance alignment in a window clears 40 % / 0.35. 77 MEDIUM calls in SECP + SPIN, all
  off-truth, 16–33 bits. Now they need the strong-flanking bar or qcov ≥ `medium_min_qcov`;
* the fallback chain could join two alignments of the *same* genome window as two exons;
* the fallback's coverage was the query *span*, not the aligned residues (0.90 vs 0.06);
* one flanking gene lifted a model the pipeline itself labels a fragment to MEDIUM.
Re-graded offline on the finished runs, the changes demote 71 of the 77 chance tandem calls
(the other 6 sit in strongly supported blocks) and all 4 other targeted calls, and touch no
on-truth call; on the melittin benchmark they demote one weak-block tandem call (confirmed
on LRZ: identical coordinate score, one off-truth tandem MEDIUM gone). The Smith-Waterman
placeholder statistics found here were fixed after all, in the audit (item 6).

**4. Measured verdicts (design §5, deliverable 3), baseline code, all adjudicated records:**
the AMBIGUOUS tier holds **0 of 137** records on a curated gene, so keeping it out of the
counts is right; the coverage demotion (`identity_coverage_decoupled`) was right 16/16; the F9
phylo verdict carries no signal here (3 of 6 `phylo_discordant` records are real orthologs,
66 chance hits are `phylo_concordant`). Those F9 numbers were measured on frame-broken
proteins (item 6); re-measured on the fixed code they do not improve: across HYAL, APYR and
PLA2 every `phylo_discordant` record (6) is a real ortholog and no off-truth record is flagged.

**5. An honest limit, not a tuning target.** SECP's second locus is nested inside a very large
host gene: the search blocks reach 23 of its 26 genes and SynVoy models 13 of them, but only
at 25–46 % identity in blocks with 2–3 flanking genes — the same statistics as the chance
hits. They stay LOW, correctly.

**Amendments** (each hashed before the runs it governs; files gitignored):
* *Amendment 1* (`family_preregistration_amendment1.tsv`, SHA-256
  `e5427d318bd1d2502e3d0e4ef845a20c71a94ebf225fde3fba48fa418f9bb27b`): an
  `unassigned_label` tier (survey/scaffold labels: family-level only) and a `seed_unrelated`
  flag.
* *Amendment 2* (`family_preregistration_amendment2.tsv`, SHA-256
  `27600b7396c1ae6f0530133662aa5f9f82925a29678e5241986667bca2077c38`): the seed address
  audit and the address-seed arms; the `seed_unrelated` flag becomes a length-aware E-value
  (> 1e-3; the raw-score floor flagged every member of a 77-aa family); an `assembly_limited`
  flag (contig < 20 kb). Both flags are reported both ways, never re-tiered.

**6. Algorithm audit (2026-09-11), `docs/ALGORITHM_AUDIT_2026-09.md`.** Tracing the search
core step by step on these runs found two defects that reach the output, both fixed (CLAUDE.md
known bug §20): every multi-exon GOI model's *protein* was translated partly out of frame (the
GFF was right; the sequences feeding seeds, tree, paralog RBH and F9 were not), and every
Smith-Waterman hit carried a fake E-value with its raw score as the bit score (six chance hits
per window, all passing every filter). Smith-Waterman now uses BLAST 11/1 with Karlin-Altschul
E-values (decision 2026-09-11), fallback identity is length-weighted, and a window's GOI loci
are ranked by aligned coverage. A synteny block E-value was measured as a replacement for the
raw flanking count (decision: measure first): AUC 0.91 vs 0.88 for on- vs off-truth records,
because the off-truth MEDIUM calls it would catch (40 of 46 in chance-level blocks) are the
ones the item-3 fixes already demote — not gated.

**7. The fixed code (LRZ snapshot `fix4` = items 3 + 6) on the tuning arms and the ladder**
(locus recall on reachable genes / final calls / off-truth final calls):

| arm | baseline | fix4 | note |
|---|---|---|---|
| HYAL | 32/32 · 36 · 3 | 32/32 · 33 · 0 | off-truth AMBIGUOUS 13 → 5 |
| APYR | 30/31 · 32 · 1 | 30/31 · 32 · 1 | the off-truth call is the real outgroup ortholog |
| KAZA | 31/32 · 63 · 1 | 31/32 · 33 · 0 | |
| DPP4 | 28/28 · 32 · 3 | 28/28 · 32 · 3 | the 3 are the real orthologs of item 1 (callcheck: seed gene); the `fix1` and first `fix4` runs died on the LRZ home quota, rerun alone |
| SPIN, curated seed | 0/32 · 27 · 27 | 0/32 · 1 · 1 | seed at the wrong address (item 2) |
| SPIN, address seed | — | 9/32 · 21 · 11 | 16 of 23 misses are layer C (AMBIGUOUS: 11 on-, 10 off-truth) |
| SECP, second locus | — | 0/26 · 0 · 0 | item 5; misses now layer B (22): without Smith-Waterman noise the chance-level models are no longer built; `fix1` still made 5 off-truth calls |
| SECP, address seed | — | 0/30 · 0 · 0 | all 30 misses layer C; 21 of the 30 members are not significantly similar to the seed (E > 1e-3) |
| PLA2, address seed | — | 14/15 · 32 · 18 | all 18 "off-truth" calls are reciprocal best hits of the seed gene, one per species at the address: a 1:1 ortholog group the curated locus lacks |

Melittin regression check (`fix4` vs the post-determinism reference, coordinate truth):
adjudicated records 95 → 10; MEDIUM-or-better calls 7 → 4: the three removed were chance
tandem calls in species without a truth gene (one of them a gene-lost species), and one weak
true family member was gained; one true bee copy drops from MEDIUM to AMBIGUOUS (its MEDIUM had come from a chance
Smith-Waterman hit posing as a tandem partner; on its own evidence it is indistinguishable from
the lost-gene records the AMBIGUOUS tier exists for); the Smith-Waterman pointer to the
divergent ant member (E 18) is gone, as decided. False positives at the gene-lost control
loci 2 → 0.

**Infrastructure note.** The LRZ `$HOME` holding the Nextflow work dirs has 27 GB of headroom;
one stage-1 genome search peaks at ~2 GB of MMseqs2 temp and a family runs up to four at once,
so 9-13 concurrent families exceeded the quota late in the run (MMseqs2 "Can not write to data
file", search aborted, no report): `fix1` DPP4 and curated SPIN, first `fix4` DPP4. Run at most
one family per job and one job at a time from `$HOME`, or move `work/` to the fast share.

**Item 8 — the hold-outs (VA, SP, APH), baseline vs the frozen `fix4`.** Five arms, each run
once with each code, as a one-at-a-time Slurm chain (jobs 5785138-47) submitted before any
hold-out number was seen; no tuning on these families, and the results below are reported as
they came. Recall is over the genes the arm can reach (the pre-registered ceiling); "off" =
final calls that overlap no truth gene of the family.

| Arm | Baseline recall · calls · off | `fix4` recall · calls · off | what changed |
|---|---|---|---|
| VA arm 1 (address seed) | 28/55 · 39 · 8 | 29/55 · 32 · 2 | the 5 chance MEDIUM calls (36-132 bp windows, 40-44 % identity, no home protein explains them) and one 73 %-identity call on an unrelated zinc-finger gene are gone; one more true gene recovered (a layer-C miss became a MEDIUM call); precision 0.79 → 0.94, and both remaining off-truth calls are best explained by the seed gene itself (candidate uncurated orthologs) |
| VA arm 2 (curated seed, address MATCH) | 26/30 · 27 · 1 | 26/30 · 27 · 1 | identical calls; one unrecovered gene moved from layer B to layer A (in `fix4` no search block covers it) |
| SP (curated seed, address arm) | 29/30 · 48 · 17 | 29/30 · 51 · 19 | same genes, better models: median on-truth call identity 53 → 69 %, `confident_goi` 2 → 21, hull rescues 28 → 15; 3 new off-truth MEDIUM calls at 40-52 % identity, and 16 of 17 (18 of 19) off-truth calls are best explained by the seed gene itself |
| APH (curated seed) | 28/31 locus · 41 · 11 | 28/31 locus · 41 · 11 | identical genes and calls; median on-truth identity 59 → 72 %, MEDIUM 8 → 4, `confident_goi` 24 → 29. The 48 address-tier genes stay unrecovered in both, as pre-registered (the three APH tandem loci share one address) |
| VA arm 3 (curated seed, address MISMATCH) | 1/55 · 35 · 17 | 0/55 · 29 · 15 | **the pre-registered failure, and it failed**: 39 of the misses are layer A (the seed's neighbourhood is not the family's), 13 truth genes are found without an anchor, and the calls land on other families' loci (16 → 13) or on orthologs of the seed gene itself (16 → 14 off-truth, `query_gene`). `fix4` drops the baseline's single accidental on-locus call |

All ten runs: exit 0, 0 disk errors, code hash verified per job (`baseline` =
`b3fcfe6+pipeline:e6556fa5a7a5`, `fix4` = `b3fcfe6+pipeline:3a610ea4352d`); chain finished
2026-09-13. Per-arm detail, seeds and coordinates are in the gitignored
`tests/benchmark_truth/family_benchmark_private.md`.

**What the hold-outs say.** Recall never dropped on an arm the baseline could reach (+1 gene on
VA arm 1, equal on the other three reachable arms), precision rose where the audit's fixes
targeted chance calls (VA arm 1: 8 off-truth → 2, 0.79 → 0.94), and the pre-registered failure
failed. The clearest effect is one the score does not show: with the frame fix the model
proteins that feed the next wave are correct, so the same loci are now modelled by the main
in-block path instead of the hull rescue, at 20 points higher identity — e.g. one SP locus goes
from `Identity=54.9 MEDIUM probable_goi` to `Identity=74.1 HIGH confident_goi` on the same
miniprot path. Median on-truth call identity: SP 53 → 69 %, APH 59 → 72 %, VA arm 1 unchanged
(58 → 59 %, that family's models are mostly single-exon). The one cost is 3 new off-truth MEDIUM
calls on SP at 40-52 % identity. No parameter was tuned on any of these families.

**Item 9 — model boundaries (2026-09-15).** Comparing every recovered gene's model with the
curated annotation exon by exon (216 genes, 1,439 curated exons): 62.6 % of exons exact on both
boundaries, 24.4 % on one, 12.2 % missing, and 147 model exons with no curated counterpart.
Internal splice junctions are right; the errors sit at the termini. Mechanisms, fixes and the 12
cases read by hand are in `docs/ALGORITHM_AUDIT_2026-09.md` §B4. Two fixes shipped — terminal
start/stop recovery (`refine_model_termini`, with `StartCodon`/`StopCodon` now on the model) and
fused-model supersession in the report dedup — plus one scorer correction: calls are matched to
truth genes by CODING overlap now, not span overlap. That correction changes one arm only, the
same way in both codes (VA address: baseline 27/55 · 39 · 9, `fix4` 28/55 · 32 · 3), so the
hold-out comparison stands; item 8's VA arm-1 row is the pre-correction scoring.

Verified on LRZ (2026-09-16/17, snapshots `fix5` and `fix6`, the latter adding the refinement
to the two rescue passes): reported gene ends that land exactly on the curated coordinate went
28/182 → 97/182 over APYR, DPP4 and SP (median distance 15 bp → 0 bp), with recall unchanged on
all three and precision unchanged on APYR and DPP4. SP's precision moves 0.61 → 0.59 and the
internal exon agreement of its reported models drops, because refined proteins change what later
search waves find in one of its three home loci; details in the audit doc §B4. Melittin is
byte-identical across all three codes.

**Still to do:** the write-up.


## 0a. Status after the 2026-09-10 session

**§3's blocker is resolved: cause found, fixed and pinned by a test.** Phase 0 is built. Four
findings change the plan, and each one is marked where it lands below.

1. **Cause of the nondeterminism** (§3). It is not "set iteration in general". The GOI parent
   collapse (`collapse_goi_queries_by_parent`, formerly inline in `process_region_block`)
   keeps one query per parent ID. The wavefront's expansion model
   `GOI_Melt|Apis_florea_fna_b0_l1_exon_ann` shares parent `GOI_Melt` with the home query, and
   both are **exactly 70 aa**. On a length tie the first one seen won, and "first" came from a
   set of strings, i.e. from the per-process hash seed. The set is rebuilt per search block,
   so every **block** searched after wave 2 used either the *A. mellifera* or the *A. florea*
   melittin, a coin flip per block. The prediction checks out exactly: the 3 GFFs identical
   across LRZ R1/R2 are *A. cerana* and *A. florea* (waves 1–2, before the expansion) plus the
   home. All 17 post-expansion genomes differ. Local replays reproduce R1 and R2 exactly, and
   R1 turns out to be a mixed draw across *Tetramorium*'s two blocks. Details and replay
   evidence: `docs/CAN_BRIDGE_ANALYSIS.md` §8.
2. **"The true *Tetramorium* melittin is absent in R2" was a scorer artifact.** R2 modelled
   the same locus (`OV788322.1:15,634,947-15,635,084`) as a `tandem_copy` **`gene`** feature
   (LOW); R1 modelled it as a fallback **mRNA** (AMBIGUOUS). `score_coordinates.py` read only
   mRNA, so R2 scored as a miss. Fixed. The run-to-run difference is real (a type/confidence
   flip at the correct locus), but locus_recall 0.12 vs 0.04 overstated it.
3. **The *Apis* home caps locus-level recall far below what §2 assumed.** 412 of the 724
   non-*Apis* genes sit at loci with no *Apis* member. 82 of those (APH 48, SP 30, HYAL 4)
   share an address with an *Apis* locus and stay reachable as `unresolvable_by_synteny`,
   so **330 of 724 are unreachable by construction** (tier `no_home_anchor`, below).
   **FAM_SP has zero** ant genes carrying an *Apis* SP locus label. One of the two *Apis* SP
   seeds reaches 30 ant genes, all through a shared flanking address, so SP's whole reachable
   set is `unresolvable_by_synteny` by construction; the other seed reaches none. SECP's
   second locus (30 genes) has no *Apis* member, so rung 3 is really one reachable locus.
   All three *Apis* PLA2 seeds carry one locus label that spans two chromosomes, so for PLA2
   the label is not an address.
4. **Seed lookup:** the protein FASTAs are keyed by gene-model id, not `protein=PROT_xxxx`.

**Public repo (decision 2026-09-10):** this doc carries family-level aggregates only. Seed
ids, locus labels and anchor genes are in the gitignored
`tests/benchmark_truth/family_benchmark_private.md`.

**The pre-registration is written and hashed** (§2.3). `family_recovery.py plan` produced
`tests/benchmark_truth/family_preregistration.tsv` (gitignored: unpublished curation, same
convention as the melittin truth), SHA-256
`64638c980de1e0a498e3d61972b0e7095537ace667c90ad8140a915f24c3c6f2`, built from
`EXPORT_MANIFEST.tsv` md5 `76f5c81de4befe318a88d8c096a3f8eb` and `LOCUS_ANCHORS.tsv` md5
`8533535e85eddad3ba35b507d0b06f4b` (Ant_Venoms_git @ 09225ab). Committing this doc commits
the hash, so the ceiling cannot quietly move after results come in. The ceiling per family
(one seed per family: the one with the most reachable genes; seed ids in the private file):

| family | rung | split | locus tier | address tier | unsearched home | no home anchor |
|---|---:|---|---:|---:|---:|---:|
| APYR | 1 | ladder | 31 | 0 | 0 | 1 |
| SPIN | 1 | ladder | 32 | 0 | 0 | 1 |
| DPP4 | 1 | ladder | 28 | 0 | 0 | 1 |
| HYAL | 2 | ladder | 28 | 4 | 0 | 1 |
| KAZA | 2 | ladder | 32 | 0 | 0 | 2 |
| SECP | 3 | tune | 30 | 0 | 0 | 30 |
| PLA2 | 4 | tune | 15 | 0 | 0 | 67 |
| VA | 4 | holdout | 55 | 0 | 30 | 88 |
| SP | 4 | holdout | 0 | 30 | 0 | 82 |
| APH | 5 | holdout | 31 | 48 | 0 | 57 |

**Phase 0 harness:** `scripts/benchmark/family_recovery.py` (`plan` / `stage` / `score`),
34 tests in `tests/test_family_recovery.py`. `stage` verified on the real data: 33
species-named targets plus the *Apis* home FASTA+GFF. **`score` has not yet run on a real family
run**: the first Phase 1 run is its end-to-end test.

**Loose ends (§6), checked:** every manifest scaffold exists in its genome FASTA (939/939,
after the `Probolomyrmex sp.` → `Probolomyrmex_sp` name fix that `species_key` does). GFF
bias does not arise: `stage` passes no `--target_gffs`, so all 33 targets are searched
identically, annotation or not. The uncommitted working tree still has to be committed
before Phase 1 (commands in the session hand-off).

---

## 1. Why this is worth doing

Every SynVoy claim so far rests on **one gene at a time** — melittin, oskar, STE2, decorin —
with a hand-built truth set per gene. That is n=1 evidence, four times over. It is also why
the over-calling problem has resisted three separate attempts (`docs/STATE_OF_THE_PROJECT.md`
F9): with one gene you cannot tell a systematic failure from a peculiarity of that gene.

Ivan's `final_dataset` is a **labelled, multi-family, multi-species truth set that already
exists**, was built independently of SynVoy, and is curated to the level SynVoy is judged at.

### What is actually in it (verified 2026-09-10, not quoted from memory)

`/home/faw/dev/projects/Ant_Venoms_git/data/final_dataset/`

| file | what it gives you |
|---|---|
| `gff/ALL_SPECIES.gff3` | 939 gene models, each tagged `family=`, `locus=`, `structural_status=`, `species=`, `is_venom=` |
| `gff/<Species>.gff3` | the same, split per species (35 files) |
| `proteins/<Species>.faa` | per-species proteins (35 files) — **the query seeds come from here** |
| `LOCUS_ANCHORS.tsv` | **172 loci with curated flanking anchor genes** (`anchor_up`, `anchor_dn`, `anchor_genes`, `anchor_n_species`) |
| `SPECIES_TAXONOMY.tsv` | 35 species: 32 ant ingroup + 3 outgroup; 10 subfamilies |
| `../genomes/` | 14 GB, 41 assembly dirs, `ACCESSIONS.tsv` maps species → accession; 23 carry a GFF |

**Excluding GR1/GR2/GR3/GR7 as asked** (192 genes — those are the paper's showcase regions and
are already benchmarked), the usable test bed is **747 genes across 12 families**:

| family | genes | species | distinct loci | shape |
|---|---:|---:|---:|---|
| `FAM_VA` | 175 | 34 | 22 | multi-copy (venom allergen / CAP) |
| `FAM_APH` | 137 | 33 | 10 | **tandem array** (acid phosphatase) |
| `FAM_SP` | 114 | 33 | 16 | multi-copy (serine protease) |
| `FAM_PLA2` | 85 | 33 | 14 | multi-copy, most fragmented |
| `FAM_SECP` | 61 | 34 | 6 | two main loci |
| `FAM_KAZA` | 35 | 34 | 3 | near single-copy |
| `FAM_HYAL` | 34 | 34 | 3 | near single-copy |
| `FAM_SPIN` | 34 | 34 | 2 | **single-copy** |
| `FAM_APYR` | 33 | 33 | 2 | **single-copy** |
| `FAM_DPP4` | 30 | 30 | 2 | **single-copy** |
| `FAM_VA2` / `FAM_VA4` | 7 / 2 | 2 / 1 | — | too small, skip |

Nearly everything is `complete-functional` (e.g. APH 137/137, SECP 61/61, HYAL 34/34), so
failures will be SynVoy's, not the truth set's.

---

## 2. The design that makes this diagnostic rather than just a score

### 2.1 Two levels of truth — and they are not the same question

The `locus=` attribute is the orthology grouping; `family=` is only family membership.

* **Family-level recall** — did SynVoy find *a* member of the family in species X?
* **Locus-level recall** — did it find the member at the *same `locus=` label*, i.e. the
  actual ortholog?

This distinction **is** the over-calling problem, now with labels. A call that is
family-correct but locus-wrong is exactly what the AMBIGUOUS tier and the F9 phylogenetic
placement check were built to catch, and here you can finally measure whether they do.

### 2.2 A difficulty ladder, so a failure localises itself

1. **Single-copy, one locus** — APYR, SPIN, DPP4: one main locus each, 29–33 species.
   Clean 1:1 orthology across 30+ ant species. **If SynVoy fails here it is a bug**, full stop.
2. **Single-copy + one rare second locus** — HYAL (main locus 29 sp, second locus 4); KAZA.
3. **Two real loci** — SECP (31 + 26 species). From an *Apis* home only the first is
   reachable (the second has no *Apis* member; §0). So this rung tests "one reachable
   locus plus a same-family decoy locus", which is still worth having for precision.
4. **Multi-locus** — VA (four main loci, 26–33 species each), SP, PLA2.
5. **Tandem array** — APH, three tandem loci at one address.

(Locus labels per rung: `tests/benchmark_truth/family_benchmark_private.md`.)

### 2.3 Pre-register what *must* fail — this is the anti-overfitting device

Frank's worry ("changing stuff would mean overfitting it") is the right one, so pin down
before running which failures are **correct behaviour**.

`LOCUS_ANCHORS.tsv` gives three anchor addresses shared by more than one locus. Synteny
**cannot** separate loci that sit at the same address — not as a limitation of SynVoy, but by
definition:

| family | loci sharing one address |
|---|---:|
| `FAM_APH` | 3 (the tandem array) |
| `FAM_HYAL` | 2 |
| `FAM_SP` | 2 (one ant locus, one *Apis* locus) |

The curator's own note on the APH tandem loci says it outright: they share the middle
locus's address. (Anchor genes and locus names: private file.)

**Score those as `family-correct / locus-unresolvable-by-synteny`, not as errors.** Any change
that "fixes" them is overfitting to a distinction the method has no signal for. This is the
line to hold. (The honest ceiling: locus-level recall on APH is capped somewhere near 1/3 of
its tandem calls unless a non-synteny signal — architecture, phylogenetic placement — does the
separating.)

**The larger pre-registered failure class (added 2026-09-10): loci the home genome does not
have.** SynVoy anchors on the *home* genome's locus. A locus with no *Apis* member, and not at
an *Apis* locus's address, has no neighbourhood to anchor on, so it is unreachable by synteny
from this home **by construction**. That covers 330 of 724 non-*Apis* genes, even with
every *Apis* locus of the family searched. It is tier
`no_home_anchor` in the harness, excluded from locus recall, and a hit there is reported as
`found_without_anchor`, never counted as a win. The four tiers are defined in
`scripts/benchmark/family_recovery.py`, and the per-family ceiling is in §0.

### 2.4 Attribute every failure to a layer before touching a parameter

SynVoy already logs enough to separate these; do not skip to tuning.

| layer | question | evidence in the run |
|---|---|---|
| **A. neighbourhood** | did any block land on the true locus? | `regions/*.bed`, `Found N discrete syntenic blocks`, compare to `anchor_up`/`anchor_dn` |
| **B. modelling** | block was right, gene not modelled? | region GFF has flanking but no `SynVoyRole=goi` there |
| **C. classification** | modelled but confidence too low? | `Confidence=`, `InferenceReason=` on the model |
| **D. adjudication** | right gene, wrong label? | ownership / dedup / `paralog_not_goi` |

A **layer-A** failure is a search bug (and `LOCUS_ANCHORS.tsv` tells you which anchors *should*
have been hit — a far sharper probe than melittin ever gave). **Layer C** is where a parameter
argument is legitimate. **Layer D** is the known-weak adjudication layer.

---

## 3. Blocking prerequisite — read this before running anything

> **RESOLVED 2026-09-10. See §0.1.** The cause was the GOI representative collapse. It is
> fixed in code, pinned by `tests/test_a0_determinism.py` (12 hash seeds), and
> `PYTHONHASHSEED=0` is now set in `nextflow.config` as defence in depth. The staged
> `hashseed_probe` is superseded; do not submit it. The right LRZ confirmation now is one
> post-fix replicate pair (§3, step 3). The original text follows, with one correction.

**The pipeline is not reproducible.** Proven 2026-09-09, LRZ job 5778368: two runs with
byte-identical flags differed in **17 of 20 per-genome GFFs**; melittin `locus_recall` came out
**0.12 in one and 0.04 in the other**. ~~The true *Tetramorium* melittin was modelled in one
replicate and absent in the other.~~ *Corrected 2026-09-10:* it was modelled in both, as an
AMBIGUOUS fallback mRNA in R1 and a LOW tandem-copy `gene` in R2, which the mRNA-only scorer
missed. Full write-up: `docs/CAN_BRIDGE_ANALYSIS.md` §8.

Deterministic through block construction (identical hits, identical blocks); diverges at
annotation. **Prime suspect:** `bin/iterative_search_runner.py:3225` iterates `unique_queries`,
a **set of strings**, into an order-dependent `found_queries`; `PYTHONHASHSEED` is set nowhere
in the repo, config, or Dockerfile.

With 10 families × ~34 species the noise would swamp every comparison, and worse, would make
"I changed a parameter and it improved" unfalsifiable — the exact overfitting trap.

**Do this first:**
1. ~~Submit the staged probe~~. Superseded 2026-09-10: the cause was reproduced and
   attributed locally (reduced-wavefront replays under pinned seeds, `CAN_BRIDGE_ANALYSIS.md` §8).
2. ✅ Done 2026-09-10: set iterations sorted at the functional sites
   (`iterative_search_runner.py`: `unique_queries`, the GOI-proxy sweep; cosmetic ones in the
   same file and `plot_synteny.py`), all of `bin/` audited with an AST pass,
   `PYTHONHASHSEED` pinned in `nextflow.config` `env {}` (it reaches container tasks too), and
   regression tests added.
3. **Still open:** one LRZ replicate pair of the 19-target melittin benchmark on the fixed,
   committed code: `scripts/lrz/determinism_postfix.sbatch`. It runs from the pulled repo,
   and **R2 uses `PYTHONHASHSEED=1`** via a `-c` override. With the seed pinned in
   `nextflow.config`, a plain pair would share seed 0 and pass vacuously. The payload refuses
   to run on a checkout without the fix, prints the seed each task actually received, and scores
   with the fixed `score_coordinates.py`. Expect 20/20 identical GFFs. Anything less means a
   second source, and the family benchmark then needs ≥3 replicates per number.

---

## 4. Proposed execution

**Home genome — settled 2026-09-10: use *Apis mellifera*.** Verified that all ten target
families have an Apis member, so one home genome covers the whole benchmark: same home, same
34 targets, only the query changes. It is also the honest test — Apis is the outgroup to all
32 ant ingroup species, so every recovery crosses the bee/ant split rather than staying inside
a subfamily. `GCF_003254395.2` (Amel_HAv3.1) is chromosome-level with a GFF, which is what
`EXTRACT_FLANKING` needs.

### Phase 0 — harness (local, no cluster)
Build `scripts/benchmark/family_recovery.py`:
* read `ALL_SPECIES.gff3` → truth table `(family, locus, species, scaffold, start, end)`
* pick seed protein from `proteins/<Species>.faa` by its `>GM_xxxxxxx` header (*corrected
  2026-09-10*: the FASTAs are keyed by gene-model id, not `protein=PROT_xxxx`)
* score a SynVoy run at **both** levels (§2.1) by coordinate overlap, reusing
  `scripts/benchmark/score_coordinates.py` (already written, already CI-frozen)
* emit the layer-A–D attribution (§2.4) per miss
* mark shared-anchor loci `unresolvable` (§2.3)

Unit-test it against a **synthetic** truth set, as `tests/test_coordinate_benchmark.py` does,
so CI does not depend on the 14 GB genome set.

### Phase 1 — the easy ladder, one seed per family
`SPIN`, `APYR`, `DPP4` (single-copy) + `HYAL`, `KAZA`. **~5 runs.** Expect near-total recall.
Anything below that is a bug worth more than the rest of this plan.

### Phase 2 — the hard ladder
`SECP`, `VA`, `SP`, `PLA2`, `APH`. **~5 runs.** This is where family-vs-locus recall will
split, and where the AMBIGUOUS tier and F9 phylo placement get their first real test.

### Phase 3 — does the seed matter?
On families that failed, re-seed from 2–3 species at different phylogenetic distances
(the 3 outgroups vs Myrmicinae ingroup). If recovery depends strongly on the seed, that is a
**property of the method** worth reporting, not a parameter to tune away.

### Phase 4 — change one thing, retry
Only now. Whatever is changed must be justified by a layer attribution from §2.4, and:

> **Split the families: tune on `SECP` + `PLA2`, hold out `VA` + `SP` + `APH` untouched until
> the end.** If a change helps the tuning families and not the held-out ones, it was
> overfitting. This is the only mechanism here that can actually answer Frank's question, and
> it only works if the held-out set is genuinely not looked at.

**Cost:** the 19-target melittin run took ~38 min. A ~34-target run is roughly 60–90 min, so
Phase 1+2 is ~10 runs ≈ 10–15 h of cluster time — comfortably one or two batch jobs. Multiply
by 3 if replicates are still needed. Keep an eye on the Ant_Venoms jobs sharing the cluster.

---

## 5. What a good outcome looks like

Not "SynVoy scores X". The deliverables are:

1. **A per-layer failure census** across 10 independent families — the first evidence that can
   distinguish a systematic SynVoy defect from gene-specific bad luck.
2. **A defensible ceiling**: which losses are unreachable by synteny (§2.3), stated up front
   rather than discovered after a tuning session.
3. **A verdict on the AMBIGUOUS tier and F9**, measured against 700+ labelled genes instead of
   19 hand-checked ones.
4. **A held-out result** that says whether any change generalises.

A finding that "SynVoy recovers single-copy venom families across 30 ant species but cannot
resolve tandem arrays, and here is the layer where it breaks" is a far stronger paper claim
than the current melittin headline — and it is honest about the ceiling.

---

## 6. Loose ends to check first

* **Species/scaffold naming** must match between Ivan's GFF3 `seqid`s and the genome FASTAs in
  `data/genomes/` — the melittin benchmark already lost time to a CM-vs-NC accession mismatch
  (`scripts/benchmark/build_melittin_loci.py` carries the alias map). Verify before Phase 1.
* Only **23 of 40** assembly dirs carry a GFF. Targets without one lose native-annotation
  borrowing; check whether that biases the ladder.
* SynVoy's uncommitted working tree (two-sided bridging, AMBIGUOUS tier, F9 phylo placement)
  should be committed, or at minimum pinned by hash, before a 10-family benchmark runs against it.
