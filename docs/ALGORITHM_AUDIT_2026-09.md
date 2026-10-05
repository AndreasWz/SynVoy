# SynVoy algorithm audit — September 2026

Scope: the GOI search core (`bin/iterative_search_runner.py`, `bin/smith_waterman_search.py`,
`bin/annotate_goi_exons.py`) and the steps that consume its output, traced on the
2026-09-10/11 family-benchmark runs (`docs/NEXT_SESSION_FAMILY_BENCHMARK.md`) and the
post-determinism-fix melittin run. Every number below was measured, not estimated: either
from the runs' kept logs/GFFs, or by re-running a pipeline step with the pipeline's own
functions and the run's own `script.sh` settings on real genome windows. Committed text is
aggregate-only; loci and identifiers stay in the gitignored benchmark notes.

## 1. What the search core does, step by step

| Step | Code | Decision it makes |
|---|---|---|
| Stage 1: map flanking + GOI proteins | `run_mmseqs_easy_search_with_retries` | MMseqs2 translated search of the whole target genome (ORFs >= 30 codons, tantan masking, composition-bias correction, `--max-seqs 300`), E <= `evalue` |
| Blocks | `identify_synteny_blocks`, `merge_synteny_blocks` | cluster flanking hit loci (< `cluster_distance`), bridge collinear gaps, keep blocks with >= `min_block_genes` distinct flanking genes, cap at `max_blocks_per_genome` |
| Window | `process_region_block` | block +- adaptive padding (50-200 kb) |
| GOI hits in window | `run_augmented_search` | MMseqs2 on query fragments + tblastn + Smith-Waterman, merged by a greedy spatial dedup that keeps the higher bit score |
| GOI loci | `split_hits_into_loci`, `_locus_rank` | keep the 2 best hit loci per window |
| Models | `detect_tandem_duplications` / miniprot (`annotate_using_miniprot`) / hit-chain fallback | tandem copies from raw hits; else a miniprot model of the region; else a chain of raw hits |
| Evidence class | `_classify_goi_evidence` | HIGH/MEDIUM/AMBIGUOUS/LOW from identity, query coverage, exon count, and the window's count of distinct flanking genes |
| Seeds, tree | `_classify_goi_for_seed_and_tree` | HIGH/MEDIUM confident/probable GOI proteins seed the next wave and, with tandem copies, feed the tree, the paralog RBH panel and the phylo-placement check |

Two properties matter for everything below:

* **Miniprot is specific; the hit-derived evidence types are not.** On 400 runs over random
  100-kb windows of a real ant genome (five family queries, shuffled and real), miniprot in
  the pipeline's sensitive mode produced no model at all (positive control on the true
  window: the full 8-exon model). The off-truth miniprot calls in the benchmark were real
  homology — orthologs at loci the curation lacks, or short fragments of family members;
  every chance call traced came from the evidence types built directly from search hits
  (`tandem_copy`, `fallback_hit_span`, `rescued_exon`, `raw_hit`).
* **Every on-truth HIGH/MEDIUM call is a miniprot model** reached by MMseqs2 or tblastn:
  149/149 in HYAL, KAZA, APYR and DPP4 still have a non-SW hit on them when the window
  search is re-run with Smith-Waterman switched off.

## 2. Findings

Severity: **A** wrong results that reach the output; **B** wrong intermediate values that
bias decisions; **C** method limitation or latent defect.

### A1. Smith-Waterman hits had fake statistics — fixed

`smith_waterman_search.py` wrote E = 0.001 for every hit and the raw alignment score in the
bit-score column. SW reports the best local alignment in each of the six reading frames of
every window, homolog or not.

* 100 % of windows got SW hits, 6.0 per window, in all eight benchmark runs (e.g. 5,942 SW
  hits in 992 KAZA windows; after the dedup 5,964 hits kept of 6,773).
* Null test (composition-preserving shuffled queries, 192 real ant-genome windows of
  100-300 kb, the pipeline's own script and filters): median true E 32-93 depending on the
  query; 36-77 % of these chance hits span >= 35 % of a short query (the tandem detector's
  bar); 90-100 % pass the short-query fallback gate, which reads "bits >= 30" off the raw
  score.
* Raw scores are about twice a bit score, so SW hits won the spatial dedup against real
  tblastn/MMseqs2 hits of the same gene (seen in HYAL on-truth windows, where the gene's
  tblastn hit disappeared once SW was on).
* The gap costs were parasail 10/1, i.e. BLAST 9/1 (parasail charges `open` for a gap's
  first residue: verified, a 1-residue gap costs `open`), a setting with low relative
  entropy (H = 0.052; BLAST 11/1: H = 0.14) where chance alignments are long.

Fix (decision 2026-09-11): BLAST-standard BLOSUM62 11/1 (parasail 12/1) and Karlin-Altschul
bit score and E-value for the window's six-frame search space (lambda 0.267, K 0.041, as
printed by blastp 2.17.0). The runner's existing E <= 1 filter for SW hits now works:
shuffled-query hits surviving it fell from 6.0 to 0.07 per window. The analytic E is
conservative here (stop codons cut a six-frame translation into short open frames; 2 % of
per-frame chance maxima reach E <= 1 where a calibrated E would give ~15 %), which errs on the
side of fewer false hits. What it costs, on the melittin truth set: SW remains the only method
that reaches five weak members of a Cardiocondyla tandem array (E 0.0027-0.25, all LOW calls
before); the curated Tetramorium melittin-family gene (E 18) and one Formica fragment (E 23)
lose their SW pointer. Both were LOW calls, and at those E-values the alignment is not evidence
of homology on its own (Karlin & Altschul 1990; Pearson 1998).

### A2. Model proteins were translated in the wrong frame — fixed

`annotate_using_miniprot` translated each CDS from its first base and trimmed the GFF3 phase
only on the first CDS; the runner built every model protein by concatenating those per-exon
translations and discarded the correct `full_protein`. Every exon behind a phase-1/2 intron
was read in the wrong frame. Coordinates, Identity and confidence were right, so nothing in
the GFF showed it — but the protein is what the pipeline uses next:

* every multi-exon GOI model with a phase-1/2 intron was affected, which in the family
  benchmark was every multi-exon model: 32/32 HIGH in APYR, 15/15 HIGH + 16/16 MEDIUM in
  DPP4, 33/33 HIGH in HYAL, 3/3 HIGH in KAZA, 30/30 HIGH in PLA2 (the one multi-exon
  melittin model has a phase-0 intron and was not); for the long genes 50-79 % of the
  coding sequence was emitted out of frame;
* a 64 %-identity hyaluronidase model: SW score against the query 1,491 as written vs 2,400
  for miniprot's own translation (47.6 % vs 64.1 % identity);
* consumers: `regions/*.faa` -> next-wave seeds (`expanded_db.faa`), `goi_for_tree.faa`
  (MAFFT/IQ-TREE), the section-1m paralog RBH ownership check, the F9 phylo-placement check,
  and flanking and rearranged-flanking models (same construction).

Concatenating CDS DNA (the discarded `full_protein`) is not correct either: it runs out of
frame after the first frameshift miniprot models (`Frameshift=1`). Fix: the model protein is
miniprot's own translation (`--trans`, the `##STA` record); each exon's sequence honours its
own phase; CDS lines in the GFF carry the real phase (they were always `0`, a GFF3 violation
that breaks any downstream CDS extraction). Tests: `tests/test_miniprot_model_translation.py`
(incl. a real miniprot run on a gene with phase-1 and phase-2 introns).

**Consequence:** the earlier adjudication verdicts — "F9 phylo placement uninformative",
paralog ownership results — were measured on these proteins for every multi-exon gene.
Re-measured on the fixed code in section 3: F9 is no better on correct proteins.

### B1. GOI precision defects found on the family benchmark — fixed

Documented in `CLAUDE.md` known bug 19: tandem copies promoted on sequence alone; "exons"
chained from two alignments of the same genome; fallback coverage measured as the query span;
fragment models lifted to MEDIUM by one flanking gene. Two more of the same kind, found in
this audit:

* **Fallback identity was the unweighted mean of per-hit identities.** A 12-aa window at
  80 % and a 100-aa alignment at 30 % were reported as a 55 % model; identity is identical
  columns over aligned columns (35 %). Now length-weighted (`_pooled_identity`), for GOI and
  flanking hit-chain models.
* **`_locus_rank` ranked a window's GOI hit loci by query span**, the bug fixed in the
  fallback gate: two short chance hits from opposite query ends outrank a real locus for one
  of the two slots kept per window. Now union coverage, also in the strict raw-hit fallback.

### B2. Identity decides where significance should — open

The classifier's sequence bars are identity and coverage thresholds that do not depend on
alignment length, and two selection steps rank by identity first (GOI candidate selection in
`process_region_block`; `_deduplicate_hits` inside the tandem detector, which also takes a
copy's identity from the lowest-E hit). Identity alone is not a measure of homology for short
alignments: the identity above which an alignment reliably implies homology rises steeply
below ~80 aligned residues (about 27-35 % at 70 residues and 48-63 % at 25, by the curves of
Sander & Schneider 1991 and Rost 1999). The
short-peptide preset's MEDIUM bar is 25 % identity. For miniprot models this is mitigated
(miniprot's own scoring is specific, section 1); for hit-derived evidence it is not, and
before A1 the hits themselves carried no usable significance. Proposal, not implemented:
gate hit-derived evidence on the hit E-value (now real for all three searches) and rank
candidates by bit score. Needs a tuning round before the hold-outs.

### B3. Synteny support has no chance model — measured, not gated

`BlockFlankingSupport` is the number of distinct flanking genes with any hit in the window
(block +- up to 200 kb). A block needs only 2 flanking genes within 150 kb, and up to 80
windows per genome are searched; nothing asks how often that happens by chance given how many
places each flanking gene hits genome-wide (gene families hit many). Gene-cluster statistics
exist for exactly this (Durand & Sankoff 2003; Hampson et al. 2003, 2005). Decision
2026-09-11: measure first — block E-values (Poisson-binomial tail of the distinct-gene count
under uniform placement, scaled by the number of windows) for every window of the local
family runs, compared with the raw count as a separator of on- and off-truth calls. Result
(section 3): AUC 0.91 vs 0.88, and the calls it would catch are already demoted by B1.

### B4. Gene-model boundaries vs curated annotations — termini fixed, exon losses open

A call that overlaps the right gene can still be the wrong gene model. Measured on the fixed code
by comparing, exon by exon, every recovered gene's model against the curated annotation of the
same gene in eight family arms (216 genes, 1,439 curated exons):

| | exons | share |
|---|---:|---:|
| both boundaries exact | 901 | 62.6 % |
| one boundary exact | 351 | 24.4 % |
| overlap only | 11 | 0.8 % |
| missing from the model | 176 | 12.2 % |
| model exons with no curated counterpart | 147 | — |

Internal splice sites are reproduced exactly wherever the model spans them — in one 11-exon gene
all ten internal junctions are byte-identical and only the two termini differ. The disagreements
are concentrated in four mechanisms, each read by hand on 12 cases (details, loci and identifiers
in the gitignored benchmark notes):

1. **Termini are short by a fixed amount** — the most common deviation. The model starts 3 codons
   inside the annotation (so it lacks the initiator Met) and ends 2 codons plus the stop codon
   early, because miniprot reports only what it aligns. `annotate_goi_exons.py` already computes
   `has_start_codon` / `has_stop_codon` per exon and `iterative_search_runner.py` writes them onto
   CDS lines when true — which is rarely: across seven arms and ~200 models `StartCodon` appears 11
   times and `StopCodon` 121. They are never aggregated to the model, and `_model_status` labels a
   model `complete` from query coverage and exon count alone, so a model that begins three codons
   inside the gene with no initiator Met is still reported `complete`, at HIGH confidence.
   **Fixed** (`refine_model_termini`): after miniprot, the terminal CDS is walked outwards in
   frame to the nearest upstream ATG and the first downstream stop codon, and `StartCodon` /
   `StopCodon` now appear on the model's own mRNA line. The walk refuses to invent sequence — it
   stops at a splice site (`AG` before the first coding base, `GT` after the last), will not cross
   an in-frame stop, is bounded by the unaligned query residues, and leaves a terminus it cannot
   resolve exactly as miniprot reported it. Verified on four real *Apis* genes with the query
   trimmed at both ends: three recovered the curated boundaries and the curated protein exactly
   (an 11-exon apyrase byte-for-byte, on the minus strand), and in the fourth a chance `GT` inside
   a valine codon blocked a legitimate 5-codon walk — the guard's price. An untrimmed query, whose
   model already reaches both ends, is left untouched. Tests: `tests/test_model_termini.py`.
2. **Small terminal exons across a long first or last intron are lost** — most of the 176 missing
   exons. In one family a 22-bp first exon 1.4 kb away is missed in nearly every species (31 of
   its 95 curated exons), which costs 6 % of the coding sequence but 60-70 % of the reported gene
   *span*; in another, two 5' micro-exons of 15 and 58 bp are missed in most species; in a third,
   two 3' exons behind a 5.1-kb intron are missed, costing 145 of 402 residues. The pipeline
   already runs a relaxed rescue pass, but it emits scattered separate LOW records instead of
   completing the model it belongs to. Fixable with moderate effort: when more than ~10 query
   residues at a terminus are unaligned, re-align that fragment against the flanking window and
   attach it to the model only with a valid splice junction.
3. **Distant terminal exons placed where they do not belong** — the mirror image: a spurious
   ~70-bp first exon 4.5 kb upstream inside a cluster of family fragments, or a genuine 26-bp
   first exon replaced by a 35-bp candidate 5.7 kb further out. Both are terminal, small, and sit
   behind a gap 10-50x the model's own median internal intron, which is the signal to flag them.
   Note what cannot work: capping miniprot's `-G` globally. The home genes' largest introns here
   are 118-1,636 bp while real introns at the same loci in the target species reach 1.4, 5.1 and
   10.2 kb, so a cap tight enough to reject the spurious exons would also reject the true ones.
   `annotate_exons_from_hit_list` currently sets `-G` to the whole search sub-region
   (`max(20000, region_end - region_start)`), i.e. no effective limit, while the adaptive
   home-derived intron limit the pipeline computes for block seeding is not passed to it at all.
   A per-model relative rule, plus the cross-species prior that the same gene is modelled in ~30
   genomes, is the defensible version.
4. **A sprawling rescue model can outrank the compact correct one.** In one case the adjudicated
   record is an 11.4-kb hull-rescue model at HIGH while the 2.6-kb in-block model that matches the
   curated gene exactly is MEDIUM; their overlap is 23 %, below the report's 80 % dedup bar, so
   both survive and the wider one wins. **Fixed** in `dedupe_goi_annotations`: records are now
   also clustered when one span sits wholly inside the other, and the contained model represents
   the pair when the wider one is at least twice as long while aligning less than 15 percentage
   points more of the query — the fused-model signature (here 8.8 kb of extra span for 8 points of
   coverage). It refuses to decide when either coverage is unrecorded, so a genuine long gene still
   beats a fragment of itself. The kept record carries `superseded_wide_model`, counted in the
   report as `wide_models_superseded`.

A fifth cause is not a pipeline defect at all. In one family every model in all 33 genomes spans
4-7x the curated gene because the *seed protein* is a predicted long isoform, twice the length of
the curated venom protein and of the family's members, whose second half has no other homolog in
the home proteome and is not homologous to its own first half. The pipeline modelled faithfully
what it was asked to find. The lesson is a query-QC step (compare the seed's length with the
home gene's other isoforms and with the family, and say so loudly), not a change to the search.

**Verified end to end on LRZ (2026-09-16/17).** `fix5` = the terminal refinement and the
fused-model rule; `fix6` = `fix5` plus the same refinement on the two rescue passes, which call
miniprot directly and had bypassed it — for 8 of 31 recovered apyrase genes the reported
coordinates came from an unrefined hull-rescue model. Three families rerun, 91 recovered genes,
reported gene ends compared with the curated ones:

| | `fix4` | `fix5` | `fix6` |
|---|---:|---:|---:|
| gene ends exactly on the curated coordinate | 28 / 182 (15 %) | 78 / 182 (43 %) | 97 / 182 (53 %) |
| genes exact at both ends | 1 | 18 | 21 |
| median distance of a gene end from the curated one | 15 bp | 9 bp | 0 bp |
| miniprot models with start and stop codon (APYR · DPP4 · SP) | 0 · 0 · 0 | 16 · 27 · 26 | 16 · 27 · 26 |

Recall and precision are unchanged on APYR (30/31) and DPP4 (28/28). Per family, exact both-end
genes went 0 → 14 (APYR), 0 → 0 (DPP4, whose other end is the lost 5' micro-exons of item 2), and
1 → 7 (SP). One cost, measured rather than assumed: on SP the refined model proteins seed later
search waves, and one of its three home loci then built 81 models instead of 44 (the other two
loci, the home GOI model and `goi_info.json` are byte-identical between codes). Recall stays
29/30, precision moves 0.61 → 0.59, and the internal exon structure of the reported SP models
agrees less often with the curated one (69.8 % → 62.0 % of curated exons exact) even as their ends
move closer (median 30 → 3 bp). The melittin benchmark is byte-identical across `fix4`, `fix5` and
`fix6` on every coordinate metric: one of its eight GOI models is a miniprot model, and that one
already had its start codon.

Consequences to keep in mind when reading any of this benchmark's numbers: recall and precision
are scored by overlap, so they are insensitive to all of the above. Reported gene *spans*, model
proteins and anything derived from them (trees, paralog checks) are not.

One of the 228 pairings shared no coding sequence with the gene it was credited to — a model on
the opposite strand inside a 10-kb first intron — because the scorer matched calls to truth genes
by span overlap. The harness now requires **coding** overlap (`load_truth_cds`, `_ov`), falling
back to span overlap only when the curated GFF3s are unavailable. Rescoring every arm with the
corrected rule changes exactly one arm, and identically in both codes, so no comparison in this
document moves: the VA address arm reads 27/55 · 39 calls · 9 off-truth for the baseline and
28/55 · 32 · 3 for `fix4` (precision 0.77 → 0.91, previously reported as 28/55 vs 29/55 and
0.79 → 0.94). Every other arm is unchanged.

### C1. Phylo-placement check (F9) — limitations

Gene divergence is 1 - identity from a *local* alignment (short high-identity windows read as
close), p-distance saturates for divergent genes, species distance is a taxonomic rank
difference rather than time, and the concordance test assumes comparable rates across
lineages — the assumption venom genes under diversifying selection violate (3 of its 6
"discordant" calls in the benchmark were real orthologs). Advisory only; re-measure after A2.

### C2. Smaller items

* `rescue_goi_hull._translate_model` concatenates CDS DNA with a first-phase trim: correct
  unless miniprot modelled a frameshift. Minor (fallback path); not changed.
* Stage 1 is MMseqs2 translated search with ORFs >= 30 codons: a flanking gene's micro-exons
  are found only inside a longer stop-to-stop frame. The GOI window search adds tblastn and SW,
  so this affects anchors, not the GOI.
* `smith_waterman_search.py --method ssearch36` silently runs the in-process parasail loop when
  parasail is installed, and `ssearch36` itself is not a translated search. Unused by the
  pipeline (auto -> parasail); that loop now also writes real statistics.
* The miniprot wrapper wrote its temp files as `/tmp/synvoy_mp_*_<pid>`; pids repeat across
  containers sharing `/tmp`, which would swap models silently. Now `mkstemp`.
* `batch_rbh_check` has no caller; it read a module-global `args` that exists only after
  `main()` (fixed), and its test fails once `mmseqs` is on PATH because its 6-aa toy peptides
  are too short to search. Always skipped before (no `mmseqs` on the test PATH).
* `cluster_grs.estimate_pvalue` is a 200-draw label-permutation test (floor p = 0.005, no
  multiplicity correction); it only orders plot regions.

## 3. Measurements after the fixes

**Synteny block E-value (B3).** Stage 1 re-run with the runner's own functions for the five
local families reproduces the runs' windows (224/224, 352/352, 433 vs 431, 998 vs 992,
1227 vs 1226). For every adjudicated record, the E-value of the window it sits in:

| records | E < 1e-6 | 1e-6 to 1e-2 | 1e-2 to 1 | E >= 1 |
|---|---:|---:|---:|---:|
| on-truth HIGH/MEDIUM (154) | 144 | 6 | 3 | 1 |
| off-truth MEDIUM (46) | 2 | 1 | 3 | 40 |
| off-truth AMBIGUOUS (94) | 27 | 34 | 26 | 7 |
| off-truth HIGH (4) | 4 | 0 | 0 | 0 |

Separation of on- from off-truth records (AUC): block E-value 0.910, raw flanking count
0.884, collinear run length 0.886. The off-truth MEDIUM calls sit in chance-level blocks —
they are the tandem-copy and fragment classes the B1 fixes already demote — and the rest of
the off-truth records sit in genuinely syntenic blocks, where no synteny statistic can help.
Not gated; a candidate for a reported attribute.

**Family benchmark on the fixed code** (LRZ snapshot `fix4`; per-arm table in
`docs/NEXT_SESSION_FAMILY_BENCHMARK.md` §0 item 7): recall unchanged on every arm the baseline
could reach (DPP4: the same 28/28 and the same three uncurated orthologs); off-truth final calls 3 → 0 (HYAL), 1 → 0 (KAZA), 27 → 1 (SPIN, wrong-address
seed); the PLA2 address seed's 18 "off-truth" calls are reciprocal best hits of the seed gene,
one per species at the address (an uncurated ortholog group). Melittin: adjudicated records
95 → 10, MEDIUM-or-better 7 → 4 (three chance tandem calls removed, one weak true member
gained, one true copy MEDIUM → AMBIGUOUS whose MEDIUM had come from a chance SW "tandem
partner"), false positives at gene-lost control loci 2 → 0.

**Hold-out confirmation (VA, SP, APH; baseline vs `fix4`, 2026-09-13).** Five arms, each run
once per code, submitted before any hold-out number was seen and with no parameter tuned on
these families (table in `docs/NEXT_SESSION_FAMILY_BENCHMARK.md` §0 item 8). Recall never
dropped: +1 gene on the VA address arm, equal on the other three reachable arms, and the
pre-registered wrong-address arm failed as predicted. Precision rose where the B1 fixes apply
(VA address arm: 8 off-truth calls → 2, 0.79 → 0.94; the 5 chance MEDIUM calls on 36-132 bp
windows are gone). The A2 frame fix shows its downstream effect here: because the model proteins
that seed the next wave are now correct, loci the baseline could only reach through the hull
rescue are modelled by the main in-block path, e.g. one SP locus goes from `Identity=54.9
MEDIUM probable_goi` to `Identity=74.1 HIGH confident_goi` on the same miniprot path. Median
on-truth call identity: SP 53 → 69 %, APH 59 → 72 % (VA address unchanged at 58 → 59 %, mostly
single-exon models). Cost: 3 new off-truth MEDIUM calls on SP at 40-52 % identity.

**Positional adjacency (measured, not gated).** Positional orthology asks more than the flanking
count: does the call sit next to the target copy of one of the GOI's nearest home neighbours
(the 2 nearest on each side, nearest flanking model within 100 kb)? On the fixed code, all 120
on-truth HIGH/MEDIUM records and all 31 on-truth AMBIGUOUS records are adjacent, as are the 19
real but uncurated orthologs; 14 of the 21 off-truth AMBIGUOUS records (the other 7 are SPIN
cluster paralogs) and the one remaining off-truth MEDIUM in a wrong-address run are not. A strict "between a left and a right
neighbour" rule fails where the GOI sits at a synteny breakpoint (29 of 31 APYR orthologs are
adjacent to one side only). What adjacency cannot do is the thing AMBIGUOUS exists for: at the
melittin locus the Bombus records (gene lost) sit ~7 kb from a home neighbour, as close as the
real Euglossa copy (3.8 kb), and in SPIN the off-truth calls are paralogs in the same cluster at
the address. It separates off-position chance hits from calls at the locus, not presence from
loss, so it does not justify promoting AMBIGUOUS; demoting non-adjacent calls would cost real
orthologs in rearranged neighbourhoods (the decorin case sits 1-2 Mb from its flanking genes).
A candidate reported attribute, not a gate.

**F9 and ownership on correct proteins (A2).** F9 does not improve: across HYAL, APYR and PLA2
all 6 `phylo_discordant` records are real orthologs and no off-truth record is flagged. The
paralog-ownership check on the SPIN address arm demoted one off-truth HIGH (right) and two
on-truth calls (wrong).

## 4. References

* Altschul SF, Gish W (1996) Local alignment statistics. Methods Enzymol 266:460-480.
* Altschul SF, Bundschuh R, Olsen R, Hwa T (2001) The estimation of statistical parameters for
  local alignment score distributions. Nucleic Acids Res 29:351-361.
* Dalquen DA, Dessimoz C (2013) Bidirectional best hits miss many orthologs in
  duplication-rich clades such as plants and animals. Genome Biol Evol 5:1800-1806.
* Durand D, Sankoff D (2003) Tests for gene clustering. J Comput Biol 10:453-482.
* Frith MC (2011) A new repeat-masking method enables specific detection of homologous
  sequences. Nucleic Acids Res 39:e23.
* Hampson S, McLysaght A, Gaut B, Baldi P (2003) LineUp: statistical detection of chromosomal
  homology with application to plant comparative genomics. Genome Res 13:999-1010.
* Hampson SE, Gaut BS, Baldi P (2005) Statistical detection of chromosomal homology using
  shared-gene density alone. Bioinformatics 21:1339-1348.
* Karlin S, Altschul SF (1990) Methods for assessing the statistical significance of molecular
  sequence features by using general scoring schemes. PNAS 87:2264-2268.
* Li H (2023) Protein-to-genome alignment with miniprot. Bioinformatics 39:btad014.
* Pearson WR (1998) Empirical statistical estimates for sequence similarity searches.
  J Mol Biol 276:71-84.
* Rost B (1999) Twilight zone of protein sequence alignments. Protein Eng 12:85-94.
* Sander C, Schneider R (1991) Database of homology-derived protein structures and the
  structural meaning of sequence alignment. Proteins 9:56-68.
* Schäffer AA et al. (2001) Improving the accuracy of PSI-BLAST protein database searches with
  composition-based statistics and other refinements. Nucleic Acids Res 29:2994-3005.
* Smith TF, Waterman MS (1981) Identification of common molecular subsequences.
  J Mol Biol 147:195-197.
* Steinegger M, Söding J (2017) MMseqs2 enables sensitive protein sequence searching for the
  analysis of massive data sets. Nat Biotechnol 35:1026-1028.
