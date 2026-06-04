# SynVoy — Realistic Progress & Soundness Assessment

**Date:** 2026-06-03
**Author of assessment:** working notes (Claude), reviewed against this session's DCN findings
**Status of project:** early development, in active use for the Ant_Venoms paper

> This is written to be honest, not encouraging. The "it's breaking apart" feeling is
> a real signal and it points at something specific and fixable — not at a doomed tool.

---

## 0. TL;DR

- **The core idea works, and is validated, for the use case it was built for:** finding
  *divergent, low-copy* genes whose flanking synteny is conserved (melittin: 11/12 vs
  2/12 for BLAST/MMseqs2 baselines). That result is real and publishable.
- **What you're feeling is scope creep, not collapse.** Every new "it doesn't work" case
  (firefly luciferase, decorin/SLRPs) is a **large, promiscuous gene family** — precisely
  the regime the architecture's three load-bearing assumptions break in. The patches
  (§1d–§1m, 26 numbered fixes, a 5,400-line search core) are whack-a-mole *because* the
  tool is being pushed past its validated domain.
- **The single most important decision is a scoping decision, not a coding one:** declare
  what SynVoy is *for*, validate that, and either (a) handle gene families as an explicit,
  first-class, honestly-limited mode, or (b) detect them and tell the user "manual curation
  required." Both are fine. Pretending the same pipeline transparently handles melittin and
  the SLRP superfamily is what's generating the edge-case avalanche.
- **There is exactly one architectural root bug worth fixing properly right now** (the
  "locus fumble" — GOI annotation is gated on flanking-seeded blocks). It hurts even
  single-copy genes when synteny is rearranged. Everything else is either in-scope-and-fine
  or out-of-scope-and-should-be-documented.

---

## 1. What SynVoy genuinely does well (don't lose sight of this)

| Strength | Evidence |
|---|---|
| **Synteny-guided recovery of divergent genes** — the core thesis | Melittin: 11/12 recovered vs 2/12 for BLAST/MMseqs2 (memory: `project_melittin_benchmark`, validated 2026-05-13). This is the headline and it holds. |
| **End-to-end operational robustness** | Easy/Pro modes, conda/docker/singularity/HPC profiles, memory tiers (laptop_safe/low_mem), genome QC gates, `-resume`, atomic JSON, 32 pytest files. This is solid engineering. |
| **Honest reporting plumbing** | §1h dedup, §1k headline metric fix, §1j paralog check, self-consistency block. The report tries hard to not lie about what was found. |
| **Shared-utility discipline** | `sequence_utils.py` as a single source of truth; helper consolidation (§13). The codebase is more maintainable than its size suggests. |

**This is a real tool with a real, measured advantage in its design domain.** For the
Ant_Venoms paper — divergent venom peptides located by conserved flanking genes — it does
the job baselines can't. Keep that framing central.

---

## 2. The architecture's three load-bearing assumptions

SynVoy's pipeline is: *locate the query's home loci → take flanking genes → find the
syntenic neighborhood in target genomes → search locally for the gene.* This rests on
three assumptions:

1. **A1 — The query maps to a small number of clean home loci.**
   (split_loci ranks the query's home hits and keeps the top `max_loci`.)
2. **A2 — Flanking-gene synteny is conserved enough to anchor a block over the gene in the
   target.** (Blocks are seeded *only* by flanking hits; the GOI is then searched *inside*
   those blocks.)
3. **A3 — A sequence found in the right neighborhood, above an identity threshold, is the
   ortholog.** (Confidence = % identity tiers; HIGH ⇒ "it's the gene".)

**For divergent, low-copy genes with conserved synteny, all three hold** → melittin/LY6 work.

**For large, promiscuous gene families, all three fail** — and they fail *together*, which is
why each family looks like a brand-new disaster.

---

## 3. Worked example: decorin (DCN) — all three assumptions break at once

Decorin is the prototype Small Leucine-Rich Proteoglycan (SLRP). Its biology is exactly
what breaks the architecture:
- It sits in a **cluster of paralogous SLRPs** (LUM/KERA/EPYC are its immediate neighbors
  in every species we checked) — so *its own flanking genes are paralogs of itself*.
- The wider family (OMD/OGN/ASPN, FMOD/PRELP, plus unrelated **LRR-domain** proteins like
  FLRT1/3, LRRC4) cross-hits decorin everywhere via the shared leucine-rich repeat.

What happened (this session, NCBI-confirmed, cow/mouse/zebrafish):

- **A1 broke.** The query hit ~23 chromosomes. `--max_loci 5` kept the 5 *highest-scoring*
  loci — DCN(chr12), PRELP(chr1), FLRT3(chr20), FLRT1(chr11), LRRC4(chr7). The
  biologically-central **OMD/OGN/ASPN cluster ranked 6th and was dropped.** Ranking is by
  alignment score, which favors tight LRR cross-hits over the true paralog cluster. *The cap
  is a runtime-control knob making a biological decision it isn't qualified to make.*
- **A2 broke.** Decorin's chr12 neighbors are rearranged across cow chromosomes, so the DCN
  home locus never formed a flanking-seeded block over cow decorin. It produced **only LOW
  `fallback_hit_span` junk at the wrong coordinates** (best on chr5 was 32.9% at 18.3 Mb;
  real decorin is at 21.01 Mb). The locus that *should* own the gene cannot find it. (This
  is the "fumble.")
- **A3 broke.** Asporin (a paralog) was recovered at ~50% inside a neighborhood and **labeled
  `GOI_DCN` / probable_goi** — i.e., a paralog reported as the gene of interest. Identity
  alone does not separate ortholog from paralog.

**Net result: the recent changes made decorin *worse than the student's original run.*** The
student's run found cow decorin at 90% — but only *accidentally*, because (without a max_loci
cap) the OMD locus's flanking cross-matched decorin's neighborhood and recovered it. §1d's
cap removed that accident; the fumble (A2) was never fixed; so decorin is now lost entirely.
§1m (locus ownership) is a correct idea but cannot re-attribute a gene that was never found,
and its panel is incomplete precisely because A1 dropped the relevant loci.

**This is the whole "breaking apart" feeling in one example: three coupled assumptions, each
patched independently, with the patches interacting badly.**

---

## 4. Biological soundness, component by component

Honest grading. "Sound" = principled and behaves correctly; "Fragile" = right idea, breaks
outside a narrow regime; "Questionable" = the biology doesn't really support it.

| Component | Grade | Notes |
|---|---|---|
| Flanking-synteny anchoring (core) | **Sound→Fragile** | Excellent for conserved synteny; degrades with rearrangement & distance. The premise (flanking order is conserved) is biologically real but has a hard horizon. |
| split_loci ranking + `max_loci` (§1d) | **Questionable** | Ranks home loci by alignment score and hard-caps. For families this drops the biologically-relevant cluster. Score-rank ≠ biological relevance. |
| GOI block search / fallback (A2) | **Fragile (root bug)** | GOI annotation gated on flanking blocks ⇒ the true locus fumbles when its synteny is rearranged. Architectural, not family-specific — hurts single-copy cases too. |
| Identity-tier confidence (A3) | **Questionable** | % identity does not separate ortholog from paralog. A 50% hit can be a close ortholog or a diverged paralog. This is the source of the asporin mislabel. |
| §1g distance auto-tune | **Questionable (near-no-op)** | The proxy (median flanking identity) is *decoupled* from GOI divergence — flanking stays 73–95% even when the GOI is wildly diverged. Fires essentially never on real data (memory: `project_1g_distance_autotune`). Currently theater. |
| §1h cross-locus dedup | **Sound** | Collapsing same-coordinate multi-locus hits is correct and well-tested. |
| §1j reciprocal-best paralog check | **Sound (the right direction)** | Reciprocal-best alignment is *the* standard local ortholog signal. Underused. |
| §1m locus ownership (RBH + synteny) | **Sound idea, dependency-limited** | Biologically the most principled recent addition (RBH against a home-paralog panel). But it's a *report-level* patch whose correctness depends on A1 having kept the right loci. Circular. |
| §1e strong-synteny rescue | **Reasonable** | Relaxed miniprot in strong-flanking/no-GOI blocks — a sensible safety net. |
| Confidence/QC/dedup reporting | **Sound** | Operationally honest. |

**Cross-cutting scientific gap:** SynVoy conflates *"a similar sequence in the syntenic
neighborhood"* with *"the ortholog."* For single-copy genes that's fine. For families,
ortholog/paralog assignment is a genuine, hard problem (it's the entire field of orthology
inference — RBH+synteny, or gene-tree/species-tree reconciliation). SynVoy gestures at the
first (§1j/§1m) but does not do it rigorously or as a first-class step.

---

## 5. The meta-pattern: scope creep & patch accumulation

- **26** numbered fix-sections in CLAUDE.md; a **5,443-line** search core; **27,899** lines
  in `bin/`. The numbered-fix density is the tell: most are reactive responses to a *new
  use case* exposing an old assumption.
- The failing cases share a signature: **luciferase (AMP-binding family, 49 home loci)** and
  **decorin (SLRP superfamily)** are both **large promiscuous families**. The working cases
  (**melittin**, **LY6**) are **divergent but low-copy**.
- Each family triggers a new patch (§1d for the luciferase 49-loci blowup, §1m for the SLRP
  mis-attribution). The patches are individually reasonable and individually tested — but
  they're patching *symptoms of the same out-of-domain use*, and they interact (the DCN
  regression is a §1d×fumble interaction).

This is the classic trajectory of a research tool that works, gets pushed past its design
envelope, and accretes edge-case handling faster than it accretes validated capability. It
is **fixable by scoping, not by more patches.**

---

## 6. Recommendations

### 6.1 Decide and document the validated scope (highest leverage, ~0 code)
Write one paragraph in the README and the paper methods: **SynVoy locates divergent,
low-copy genes via conserved flanking synteny.** State the regime it is validated in
(melittin/LY6) and the regime it is *not* (large paralog families / repeat-domain
superfamilies, where ortholog/paralog assignment is unresolved). This single act stops the
whack-a-mole: a failing SLRP run becomes "expected, out of scope," not "the tool is broken."

### 6.2 Paper strategy: lead with the strength, don't over-claim generality
The paper is *"How Ant Genomes Repeatedly Reinvent Venom."* The tool's job there is finding
divergent venom genes by synteny — its validated wheelhouse. Present it as a **specialist
method for synteny-recoverable divergent genes**, benchmarked on melittin, *not* as a general
ortholog finder. Reviewers will respect a sharp, honest claim far more than a broad fragile
one. Do **not** put decorin/SLRPs in the paper as a success case.

### 6.3 Gene families: pick (a) or (b), explicitly — do not keep patching
- **(b) Detect-and-warn (recommended, near-term).** Early in the run, detect a gene-family
  query (many home loci, high paralog cross-hit density, repeat-domain Pfam) and emit a loud
  banner: *"Gene-family query detected; ortholog/paralog assignment is not reliable; treat
  GOI calls as candidates requiring manual curation."* Make the report's `goi_class` honest
  (`ambiguous_family_member` by default for these). This is achievable in days and is the
  honest behavior.
- **(a) First-class paralog mode (longer-term, only if needed for the science).** Promote
  §1j/§1m from report patches to a real assignment step: derive the *complete* home-paralog
  panel (don't let `max_loci` truncate it), do RBH+synteny ownership *before* classification,
  and assign each target to its best paralog (ortholog candidate) with an explicit
  paralog label. This is real orthology inference and should be treated as a feature with its
  own validation set, not an edge-case fix.

### 6.4 Fix the one true architectural root bug: the "fumble" (A2)
Make a home locus able to model *its own* GOI even when its flanking is rearranged in the
target. Concretely: when GOI hits exist on a chromosome that has no flanking-seeded block,
open a GOI-seeded search window there (generalize the existing `goi_proxy` path to the real
GOI query), or run the GOI annotation on cluster_grs-style windows. This helps **single-copy
divergent genes too** (the design domain), so it's worth doing properly. It is the highest-
value *code* fix on the board.

### 6.5 Re-ground or remove the questionable knobs
- **§1g distance auto-tune:** either re-derive the threshold *relative to the observed
  flanking median* (GOI ≈ flanking − Δ) so it actually fires, or remove it and stop claiming
  distance adaptation the tool doesn't really do.
- **`max_loci` ranking:** rank by *biological* relevance (e.g., include the full reciprocal-
  best paralog set, or cluster-aware selection), not raw alignment score — or at minimum,
  never let the cap silently drop a locus that the panel/ownership step needs.
- **Identity-tier confidence:** stop treating high identity as "it's the ortholog." Gate
  HIGH on RBH-consistency (best home match is the GOI itself), not identity alone.

### 6.6 Build a real regression panel across the difficulty axes
Right now each hard case is discovered ad hoc. Define a fixed benchmark set spanning:
single-copy divergent (melittin), micro-exon, small family (LY6), large family (SLRP/DCN,
as an *expected-hard/limitations* case). Run all on every change. This converts "new use
case breaks everything" into "the SLRP regression behaves as documented."

### 6.7 Stop-the-bleeding rule
Adopt a one-line policy: **no new numbered §1x patch without first classifying the failure
as in-scope (fix root cause) or out-of-scope (document + warn).** The 26-and-counting fix
list is the symptom to treat.

---

## 7. Concrete near-term plan (2–3 focused work items, in order)

1. **Scope statement** (README + paper methods + a `goi_class` honesty pass). ~half a day.
   Immediately stops the existential "everything breaks" feeling.
2. **Gene-family detect-and-warn (6.3b).** A few days. Turns SLRP/luciferase from "bugs" into
   "honestly-labeled candidates."
3. **Fix the fumble (6.4).** The one architectural fix that helps the *core* domain. Validate
   on a single-copy rearranged case + LY6 + melittin (must not regress).

Everything else (full paralog mode, §1g/max_loci re-grounding) is deferrable and should only
be picked up if the *science* (not the desire for generality) demands it.

---

## 8. Bottom line

SynVoy is **not breaking apart — it is succeeding at its real job and being asked to do a
different, harder one.** The melittin result is genuine and paper-worthy. The decorin pain is
the architecture honestly reporting that large paralog families are outside its validated
envelope. Draw that envelope explicitly, fix the one root bug inside it (the fumble), label
the outside honestly, and the whack-a-mole stops. That is a far healthier position than it
feels like right now.
