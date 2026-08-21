# TP53 / TP63 / TP73 paralog-discrimination assessment

Assessment of `docs/TODO.md §2c` — *does SynVoy correctly separate the TP53 paralog
family into its three syntenic loci?* Source run: `results/tp53_15042026/` (query =
human TP53, easy mode, 6 vertebrate target genomes). Analysis is static over the
captured per-locus homology TSVs — no rerun. Date: 2026-05-27.

> **Status update (2026-07-26): the two fixes this document asks for have landed;
> the re-score it recommends has not been done.**
>
> - **§1j reciprocal-best check — implemented.** `bin/reciprocal_best_paralog_check.py`
>   SW-aligns each recovered target against every home paralog;
>   `generate_report.py` flags a call whose best-matching paralog differs from the
>   per-cell modal one as `paralog_confusion` in the `self_consistency` block.
>   Single-paralog queries are an automatic no-op. Toggle: `--disable_paralog_check`.
> - **Locus ownership (§1m) — implemented, and goes further than §1j.** The home-paralog
>   panel is built from the whole home proteome (`bin/build_home_paralog_panel.py`), so a
>   call whose best home match is a *family paralog* rather than the GOI is relabelled
>   `paralog_not_goi` and **excluded from the headline HIGH/MEDIUM counts** — it is
>   demoted, not merely flagged. Toggle: `--disable_locus_ownership`.
> - **§2a eggNOG OG-gating — not implemented.** No eggNOG code exists in the repo.
> - **The recommended re-run was never scored.** The TP53 target genomes were removed
>   from disk, so reproducing §2's MEDIUM-bleed numbers needs a refetch first. Until
>   then, treat every number below as describing the **pre-§1j/§1m** behaviour: the
>   MEDIUM family bleed it measures is precisely what those two fixes target, so the
>   figures here are an upper bound on the problem, not the current state.

## 1. Home-locus separation — clean ✅

`split_loci.py` split the TP53 query into **three home loci that land on the correct
human paralog positions**:

| Locus | Home coords | Paralog | Cytoband (expected) |
|-------|-------------|---------|---------------------|
| locus_1 | NC_000017.11:7.67 Mb | **TP53** | 17p13.1 ✓ |
| locus_2 | NC_000001.11:3.72 Mb | **TP73** | 1p36.3 ✓ |
| locus_3 | NC_000003.12:189.8 Mb | **TP63** | 3q28 ✓ |

All three survive the §1d locus cap (3 ≤ `max_loci=5`). This is the prerequisite for
paralog discrimination and it works.

## 2. Target ortholog assignment — strong at HIGH, family bleed at MEDIUM 🟡

Counting every HIGH/MEDIUM GOI call across the 6 target genomes by the target genome's
own annotated product (RefSeq products = ground truth), grouped by locus:

| Locus (home paralog) | HIGH-confidence calls | MEDIUM-confidence calls |
|----------------------|-----------------------|--------------------------|
| locus_1 (TP53) | **TP53 ×1** (pure) | TP53 ×1 (pure) |
| locus_2 (TP73) | **TP73 ×6**, TP63 ×1 | TP53 ×1, TP63 ×2, TP73 ×1 |
| locus_3 (TP63) | **TP63 ×1** (pure) | TP73 ×6, TP63 ×7, TP53 ×1 |

**At HIGH confidence the locus→paralog assignment is essentially correct** (13/14 HIGH
calls name the locus's own paralog; the lone exception is one TP63 product at the TP73
locus). Per-genome best calls confirm this: e.g. mouse → Trp53@locus_1 (79%),
Trp73@locus_2 (88%), p63@locus_3; `GCF_049306965.1` → p53/p73/p63 each at the right locus.

**At MEDIUM confidence there is clear family cross-talk** — the TP63 locus in particular
surfaces as many TP73-annotated genes (6) as TP63 ones (7). The iterative search casts a
wide net at lower confidence and pulls in sibling paralogs from the broader p53 family.

Sensitivity is paralog-dependent and tracks conservation: **TP73 (most conserved) is
recovered HIGH in 5/6 genomes, TP63 in ~4/6, TP53 (fastest-evolving, shortest) is the
weakest** — only 2/6 genomes yield a real p53 call; elsewhere a syntenic neighbour
(SRY-box, growth-hormone-regulator) is the MEDIUM best-in-block, a textbook §1e case.

## 3. Verdict & recommendation

- **Publishable as a second benchmark *if scoped to HIGH-confidence calls*** — there the
  story is clean and is exactly the discrimination standard homology tools fail: each
  recovered ortholog is assigned to the correct paralogous locus, not lumped into one
  "TP53-family" hit.
- **Not yet clean across HIGH+MEDIUM** — the MEDIUM family bleed must be characterised or
  suppressed before claiming SynVoy "uniquely resolves" the paralog problem at all
  confidence levels. The two fixes that target this directly:
  - **§1j reciprocal-best check** — flag a call as `paralog_confusion` when the recovered
    target is closer to a *different* home paralog than to the locus's own. This would
    catch the locus_3 TP73 bleed mechanically.
  - **§2a eggNOG OG-gating** — demote calls whose orthologous group doesn't match the
    locus paralog's OG.
- **Next step:** re-run with `-profile preset_paralog_discrimination` (now auto-suggested
  by §1f) + the landed §1e, then land §1j and re-score. That run becomes the paper figure.
