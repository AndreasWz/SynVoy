"""From the pipeline's output folder to what the result page shows.

A finished search leaves ``synvoy_report.json`` and figures (docs/OUTPUT.md). This
module reads them and writes a small, self-contained *bundle*:

    result.json          everything the page needs, in the site's plain wording
    files/               the figures and downloads it links to

The same bundle is built for a search run on the server and for a finished run turned
into a gallery example, so both are shown by the same page code.

Wording. The report labels each call HIGH / MEDIUM / AMBIGUOUS and gives it a class.
The page says:

    found       HIGH       counted as the gene
    probable    MEDIUM     counted as the gene
    candidate   AMBIGUOUS  right neighbourhood, sequence does not confirm it
    set_aside   class paralog_not_goi or identity_coverage_decoupled: the report
                itself leaves these out of its headline

``self_check`` in the result compares the counts made here with the report's own
headline numbers, so a new rule in the pipeline cannot silently change what the page
claims.

A genome without a call is a result as well. For every genome the bundle says how many
of the home gene's neighbours the search placed there (``neighbours``) and, where no
call was made, how many weak matches of the gene itself it saw (``leads``). Both are
counted from the gene models the search wrote (``plot_inputs_*/<genome>.gff``).
"""

import json
import os
import re
import shutil
import zipfile
from pathlib import Path

REPORT_NAME = "synvoy_report.json"
SCHEMA = 1
MAX_FILE_BYTES = 40 * 1024 * 1024

_FASTA_SUFFIX_RE = re.compile(r"\.(fna|fa|fasta)(\.gz)?$")
_LOCUS_FILE_RE = re.compile(r"^synteny_block_locus_(\d+)_(.+)$")
_SET_ASIDE = {"paralog_not_goi": "related_gene", "identity_coverage_decoupled": "short_match"}
_LABELS = {"HIGH": "found", "MEDIUM": "probable", "AMBIGUOUS": "candidate"}
_LABEL_RANK = {"found": 0, "probable": 1, "candidate": 2, "set_aside": 3}
_BLOCK_RE = re.compile(r"_b(\d+)_")            # the search block in a model's ID
ALL_FILES_PREFIX = "synvoy_"                    # files/synvoy_<name>.zip: everything in one

# Suffix after "synteny_block_locus_<N>_"  ->  key in result.json
_LOCUS_FILES = (
    ("anchor_grid.svg", "grid"),
    ("anchor_grid.html", "grid_page"),
    ("anchor_grid_with_fragments.svg", "grid_fragments"),
    ("synteny_plot_view.svg", "ribbons"),
    ("synteny_plot.html", "ribbons_page"),
    ("gene_positions.svg", "positions"),
    ("tree.html", "tree_page"),
    ("anchor_grid_print.svg", "print"),
    ("anchor_grid_print.legend.txt", "print_legend"),
    ("anchor_grid_print.source.tsv", "print_values"),
    ("gene_names.tsv", "gene_names"),
)


class ResultError(Exception):
    """The output folder does not hold a usable result."""


def genome_id(name):
    """'Apis_florea.fna' -> 'Apis_florea' (the identifier used in the report records)."""
    return _FASTA_SUFFIX_RE.sub("", str(name))


def _number(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if number != number or number in (float("inf"), float("-inf")):
        return None
    return number


def _read_json(path):
    try:
        with open(str(path)) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def _inside(path, root):
    try:
        Path(path).resolve().relative_to(Path(root).resolve())
        return True
    except (OSError, ValueError):
        return False


def _copy(source, root, files_dir, name):
    """Copy one regular file of the output folder into the bundle. Returns the name
    relative to the bundle, or None when the file is missing or not acceptable."""
    source = Path(source)
    if not source.is_file() or not _inside(source, root):
        return None
    if source.stat().st_size > MAX_FILE_BYTES:
        return None
    shutil.copyfile(str(source), str(files_dir / name))
    return "files/" + name


def classify_call(record):
    """(label, reason) for one record of ``goi_dedup.records``."""
    reason = _SET_ASIDE.get(record.get("goi_class") or "")
    if reason:
        return "set_aside", reason
    label = _LABELS.get((record.get("confidence") or "").upper())
    return (label, "") if label else (None, "")


def read_calls(report):
    """Calls of the report grouped by genome: {genome_id: [call, ...]}."""
    calls = {}
    records = (report.get("goi_dedup") or {}).get("records") or []
    for record in records:
        if not isinstance(record, dict) or not record.get("genome"):
            continue
        label, reason = classify_call(record)
        if label is None:
            continue
        coverage = _number(record.get("query_cov"))
        owner = str(record.get("owning_gene") or "")
        loci = record.get("provenance") or []
        call = {
            "label": label,
            "reason": reason,
            "confidence": str(record.get("confidence") or ""),
            "class": str(record.get("goi_class") or ""),
            "identity": _number(record.get("identity")),
            "coverage": coverage,
            "chrom": str(record.get("chrom") or ""),
            "start": int(_number(record.get("start")) or 0),
            "end": int(_number(record.get("end")) or 0),
            "model": str(record.get("model_status") or ""),
            "locus": str(record.get("owning_locus") or (loci[0] if loci else "")),
            "by_rescue": record.get("goi_class") == "synteny_hull_rescue"
                         or record.get("original_goi_class") == "synteny_hull_rescue",
            "related_gene": re.sub(r"^gene-", "", owner) if reason == "related_gene" else "",
            "low_quality": bool(record.get("low_quality_assembly")),
        }
        calls.setdefault(genome_id(record["genome"]), []).append(call)
    for genome_calls in calls.values():
        genome_calls.sort(key=lambda c: (_LABEL_RANK[c["label"]], -(c["identity"] or 0.0),
                                         c["chrom"], c["start"]))
    return calls


def genome_status(calls):
    labels = {c["label"] for c in calls}
    if labels & {"found", "probable"}:
        return "found"
    if "candidate" in labels:
        return "candidate"
    return "none"


def summarise(report, names=None, order=None):
    """Per-genome table and the counts for the top of the page.

    ``names`` maps a genome id to ``{"species", "common_name"}``; ``order`` lists the
    genome ids as they should appear (closest relatives first).
    """
    names = names or {}
    calls = read_calls(report)
    qc = {}
    for entry in report.get("genome_qc") or []:
        if isinstance(entry, dict) and entry.get("genome"):
            qc[genome_id(entry["genome"])] = entry
    qc_summary = report.get("qc_summary") or {}
    failed_qc = {genome_id(e.get("genome", "")) for e in qc_summary.get("failed_genomes") or []
                 if isinstance(e, dict)}
    kept_after_fail = {genome_id(g) for g in
                       qc_summary.get("failed_qc_genomes_with_downstream_results") or []}
    low_quality = {genome_id(g) for g in
                   (report.get("summary") or {}).get("low_quality_genomes_searched") or []}

    ids = list(order or [])
    for extra in sorted(set(qc) | set(calls) | set(names)):
        if extra not in ids:
            ids.append(extra)

    genomes = []
    for gid in ids:
        genome_calls = calls.get(gid, [])
        searched = not (gid in failed_qc and gid not in kept_after_fail
                        and gid not in low_quality and not genome_calls)
        entry = qc.get(gid) or {}
        genomes.append({
            "id": gid,
            "species": (names.get(gid) or {}).get("species") or gid.replace("_", " "),
            "common_name": (names.get(gid) or {}).get("common_name") or "",
            "status": genome_status(genome_calls) if searched else "not_searched",
            "low_quality": gid in low_quality or any(c["low_quality"] for c in genome_calls),
            "related_gene_only": bool(genome_calls) and all(
                c["label"] == "set_aside" for c in genome_calls),
            "size_mb": round((_number(entry.get("total_len")) or 0) / 1e6, 1) or None,
            "calls": genome_calls,
        })

    all_calls = [c for g in genomes for c in g["calls"]]
    counted = [c for c in all_calls if c["label"] in ("found", "probable")]
    identities = [c["identity"] for c in counted if c["identity"] is not None]
    report_summary = report.get("summary") or {}
    summary = {
        "genomes": len(genomes),
        "found": sum(1 for g in genomes if g["status"] == "found"),
        "clear": sum(1 for g in genomes if any(c["label"] == "found" for c in g["calls"])),
        "candidate": sum(1 for g in genomes if g["status"] == "candidate"),
        "none": sum(1 for g in genomes if g["status"] == "none"),
        "not_searched": sum(1 for g in genomes if g["status"] == "not_searched"),
        "copies": len(counted),
        "copies_clear": sum(1 for c in counted if c["label"] == "found"),
        "copies_probable": sum(1 for c in counted if c["label"] == "probable"),
        "candidates": sum(1 for c in all_calls if c["label"] == "candidate"),
        "set_aside": sum(1 for c in all_calls if c["label"] == "set_aside"),
        "weak_leads": int(_number(report_summary.get("low_confidence_goi")) or 0),
        "identity_min": min(identities) if identities else None,
        "identity_max": max(identities) if identities else None,
    }
    # The report's own headline counts (after its paralog check). They must agree.
    expected_clear = report_summary.get("high_confidence_goi")
    expected_probable = report_summary.get("medium_confidence_goi")
    self_check = {
        "agrees": (expected_clear in (None, summary["copies_clear"])
                   and expected_probable in (None, summary["copies_probable"])),
        "report_high": expected_clear,
        "report_medium": expected_probable,
    }
    return genomes, summary, self_check


def _attributes(text):
    out = {}
    for part in text.strip().split(";"):
        key, sep, value = part.partition("=")
        if sep:
            out[key.strip()] = value.strip()
    return out


def placed_in_genome(gff_path):
    """What the search placed in one genome for one home locus, read from the gene
    models it wrote. Neighbours (``SynVoyRole=flanking``) are counted as home genes:
    ``found`` anywhere in the genome, ``together`` in the best single place (a search
    block), ``places`` = blocks that hold any. ``leads`` are the models of the gene
    itself (``SynVoyRole=goi``), whatever their confidence."""
    genes, by_place, leads, best = set(), {}, 0, None
    try:
        with open(str(gff_path), errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] not in ("mRNA", "gene"):
                    continue
                attributes = _attributes(fields[8])
                role = attributes.get("SynVoyRole")
                if role == "flanking":
                    gene = attributes.get("SynVoy_Parent") or attributes.get("Name")
                    if not gene:
                        continue
                    genes.add(gene)
                    block = _BLOCK_RE.search(attributes.get("ID", ""))
                    if block:
                        by_place.setdefault((fields[0], block.group(1)), set()).add(gene)
                elif role == "goi":
                    leads += 1
                    identity = _number(attributes.get("Identity"))
                    if identity is not None and (best is None or identity > best[0]):
                        best = (identity, _number(attributes.get("QueryCoverage")))
    except OSError:
        return None
    return {
        "found": len(genes),
        "together": max([len(found) for found in by_place.values()] + [0]),
        "places": len(by_place),
        "leads": leads,
        "lead_identity": best[0] if best else None,
        "lead_coverage": best[1] if best else None,
    }


def _searched_neighbours(out_dir, number, fallback):
    """How many neighbouring genes of home locus ``number`` the search used as
    landmarks: the genes with a protein in its start database. (The list of the
    neighbourhood also holds genes without a protein, which cannot be searched for.)"""
    path = Path(out_dir) / "intermediate" / "initial_db" / ("initial_db_locus_%d.faa" % number)
    genes = set()
    try:
        with open(str(path), errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    name = line[1:].split()[0].split("|")[0] if line[1:].strip() else ""
                    if name and not name.startswith("GOI"):
                        genes.add(name)
    except OSError:
        return fallback
    return len(genes) or fallback


def neighbourhoods(out_dir, loci):
    """{genome id: what was placed there}, for the home locus where the most
    neighbours were found side by side. ``loci`` is the list made by ``_home_loci``."""
    out_dir = Path(out_dir)
    best = {}
    for locus in loci:
        folder = out_dir / ("plot_inputs_synteny_block_locus_%d" % locus["n"])
        if not folder.is_dir():
            continue
        total = _searched_neighbours(out_dir, locus["n"], locus.get("neighbours"))
        for path in sorted(folder.glob("*.gff")):
            stem = path.name[:-len(".gff")]
            if not _FASTA_SUFFIX_RE.search(stem) or not _inside(path, out_dir):
                continue                        # the home annotation sits in the same folder
            placed = placed_in_genome(path)
            if placed is None:
                continue
            placed.update({"locus": locus["n"], "total": total})
            gid = genome_id(stem)
            rank = (placed["together"], placed["found"], -locus["n"])
            if gid not in best or rank > best[gid][0]:
                best[gid] = (rank, placed)
    return {gid: placed for gid, (_, placed) in best.items()}


def _home_loci(out_dir):
    loci = []
    folder = Path(out_dir) / "intermediate" / "split_loci"
    for path in sorted(folder.glob("locus_*.bed")) if folder.is_dir() else []:
        match = re.match(r"^locus_(\d+)\.bed$", path.name)
        if not match:
            continue
        try:
            fields = path.read_text().splitlines()[0].split("\t")
            locus = {"n": int(match.group(1)), "chrom": fields[0],
                     "start": int(fields[1]), "end": int(fields[2])}
        except (OSError, IndexError, ValueError):
            continue
        flanking = Path(out_dir) / "intermediate" / "flanking" / (
            "synteny_block_locus_%d.bed" % locus["n"])
        try:
            rows = [l.split("\t") for l in flanking.read_text().splitlines() if l.strip()]
            locus["neighbours"] = sum(1 for r in rows if len(r) > 3 and not r[3].startswith("GOI"))
        except OSError:
            locus["neighbours"] = None
        loci.append(locus)
    return sorted(loci, key=lambda l: l["n"])


def _zip(path, members):
    """Write ``members`` ([(name in archive, file)]) with fixed dates, so the same
    input gives the same archive."""
    with zipfile.ZipFile(str(path), "w", zipfile.ZIP_DEFLATED) as archive:
        for name, source in sorted(members):
            info = zipfile.ZipInfo(name, date_time=(2000, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            with open(str(source), "rb") as fh:
                archive.writestr(info, fh.read())


def _gene_model_files(out_dir):
    """Per-genome gene models written by the search and by the rescue steps. The copy
    of the home annotation that sits in the same folder is left out."""
    members = []
    out_dir = Path(out_dir)
    for folder in sorted(out_dir.glob("plot_inputs_synteny_block_locus_*")):
        locus = folder.name.replace("plot_inputs_synteny_block_", "")
        for path in sorted(folder.glob("*.gff")):
            if _FASTA_SUFFIX_RE.search(path.name[:-len(".gff")]) and path.is_file() \
                    and path.stat().st_size <= MAX_FILE_BYTES:
                members.append(("%s/%s" % (locus, path.name), path))
    rescue = out_dir / "rescue"
    if rescue.is_dir():
        for path in sorted(rescue.glob("locus_*/*")):
            if path.suffix in (".gff", ".faa") and path.is_file() \
                    and path.stat().st_size <= MAX_FILE_BYTES:
                members.append(("%s/rescue/%s" % (path.parent.name, path.name), path))
    return [(name, path) for name, path in members if _inside(path, out_dir)]


_RESULT_WORD = {"found": "found", "probable": "probable", "candidate": "candidate",
                "set_aside": "set aside"}
TABLE_COLUMNS = ("species", "genome", "result", "confidence", "identity_pct",
                 "query_coverage_pct", "scaffold", "start", "end",
                 "neighbours_side_by_side", "neighbours_anywhere", "neighbours_searched",
                 "places_with_neighbours", "weak_matches", "note")


def genome_table(genomes):
    """The table of the result page as rows of text: one row per call, and one row
    for every genome without a call."""
    rows = []
    for genome in genomes:
        seen = genome.get("neighbours") or {}
        leads = genome.get("leads") or {}
        common = {
            "species": genome["species"], "genome": genome["id"],
            "neighbours_side_by_side": seen.get("together", ""),
            "neighbours_anywhere": seen.get("found", ""),
            "neighbours_searched": seen.get("total") or "",
            "places_with_neighbours": seen.get("places", ""),
            "weak_matches": leads.get("n", ""),
        }
        calls = genome["calls"]
        if not calls:
            rows.append(dict(common, result="not searched" if genome["status"] == "not_searched"
                             else "not found",
                             note="assembly too fragmented" if genome["status"] == "not_searched"
                             else ""))
        for call in calls:
            note = {"related_gene": "matches the related gene %s better" % call["related_gene"]
                    if call["related_gene"] else "matches a related gene better",
                    "short_match": "matches only a short piece of the query"}.get(call["reason"], "")
            rows.append(dict(
                common, result=_RESULT_WORD[call["label"]], confidence=call["confidence"],
                identity_pct="" if call["identity"] is None else "%.1f" % call["identity"],
                query_coverage_pct="" if call["coverage"] is None else "%.0f" % (call["coverage"] * 100),
                scaffold=call["chrom"], start=call["start"] or "", end=call["end"] or "", note=note))
    lines = ["\t".join(TABLE_COLUMNS)]
    for row in rows:
        lines.append("\t".join(str(row.get(column, "")) for column in TABLE_COLUMNS))
    return "\n".join(lines) + "\n"


def short_name(query):
    """What a search goes by: the gene symbol when there is one, else the name."""
    return str(query.get("gene") or query.get("name") or query.get("label") or "")


def _read_me(result, names):
    query, home, summary = result["query"], result["home"], result["summary"]
    pipeline, timing = result["pipeline"], result["timing"]
    title = "%s in %d genome%s" % (short_name(query) or "The gene",
                                   summary["genomes"], "" if summary["genomes"] == 1 else "s")
    software = "SynVoy " + (pipeline.get("label") or "pre-release")
    if pipeline.get("commit"):
        software += ", commit %s%s" % (pipeline["commit"],
                                       " with local changes" if pipeline.get("modified") else "")
    facts = [
        ("Query", ", ".join(str(x) for x in (
            query.get("gene"), query.get("name"), query.get("accession"),
            "%s residues" % query["length"] if query.get("length") else "",
            query.get("organism")) if x)),
        ("Home genome", ", ".join(x for x in (home.get("species"), home.get("assembly")) if x)),
        ("Genomes", "%d searched: the gene found in %d, a candidate in %d, no call in %d%s" % (
            summary["genomes"] - summary["not_searched"], summary["found"], summary["candidate"],
            summary["none"], ", %d not searched" % summary["not_searched"]
            if summary["not_searched"] else "")),
        ("Software", software),
        ("Settings", ", ".join(x for x in (
            "profile " + pipeline["profile"] if pipeline.get("profile") else "",
            pipeline.get("preset") or "") if x)),
        ("Finished", (timing.get("finished_at") or "").replace("T", " ").replace("Z", " UTC")),
    ]
    files = [
        ("genome_by_genome.tsv", "one row per call and one per genome without a call: result, "
                                 "identity, coverage, place, and how many neighbouring genes "
                                 "were found there"),
        ("synvoy_report.json", "every call of the search with its evidence, as the pipeline "
                               "wrote it (docs/OUTPUT.md in the SynVoy repository)"),
        ("query.faa", "the query protein"),
        ("figures/", "the figures of the result page; *_anchor_grid_print.svg is the grid at "
                     "print size, with a draft of its legend (.legend.txt) and the value of "
                     "every cell (.source.tsv)"),
        ("gene_models/", "per home locus, one GFF per genome with the gene models the search "
                         "wrote (the gene and its neighbours); rescue/ holds the models of "
                         "the rescue steps"),
    ]
    lines = ["SynVoy result: " + title, "=" * (15 + len(title)), ""]
    lines += ["%-12s %s" % (key, value) for key, value in facts if value]
    lines += ["", "Files", "-----"]
    lines += ["%-22s %s" % (name, text) for name, text in files
              if any(member == name or member.startswith(name) for member in names)]
    lines += ["", "A genome without a call was searched: its row in genome_by_genome.tsv says how",
              "much of the neighbourhood was found. A missing gene model is not evidence that",
              "the gene is absent.", ""]
    for note in result.get("notes") or []:
        lines += ["Note: " + note["text"], ""]
    return "\n".join(lines)


def build_bundle(out_dir, bundle_dir, meta, figures_dir=None):
    """Write the bundle for one finished run and return the content of result.json.

    ``meta`` carries what the output folder does not know: ``id``, ``query``, ``home``,
    ``targets`` (list of ``{"stem", "species", "common_name"}``), ``pipeline``,
    ``timing``, ``scope``, ``failed_tasks`` and ``notes`` (warnings to show on the page). ``figures_dir`` replaces the folder the
    figures are taken from (re-drawn figures for a gallery example).

    The bundle is built next to its final place and moved there in one step, so a
    half-written result is never visible.
    """
    out_dir, bundle_dir = Path(out_dir), Path(bundle_dir)
    report = _read_json(out_dir / REPORT_NAME)
    if not isinstance(report, dict) or not isinstance(report.get("goi_dedup"), dict):
        raise ResultError("%s is missing or not a SynVoy report" % REPORT_NAME)

    targets = meta.get("targets") or []
    names = {t["stem"]: {"species": t.get("species"), "common_name": t.get("common_name")}
             for t in targets if t.get("stem")}
    genomes, summary, self_check = summarise(report, names, [t["stem"] for t in targets
                                                             if t.get("stem")])
    summary["failed_tasks"] = int(meta.get("failed_tasks") or 0)

    building = bundle_dir.with_name(bundle_dir.name + ".building")
    if building.exists():
        shutil.rmtree(str(building))
    files_dir = building / "files"
    files_dir.mkdir(parents=True)

    figure_root = Path(figures_dir) if figures_dir else out_dir
    loci = {}
    for path in sorted(figure_root.iterdir()) if figure_root.is_dir() else []:
        match = _LOCUS_FILE_RE.match(path.name)
        if not match:
            continue
        number, rest = int(match.group(1)), match.group(2)
        for suffix, key in _LOCUS_FILES:
            if rest == suffix:
                copied = _copy(path, figure_root, files_dir, path.name)
                if copied:
                    loci.setdefault(number, {})[key] = copied
    for number in list(loci):
        newick = out_dir / ("locus_%d_tree.nwk" % number)
        try:
            # With the tree step switched off the pipeline writes a one-leaf stand-in.
            real_tree = newick.is_file() and "GOI_placeholder" not in newick.read_text()[:2000]
        except OSError:
            real_tree = False
        if real_tree:
            copied = _copy(newick, out_dir, files_dir, newick.name)
            if copied:
                loci[number]["tree_newick"] = copied

    home_loci = {l["n"]: l for l in _home_loci(out_dir)}
    locus_list = []
    for number in sorted(set(loci) | set(home_loci)):
        entry = {"n": number, "files": loci.get(number, {})}
        entry.update({k: v for k, v in home_loci.get(number, {}).items() if k != "n"})
        locus_list.append(entry)
    summary["loci"] = len(locus_list)

    # What the search placed in each genome, whether or not it made a call there.
    placed = neighbourhoods(out_dir, [{"n": l["n"], "neighbours": l.get("neighbours")}
                                      for l in locus_list])
    for genome in genomes:
        seen = placed.get(genome["id"])
        genome["neighbours"] = None if seen is None else {
            key: seen[key] for key in ("found", "together", "places", "total", "locus")}
        genome["leads"] = None
        if seen is not None and genome["status"] == "none" and seen["leads"]:
            genome["leads"] = {"n": seen["leads"], "identity": seen["lead_identity"],
                               "coverage": seen["lead_coverage"]}

    downloads = []

    def offer(key, relative):
        if relative:
            size = (building / relative).stat().st_size
            downloads.append({"key": key, "file": relative, "bytes": size})

    offer("report", _copy(out_dir / REPORT_NAME, out_dir, files_dir, REPORT_NAME))
    members = _gene_model_files(out_dir)
    if members:
        _zip(files_dir / "gene_models.zip", members)
        offer("gene_models", "files/gene_models.zip")
    fasta = (meta.get("query") or {}).get("fasta")
    if fasta:
        (files_dir / "query.faa").write_text(fasta)
        offer("query", "files/query.faa")
    else:
        offer("query", _copy(out_dir / "intermediate" / "query" / "normalized_query.faa",
                             out_dir, files_dir, "query.faa"))

    params = _read_json(out_dir / "intermediate" / "locate_gene" / "effective_params.json") or {}
    # The report's ``gene_family_advisory`` is not shown: its locus count is the number
    # of (genome, locus, chromosome) cells, so it fires on single-copy genes too.
    notes = []
    for text in report.get("warnings") or []:
        if isinstance(text, str) and text.strip():
            notes.append({"kind": "report", "text": text.strip()[:600]})
    for note in meta.get("notes") or []:
        notes.append({"kind": str(note.get("kind") or "note"), "text": str(note.get("text") or "")[:600]})

    (files_dir / "genome_by_genome.tsv").write_text(genome_table(genomes))
    offer("table", "files/genome_by_genome.tsv")

    query = dict(meta.get("query") or {})
    query.pop("fasta", None)
    result = {
        "schema": SCHEMA,
        "id": meta.get("id") or "",
        "scope": meta.get("scope") or "",
        "query": query,
        "home": meta.get("home") or {},
        "pipeline": dict(meta.get("pipeline") or {}, preset=str(params.get("preset_applied") or "")),
        "timing": meta.get("timing") or {},
        "summary": summary,
        "self_check": self_check,
        "notes": notes,
        "genomes": genomes,
        "loci": locus_list,
        "downloads": downloads,
    }

    # Everything in one archive, so nobody has to click through the single files.
    archive = [("synvoy_report.json", files_dir / REPORT_NAME),
               ("query.faa", files_dir / "query.faa"),
               ("genome_by_genome.tsv", files_dir / "genome_by_genome.tsv")]
    for locus in locus_list:
        for relative in sorted(set(locus["files"].values())):
            archive.append(("figures/" + relative.split("/", 1)[1], building / relative))
    archive += [("gene_models/" + name, path) for name, path in members]
    archive = [(name, path) for name, path in archive if Path(path).is_file()]
    read_me = building / "README.txt"
    read_me.write_text(_read_me(result, [name for name, _ in archive]))
    label = re.sub(r"[^A-Za-z0-9-]+", "_", short_name(query))
    zip_name = ALL_FILES_PREFIX + (label.strip("_-")[:40].rstrip("_-") or "result") + ".zip"
    _zip(files_dir / zip_name, archive + [("README.txt", read_me)])
    read_me.unlink()
    downloads.insert(0, {"key": "all", "file": "files/" + zip_name,
                         "bytes": (files_dir / zip_name).stat().st_size})

    temporary = building / "result.json.tmp"
    with open(str(temporary), "w") as fh:
        json.dump(result, fh, indent=1, sort_keys=True)
        fh.write("\n")
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(str(temporary), str(building / "result.json"))
    if bundle_dir.exists():
        shutil.rmtree(str(bundle_dir))
    os.replace(str(building), str(bundle_dir))
    return result


def bundle_files(result):
    """Every file a result page may ask for: the only names the API will serve."""
    names = set()
    for locus in result.get("loci") or []:
        names.update((locus.get("files") or {}).values())
    for item in result.get("downloads") or []:
        names.add(item.get("file"))
    return {n for n in names if isinstance(n, str) and n.startswith("files/")}


def card(result):
    """The few numbers a gallery card shows."""
    summary = result.get("summary") or {}
    query = result.get("query") or {}
    return {
        "name": query.get("name") or query.get("label") or "",
        "length": query.get("length"),
        "scope": result.get("scope") or "",
        "runtime_s": (result.get("timing") or {}).get("runtime_s"),
        "genomes": summary.get("genomes"),
        "found": summary.get("found"),
        "clear": summary.get("clear"),
        "candidate": summary.get("candidate"),
        "none": summary.get("none"),
        "identity_min": summary.get("identity_min"),
        "identity_max": summary.get("identity_max"),
    }
