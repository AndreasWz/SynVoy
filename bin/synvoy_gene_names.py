"""Display names for home genes in SynVoy figures.

A gene ID such as ``gene-LOC726817`` or ``gene-Dmel_CG4491`` says nothing to a
reader. This module turns each home gene into a short label plus the full name
behind it, and every figure draws its labels and its "Gene names" table from
the same table.

Resolution order, per gene:
  1. names table supplied by the user (``--gene_names_tsv``)   -> source ``user``
     (a table written by an earlier run keeps its own source values).
  2. gene symbol in the home GFF (``gene=`` / ``Name=`` ...)    -> ``gff_symbol``
  3. current symbol in NCBI Gene, looked up by the ``GeneID``
     in the GFF's ``Dbxref``; only for genes still unnamed      -> ``ncbi_gene``
  4. abbreviation of the product name (see ``abbreviate_product``)
                                                                -> ``product_name``
  5. the gene ID itself, when nothing names the gene            -> ``id_only``

Labels from step 4 are not official gene symbols. Figures mark them with
``DERIVED_MARK`` and say so under the names table.

Import-only helper (like ``synvoy_taxa``); the plotting scripts own the CLI.
"""

import json
import os
import re
from html import escape as _esc
import urllib.error
import urllib.request
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Tuple
from urllib.parse import unquote

SOURCE_USER = "user"
SOURCE_GFF = "gff_symbol"
SOURCE_NCBI = "ncbi_gene"
SOURCE_PRODUCT = "product_name"
SOURCE_ID = "id_only"
SOURCES = (SOURCE_USER, SOURCE_GFF, SOURCE_NCBI, SOURCE_PRODUCT, SOURCE_ID)

SOURCE_TEXT = {
    SOURCE_USER: "names table",
    SOURCE_GFF: "home GFF",
    SOURCE_NCBI: "NCBI Gene",
    SOURCE_PRODUCT: "product name",
    SOURCE_ID: "no name annotated",
}

# Appended to a label that SynVoy derived from the product name.
DERIVED_MARK = "*"
DERIVED_NOTE = (f"{DERIVED_MARK} = label derived by SynVoy from the annotated product "
                "name; not an official gene symbol")

# A label longer than this is cut: a rotated column header has to stay inside
# its band. Twelve characters is the length of an NCBI ``LOC`` id.
MAX_LABEL_CHARS = 12
# A product name that is one word is more readable whole than cut down to three
# letters ("junctophilin-1", not "JUN1"), so it may be a little longer.
MAX_WHOLE_NAME_CHARS = 16
# An acronym of a one- or two-word name shorter than this says too little ("AK").
MIN_ACRONYM_CHARS = 4

NCBI_GENE_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/"
NCBI_TIMEOUT_S = 20
NCBI_CHUNK = 100

TSV_COLUMNS = ("gene_id", "label", "full_name", "source", "ncbi_gene_id")

_GFF_SYMBOL_KEYS = ("gene", "gene_name", "gene_symbol", "symbol", "Name")
_GFF_PRODUCT_KEYS = ("description", "product", "Note")
_TRANSCRIPT_TYPES = {"mRNA", "transcript", "CDS", "lnc_RNA", "ncRNA", "tRNA", "rRNA",
                     "snRNA", "snoRNA", "miRNA", "misc_RNA", "antisense_RNA",
                     "guide_RNA", "primary_transcript"}

_GENERIC_ID_RES = (
    re.compile(r"^LOC\d+$", re.IGNORECASE),
    re.compile(r"^[A-Za-z]{1,8}\d*_\d+$"),                         # XY_0012345
    re.compile(r"^[A-Za-z][A-Za-z0-9]{1,11}_[A-Za-z]{0,4}\d{3,}$"),  # locus tag: Dmel_CG4491
    re.compile(r"^(ENS[A-Z]*[GTP]|FBgn|FBtr|WBGene)\d+(\.\d+)?$"),  # database accessions
    re.compile(r"^[A-Z]{1,2}_\d+(\.\d+)?$"),                       # XM_006567126.3, NP_...
    re.compile(r"^(gene|rna|cds|g|jg|evm|mRNA)[-_.:]?\d+(\.t?\d+)?$"),            # g123.t1
)

_UNINFORMATIVE_RE = re.compile(
    r"uncharacteri[sz]ed|hypothetical|unknown|unnamed|predicted protein|"
    r"domain of unknown function|^protein$|^gene$|\bLOC\d+\b", re.IGNORECASE)

# Words that carry no information in a product name.
_FILLER = {"protein", "proteins", "of", "and", "the", "a", "an", "in", "to", "for",
           "with", "type", "family", "member", "homolog", "homologue", "isoform",
           "kda", "putative", "probable", "predicted"}
# Trailing capitals that are chemistry, not a gene symbol ("V-type proton ATPase").
_BIOCHEM_ABBREVIATIONS = {"DNA", "RNA", "ATP", "ADP", "AMP", "GTP", "GDP", "NAD", "NADH",
                          "NADP", "NADPH", "FAD", "CoA", "ATPase", "GTPase"}
_ROMAN_RE = re.compile(r"^(?=[IVX])X{0,2}(IX|IV|V?I{0,3})$")
# Compartment qualifiers NCBI appends after a comma ("..., mitochondrial").
_LOCATION_RE = re.compile(
    r",\s*(mitochondrial|cytoplasmic|chloroplastic|peroxisomal|nuclear|"
    r"axonemal|ciliary|muscle|cytosolic|lysosomal|vacuolar)(\s+\S+)?$", re.IGNORECASE)


def is_generic_id(label: str) -> bool:
    """True for an identifier that is not a gene name (LOC id, locus tag, accession)."""
    txt = (label or "").strip()
    for prefix in ("gene-", "gene:", "rna-", "transcript:"):
        if txt.startswith(prefix):
            txt = txt[len(prefix):]
    if not txt:
        return True
    return any(rx.match(txt) for rx in _GENERIC_ID_RES)


def is_uninformative_product(product: str) -> bool:
    """True when a product description does not name the gene."""
    txt = (product or "").strip()
    return not txt or bool(_UNINFORMATIVE_RE.search(txt))


def clean_product(product: str) -> str:
    """Strip transcript, isoform, evidence and compartment qualifiers from a product name."""
    txt = unquote(str(product or "")).strip()
    txt = re.sub(r"\s*\[Source:[^\]]*\]", "", txt)                 # Ensembl description tail
    txt = re.sub(r"^LOW QUALITY PROTEIN:\s*", "", txt, flags=re.IGNORECASE)
    txt = re.sub(r",?\s*transcript variant\s+\S+$", "", txt, flags=re.IGNORECASE)
    txt = re.sub(r",?\s*isoform\s+\S+$", "", txt, flags=re.IGNORECASE)
    txt = re.sub(r"\s*[\[(][^\])]*[\])]", "", txt)                 # [ubiquinone], (strain X)
    txt = _LOCATION_RE.sub("", txt)
    txt = re.sub(r"^(putative|probable|predicted)\s+", "", txt, flags=re.IGNORECASE)
    return re.sub(r"\s+", " ", txt).strip(" ,;")


def _is_symbol_token(token: str) -> bool:
    """A word that is itself a gene symbol: MYLIP, DPY19L1, Mcm7, Runt, SgAbd-1.
    A designator such as S14, B9 or IIH is not one, nor is a word that only
    carries a chemical locant (O-methyltransferase, C-mannosyltransferase)."""
    core = re.sub(r"[^A-Za-z0-9]", "", re.sub(r"^[A-Za-z]-", "", token))
    letters = re.sub(r"[^A-Za-z]", "", core)
    if len(core) < 3 or len(letters) < 2 or not re.search(r"[A-Z]", core):
        return False
    if _ROMAN_RE.match(re.sub(r"[A-H]$", "", core)) or core in _BIOCHEM_ABBREVIATIONS:
        return False
    return True


def _is_filler(part: str) -> bool:
    # "a" is the article; "A" is a designator (protein kinase A).
    return part.lower() in _FILLER and not (len(part) == 1 and part.isupper())


def _is_kept_token(part: str) -> bool:
    """A piece that goes into an abbreviation whole: a number, a single letter, a
    Roman numeral, or something that already is an abbreviation (AMP, CD109, gp1)."""
    if len(part) == 1 or re.search(r"\d", part) or _ROMAN_RE.match(part):
        return True
    return part.isupper() and len(part) <= 5


def abbreviate_product(product: str) -> str:
    """Abbreviation of a product name, built by fixed rules ('' when it names nothing).

    1. A gene symbol at the end of the name is used as it stands
       ("E3 ubiquitin-protein ligase MYLIP" -> MYLIP, "segmentation protein Runt" -> Runt).
    2. A short name with one real word is kept whole: up to MAX_WHOLE_NAME_CHARS
       characters when it is a single word ("protein lozenge" -> lozenge,
       "formin-2", "junctophilin-1"), up to MAX_LABEL_CHARS with a designator
       after it ("annexin B9").
    3. Otherwise the pieces are joined in upper case. The name is split at spaces,
       hyphens and slashes; filler words (protein, of, and, ...) are dropped;
       numbers, single letters and existing abbreviations (AMP, CD109) are kept
       whole; every other word gives its first letter, "like" gives L
       ("pancreatic triacylglycerol lipase" -> PTL, "... protein 629-like" -> ZF629L).
       If that leaves fewer than MIN_ACRONYM_CHARS characters and the name has only
       one or two words, the first word gives its first three letters instead
       ("arginine kinase" -> ARGK, "cytochrome b5" -> CYTB5).
    """
    name = clean_product(product)
    if is_uninformative_product(name):
        return ""

    words = [w.strip(",") for w in name.split(" ")]
    kept_words = [w for w in words if w and not _is_filler(w)]
    if len(kept_words) > 1:
        last = kept_words[-1]
        base = re.sub(r"-like$", "", last, flags=re.IGNORECASE)
        if _is_symbol_token(base):
            return last if len(last) <= MAX_LABEL_CHARS else (base + "L")[:MAX_LABEL_CHARS]

    parts = [re.sub(r"[^A-Za-z0-9]", "", p) for p in re.split(r"[\s/,\-]+", name)]
    parts = [p for p in parts if p and not _is_filler(p)]
    if not parts:
        return ""
    real = [p for p in parts if not _is_kept_token(p) and p.lower() != "like"]

    whole = " ".join(kept_words)
    whole_limit = MAX_LABEL_CHARS if " " in whole else MAX_WHOLE_NAME_CHARS
    if len(whole) <= whole_limit and (len(real) == 1 or len(parts) == 1):
        return whole if len(whole) >= 3 or len(name) > MAX_WHOLE_NAME_CHARS else name

    def joined(expand_first):
        out, first = [], expand_first
        for p in parts:
            if _is_kept_token(p):
                out.append(p.upper())
            elif p.lower() == "like":
                out.append("L")
            else:
                out.append(p[:3].upper() if first else p[0].upper())
                first = False
        return "".join(out)

    label = joined(False)
    if len(label) < MIN_ACRONYM_CHARS and len(real) <= 2:
        label = joined(True)
    label = label[:MAX_LABEL_CHARS]
    return label if len(label) >= 3 or len(name) > MAX_WHOLE_NAME_CHARS else name


# ----------------------------------------------------------------------
# Home GFF
# ----------------------------------------------------------------------

def _attrs(field: str) -> Dict[str, str]:
    out = {}
    for kv in field.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k.strip()] = unquote(v.strip())
    return out


def load_gff_gene_info(gff_path: Optional[str],
                       wanted: Optional[Iterable[str]] = None) -> Dict[str, dict]:
    """Naming attributes of the genes in a GFF: ``{gene_id: {symbol, product, ncbi_gene_id}}``.

    ``symbol`` is '' when the GFF only repeats an identifier (``gene=LOC726817``,
    a locus tag). ``product`` is the gene's ``description`` or, failing that, the
    ``product`` of one of its transcripts. ``wanted`` limits the result to those
    gene IDs (the file is still read once).
    """
    info: Dict[str, dict] = {}
    if not gff_path or gff_path == "NO_GFF" or not os.path.exists(gff_path):
        return info
    wanted = set(wanted) if wanted is not None else None
    parent_of: Dict[str, str] = {}
    products: Dict[str, List[str]] = defaultdict(list)
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            ftype = p[2]
            if ftype in ("gene", "pseudogene"):
                a = _attrs(p[8])
                gid = a.get("ID", "")
                if not gid or (wanted is not None and gid not in wanted):
                    continue
                locus_tag = a.get("locus_tag", "")
                symbol = ""
                for key in _GFF_SYMBOL_KEYS:
                    cand = a.get(key, "").strip()
                    if cand and cand != locus_tag and not is_generic_id(cand):
                        symbol = cand
                        break
                m = re.search(r"GeneID:(\d+)", a.get("Dbxref", ""))
                info[gid] = {
                    "symbol": symbol,
                    "product": next((a[k] for k in _GFF_PRODUCT_KEYS if a.get(k)), ""),
                    "ncbi_gene_id": m.group(1) if m else "",
                }
            elif ftype in _TRANSCRIPT_TYPES:
                a = _attrs(p[8])
                parent = a.get("Parent", "").split(",")[0]
                if not parent:
                    continue
                if ftype != "CDS" and a.get("ID"):
                    parent_of[a["ID"]] = parent
                if a.get("product"):
                    products[parent].append(a["product"])
    for parent, names in products.items():
        gid = parent_of.get(parent, parent)
        entry = info.get(gid)
        if entry is None or (entry["product"] and not is_uninformative_product(entry["product"])):
            continue
        informative = [n for n in names if not is_uninformative_product(n)]
        if informative or not entry["product"]:
            entry["product"] = (informative or names)[0]
    return info


# ----------------------------------------------------------------------
# NCBI Gene
# ----------------------------------------------------------------------

def fetch_ncbi_gene_names(ncbi_gene_ids: Iterable[str],
                          timeout: float = NCBI_TIMEOUT_S) -> Dict[str, Tuple[str, str]]:
    """Current ``(symbol, description)`` of NCBI Gene records, keyed by GeneID.

    Raises ``OSError`` when NCBI cannot be reached or answers with something
    other than a gene report, so the caller can say that the lookup failed.
    """
    ids = sorted({str(i) for i in ncbi_gene_ids if str(i).isdigit()}, key=int)
    out: Dict[str, Tuple[str, str]] = {}
    for i in range(0, len(ids), NCBI_CHUNK):
        chunk = ids[i:i + NCBI_CHUNK]
        url = f"{NCBI_GENE_URL}{','.join(chunk)}?page_size={NCBI_CHUNK}"
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        api_key = os.environ.get("NCBI_API_KEY", "").strip()
        if api_key:
            req.add_header("api-key", api_key)
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                payload = json.load(resp)
        except (urllib.error.URLError, ValueError, TimeoutError) as exc:
            raise OSError(f"NCBI Gene lookup failed: {exc}") from exc
        for report in payload.get("reports", []) or []:
            gene = report.get("gene") or {}
            gid = str(gene.get("gene_id", ""))
            if gid:
                out[gid] = (str(gene.get("symbol", "") or ""),
                            str(gene.get("description", "") or ""))
    return out


# ----------------------------------------------------------------------
# Names table
# ----------------------------------------------------------------------

def read_names_tsv(path: Optional[str]) -> Dict[str, dict]:
    """Read a names table. A row needs ``gene_id`` and ``label``; ``source`` defaults
    to ``user``, so a hand-written two-column table is enough."""
    table: Dict[str, dict] = {}
    if not path or not os.path.exists(path):
        return table
    header = None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if not line.strip() or line.startswith("#"):
                continue
            cells = line.split("\t")
            if header is None and cells[0].strip().lower() == "gene_id":
                header = [c.strip().lower() for c in cells]
                continue
            cols = header or list(TSV_COLUMNS)
            row = {c: (cells[i].strip() if i < len(cells) else "") for i, c in enumerate(cols)}
            if not row.get("gene_id") or not row.get("label"):
                continue
            source = row.get("source") or SOURCE_USER
            table[row["gene_id"]] = {
                "label": row["label"],
                "full_name": row.get("full_name", ""),
                "source": source if source in SOURCES else SOURCE_USER,
                "ncbi_gene_id": row.get("ncbi_gene_id", ""),
            }
    return table


def write_names_tsv(path: str, table: Dict[str, dict], order: Iterable[str],
                    lookup_note: str = "") -> None:
    """Write the names table; the file can be edited and passed back with
    ``--gene_names_tsv``."""
    lines = [
        "# SynVoy gene names: the label each home gene carries in the figures.",
        "# Edit 'label' / 'full_name' and set 'source' to 'user' to fix a name, then pass",
        "# this file back with --gene_names_tsv. Rows with source 'product_name' are",
        f"# abbreviations made by SynVoy and are marked with '{DERIVED_MARK}' in the figures.",
    ]
    if lookup_note:
        lines.append(f"# NCBI Gene lookup: {lookup_note}")
    lines.append("\t".join(TSV_COLUMNS))
    seen = set()
    for gid in list(order) + sorted(table):
        if gid in seen or gid not in table:
            continue
        seen.add(gid)
        row = table[gid]
        lines.append("\t".join(str(row.get(c, "") if c != "gene_id" else gid).replace("\t", " ")
                               for c in TSV_COLUMNS))
    tmp = f"{path}.tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    os.replace(tmp, path)


def _full_name(product: str) -> str:
    """Product name for the table: '' when it does not name the gene."""
    return "" if is_uninformative_product(product) else clean_product(product)


def _bare_id(gene_id: str) -> str:
    for prefix in ("gene-", "gene:"):
        if gene_id.startswith(prefix):
            return gene_id[len(prefix):]
    return gene_id


def resolve_gene_names(gene_ids: Iterable[str], gff_info: Dict[str, dict],
                       overrides: Optional[Dict[str, dict]] = None,
                       allow_network: bool = True,
                       fetch=fetch_ncbi_gene_names) -> Tuple[Dict[str, dict], str]:
    """Name every gene in ``gene_ids`` (home genomic order).

    Returns ``(table, lookup_note)``. ``table[gene_id]`` has ``label``,
    ``full_name``, ``source`` and ``ncbi_gene_id``. ``lookup_note`` says what the
    NCBI Gene lookup did: it runs only for genes neither the names table nor the
    GFF names, and a failure leaves those genes to the GFF product name.
    """
    overrides = overrides or {}
    order = list(dict.fromkeys(gene_ids))
    table: Dict[str, dict] = {}
    pending = []
    for gid in order:
        if gid in overrides:
            table[gid] = dict(overrides[gid])
            continue
        entry = gff_info.get(gid) or {}
        if entry.get("symbol"):
            table[gid] = {"label": entry["symbol"], "source": SOURCE_GFF,
                          "full_name": _full_name(entry.get("product", "")),
                          "ncbi_gene_id": entry.get("ncbi_gene_id", "")}
        else:
            pending.append(gid)

    ncbi: Dict[str, Tuple[str, str]] = {}
    wanted = sorted({gff_info[g]["ncbi_gene_id"] for g in pending
                     if (gff_info.get(g) or {}).get("ncbi_gene_id")}, key=int)
    if not wanted:
        note = "not needed (no unnamed gene carries an NCBI GeneID)"
    elif not allow_network:
        note = f"skipped (--no_network); {len(wanted)} genes left to the GFF"
    else:
        try:
            ncbi = fetch(wanted)
            note = f"ok ({len(ncbi)} of {len(wanted)} records returned)"
        except OSError as exc:
            note = f"FAILED ({exc}); {len(wanted)} genes left to the GFF"

    for gid in pending:
        entry = gff_info.get(gid) or {}
        ncbi_id = entry.get("ncbi_gene_id", "")
        symbol, description = ncbi.get(ncbi_id, ("", ""))
        product = description if not is_uninformative_product(description) \
            else entry.get("product", "")
        # A symbol that only repeats the gene ID (a yeast ORF name) names nothing.
        if symbol and not is_generic_id(symbol) and symbol.lower() != _bare_id(gid).lower():
            table[gid] = {"label": symbol, "source": SOURCE_NCBI,
                          "full_name": _full_name(product), "ncbi_gene_id": ncbi_id}
            continue
        label = abbreviate_product(product)
        if label:
            table[gid] = {"label": label, "source": SOURCE_PRODUCT,
                          "full_name": clean_product(product), "ncbi_gene_id": ncbi_id}
        else:
            table[gid] = {"label": _bare_id(gid), "source": SOURCE_ID,
                          "full_name": "", "ncbi_gene_id": ncbi_id}

    _number_duplicates(table, [g for g in order if g not in overrides])
    return table, note


def _number_duplicates(table: Dict[str, dict], order: List[str]) -> None:
    """Two genes must not share a label: derived labels that collide (with each other
    or with a symbol) get -1, -2, ... in home genomic order."""
    by_label = defaultdict(list)
    for gid in order:
        by_label[table[gid]["label"].lower()].append(gid)
    for gids in by_label.values():
        derived = [g for g in gids if table[g]["source"] == SOURCE_PRODUCT]
        if len(gids) < 2 or not derived:
            continue
        for n, gid in enumerate(derived, start=1):
            base = table[gid]["label"][:MAX_WHOLE_NAME_CHARS - 2]
            table[gid]["label"] = f"{base}-{n}"


def display_label(row: Optional[dict], fallback: str = "") -> str:
    """Label as drawn: derived labels carry ``DERIVED_MARK``."""
    if not row:
        return fallback
    label = row.get("label") or fallback
    return label + DERIVED_MARK if row.get("source") == SOURCE_PRODUCT else label


# ----------------------------------------------------------------------
# The "Gene names" table drawn under every figure
# ----------------------------------------------------------------------

TABLE_FONT = 10.5
TABLE_ROW_H = 15
TABLE_NAME_MAX_W = 300     # a longer full name is cut; the hover title has all of it
_TABLE_GAP = 12
_TABLE_COL_GAP = 34
_TABLE_SWATCH = 9


def table_row(named: Optional[dict], bare_id: str, colour: str, is_goi: bool = False) -> dict:
    """One row of the gene-names table. ``named`` is the gene's row of the names
    table, or None for a GOI drawn from the query hit alone."""
    if named is None:
        return {"label": "GOI", "full": "query hit; no annotated home gene", "id": "",
                "source": "", "derived": False, "named": False, "is_goi": True,
                "colour": colour}
    return {
        "label": display_label(named),
        "full": named.get("full_name", ""),
        "id": "" if bare_id == named.get("label") else bare_id,
        "source": SOURCE_TEXT.get(named.get("source"), ""),
        "derived": named.get("source") == SOURCE_PRODUCT,
        "named": named.get("source") != SOURCE_ID,
        "is_goi": is_goi,
        "colour": colour,
    }


def _clip_to_width(text: str, size: float, max_w: float, text_width) -> str:
    if text_width(text, size) <= max_w:
        return text
    while text and text_width(text + "…", size) > max_w:
        text = text[:-1]
    return text.rstrip() + "…"


def _darken(hexc: str, factor: float = 0.65) -> str:
    h = hexc.lstrip("#")
    if len(h) != 6:
        return "#555555"
    return "#" + "".join(f"{int(int(h[i:i + 2], 16) * factor):02x}" for i in (0, 2, 4))


def names_table_svg(rows: List[dict], x: float, y: float, avail_w: float, text_width,
                    lookup_note: str = "", goi_stroke: str = "#b91c1c") -> Tuple[List[str], float]:
    """SVG of the gene-names table with its top-left corner at (x, y).

    One entry per home gene: colour swatch, label, full name, gene ID, where the
    label comes from. Entries run down as many columns as fit into ``avail_w``.
    ``text_width(text, size, bold=False)`` measures text in the figure's font.
    Returns ``(svg_parts, height)``; nothing is drawn (height 0) when no gene has
    a name, since a list of bare IDs explains nothing.
    """
    if not any(r["named"] for r in rows):
        return [], 0
    fs, rh, gap = TABLE_FONT, TABLE_ROW_H, _TABLE_GAP
    w_label = max(text_width(r["label"], fs, True) for r in rows)
    w_id = max(text_width(r["id"], fs) for r in rows)
    w_src = max(text_width(r["source"], fs) for r in rows)
    # One entry has to fit the canvas: the full name gives way first.
    fixed = _TABLE_SWATCH + 6 + w_label + w_id + w_src + gap * 3
    w_full = max(0.0, min(TABLE_NAME_MAX_W, max(text_width(r["full"], fs) for r in rows),
                          avail_w - fixed))
    widths = [w for w in (w_label, w_full, w_id, w_src) if w > 0]
    entry_w = _TABLE_SWATCH + 6 + sum(widths) + gap * (len(widths) - 1)
    n_cols = max(1, min(len(rows), int((avail_w + _TABLE_COL_GAP) // (entry_w + _TABLE_COL_GAP))))
    per_col = -(-len(rows) // n_cols)

    head = "Gene names"
    sub = "label   ·   name   ·   gene ID   ·   where the label comes from"
    notes = []
    if any(r["derived"] for r in rows):
        notes.append(DERIVED_NOTE)
    if lookup_note.startswith(("FAILED", "skipped")):
        notes.append("NCBI Gene was not queried"
                     + (" (not reachable)" if lookup_note.startswith("FAILED") else "")
                     + ": genes without a symbol in the home GFF carry an abbreviation "
                       "or their ID")

    out = ['<g class="gene-names-table">']
    out.append(f'<line x1="{x:.1f}" y1="{y + 6:.1f}" x2="{x + avail_w:.1f}" y2="{y + 6:.1f}" '
               f'stroke="#e1e6ef" stroke-width="1"/>')
    y0 = y + 24
    out.append(f'<text x="{x:.1f}" y="{y0:.1f}" font-size="12" font-weight="700" '
               f'fill="#1a1d26">{head}</text>')
    out.append(f'<text x="{x + text_width(head, 12, True) + 14:.1f}" y="{y0:.1f}" '
               f'font-size="9.5" fill="#8c95a6">{_esc(sub)}</text>')
    top = y0 + 8
    for k, r in enumerate(rows):
        cx = x + (k // per_col) * (entry_w + _TABLE_COL_GAP)
        cy = top + (k % per_col) * rh + rh - 4
        title = " · ".join(t for t in (r["label"], r["full"], r["id"], r["source"]) if t)
        cell = [f'<title>{_esc(title)}</title>',
                f'<rect x="{cx:.1f}" y="{cy - _TABLE_SWATCH + 1:.1f}" width="{_TABLE_SWATCH}" '
                f'height="{_TABLE_SWATCH}" rx="2" fill="{r["colour"]}" '
                f'stroke="{goi_stroke if r["is_goi"] else _darken(r["colour"])}" '
                f'stroke-width="{1.2 if r["is_goi"] else 0.6}"/>']
        tx = cx + _TABLE_SWATCH + 6
        cell.append(f'<text x="{tx:.1f}" y="{cy:.1f}" font-size="{fs}" font-weight="700" '
                    f'fill="{goi_stroke if r["is_goi"] else "#1a1d26"}">{_esc(r["label"])}</text>')
        tx += w_label + gap
        if w_full:
            cell.append(f'<text x="{tx:.1f}" y="{cy:.1f}" font-size="{fs}" fill="#42495a">'
                        f'{_esc(_clip_to_width(r["full"], fs, w_full, text_width))}</text>')
            tx += w_full + gap
        if w_id:
            cell.append(f'<text x="{tx:.1f}" y="{cy:.1f}" font-size="{fs}" fill="#6b7280">'
                        f'{_esc(r["id"])}</text>')
            tx += w_id + gap
        if w_src:
            cell.append(f'<text x="{tx:.1f}" y="{cy:.1f}" font-size="{fs}" fill="#8c95a6" '
                        f'font-style="italic">{_esc(r["source"])}</text>')
        out.append("<g>" + "".join(cell) + "</g>")
    bottom = top + per_col * rh
    for i, note in enumerate(notes):
        out.append(f'<text x="{x:.1f}" y="{bottom + 16 + i * 13:.1f}" font-size="9.5" '
                   f'fill="#8c95a6">{_esc(note)}</text>')
    out.append('</g>')
    return out, (bottom + 12 + len(notes) * 13 + 10) - y
