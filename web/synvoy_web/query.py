"""What a visitor may search with: a pasted protein sequence or a protein accession.

The sequence rules are the pipeline's own (``bin/normalize_query.py``): case is
folded, whitespace, digits, alignment gaps and one trailing ``*`` are formatting,
anything else that is not an amino-acid letter is refused. One thing is stricter
here: DNA is refused instead of translated, because the site asks for a protein and a
silent translation of the wrong strand would be searched without anyone noticing.
``tests/test_query.py`` checks that the alphabets still agree with the pipeline.
"""

import json
import re
import time
from urllib.parse import quote

from .net import FetchError, http_get

PROTEIN_ALPHABET = frozenset("ACDEFGHIKLMNPQRSTVWYBJOUXZ")
NUC_ALPHABET = frozenset("ACGTUNRYWSKMBDHV")
_FORMATTING = frozenset("-.")

# https://www.uniprot.org/help/accession_numbers
UNIPROT_RE = re.compile(r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$")
UNIPROT_ISOFORM_RE = re.compile(r"^(.+)-[0-9]{1,3}$")
NCBI_REFSEQ_RE = re.compile(r"^[XNWYA]P_[0-9]{5,10}(?:\.[0-9]{1,3})?$")
NCBI_GENBANK_RE = re.compile(r"^[A-Z]{3}[0-9]{5,8}(?:\.[0-9]{1,3})?$")

MAX_PASTE_CHARS = 60000
NAME_RE = re.compile(r"[^A-Za-z0-9_.\-]+")
LOOKUP_TTL_DAYS = 30


class QueryError(Exception):
    """The request cannot be accepted as it is. ``message`` is written for the visitor."""

    def __init__(self, code, message, field=None, status=400):
        Exception.__init__(self, message)
        self.code = code
        self.message = message
        self.field = field
        self.status = status


def safe_name(text, fallback="query", limit=40):
    """A label that is safe in a FASTA header, a file name and a figure."""
    name = NAME_RE.sub("_", (text or "").strip()).strip("_.-")[:limit].strip("_.-")
    return name or fallback


def is_nucleotide(seq):
    if not seq or not all(ch in NUC_ALPHABET for ch in seq):
        return False
    plain = sum(1 for ch in seq if ch in "ACGTUN")
    return plain >= 0.9 * len(seq)


def clean_protein(text, min_length, max_length):
    """Return ``(header, sequence)`` for pasted text, or raise QueryError.

    ``header`` is the text of the FASTA ``>`` line without the ``>`` ('' if none).
    """
    if not isinstance(text, str):
        raise QueryError("sequence_missing", "Paste a protein sequence.", "sequence")
    if len(text) > MAX_PASTE_CHARS:
        raise QueryError("sequence_too_long",
                         "The pasted text is too long. Paste one protein sequence.", "sequence")
    headers, body = [], []
    for line in text.replace("\r", "\n").split("\n"):
        line = line.strip()
        if line.startswith(">"):
            headers.append(line[1:].strip())
        elif line:
            body.append(line)
    if len(headers) > 1:
        raise QueryError("several_sequences",
                         "The text holds %d sequences. Paste one protein sequence." % len(headers),
                         "sequence")
    raw = "".join(body).upper()
    seq = "".join(ch for ch in raw if not (ch.isspace() or ch.isdigit() or ch in _FORMATTING))
    if seq.endswith("*"):
        seq = seq[:-1]
    if not seq:
        raise QueryError("sequence_missing", "Paste a protein sequence.", "sequence")
    bad = sorted({ch for ch in seq if ch not in PROTEIN_ALPHABET})
    if bad:
        shown = " ".join(bad[:6])
        if "*" in bad:
            message = ("The sequence has a stop sign (*) in the middle. "
                       "Paste one complete protein without stops.")
        else:
            message = "The sequence has characters that are not amino-acid letters: %s" % shown
        raise QueryError("sequence_invalid", message, "sequence")
    if is_nucleotide(seq):
        raise QueryError("sequence_is_dna",
                         "This looks like DNA or RNA. Paste the protein sequence instead.",
                         "sequence")
    if len(seq) < min_length:
        raise QueryError("sequence_too_short",
                         "The sequence has %d residues. The shortest accepted is %d."
                         % (len(seq), min_length), "sequence")
    if len(seq) > max_length:
        raise QueryError("sequence_too_long",
                         "The sequence has %d residues. The longest accepted is %d."
                         % (len(seq), max_length), "sequence")
    return (headers[0] if headers else ""), seq


def classify_accession(text):
    """Return ``(source, accession)`` with source 'uniprot' or 'ncbi', or raise."""
    if not isinstance(text, str):
        raise QueryError("accession_missing", "Enter a protein accession.", "accession")
    acc = text.strip().upper()
    if not acc:
        raise QueryError("accession_missing", "Enter a protein accession.", "accession")
    if len(acc) > 20:
        raise QueryError("accession_unknown",
                         "This is not a UniProt or NCBI protein accession.", "accession")
    if NCBI_REFSEQ_RE.match(acc) or NCBI_GENBANK_RE.match(acc):
        return "ncbi", acc
    isoform = UNIPROT_ISOFORM_RE.match(acc)
    if isoform and UNIPROT_RE.match(isoform.group(1)):
        raise QueryError("accession_isoform",
                         "Isoform accessions are not supported. Use %s, or paste the "
                         "sequence of the isoform." % isoform.group(1), "accession")
    base = acc.split(".")[0]
    if UNIPROT_RE.match(base):
        return "uniprot", base
    raise QueryError("accession_unknown",
                     "This is not a UniProt or NCBI protein accession. Examples: P01501, "
                     "XP_031329057.1. You can also paste the sequence.", "accession")


def _one_line(text, limit=300):
    return re.sub(r"[\x00-\x1f\x7f]+", " ", text or "").strip()[:limit]


def _parse_fasta_record(text):
    """(header, sequence) of the single record in ``text``; ValueError otherwise."""
    lines = [l.strip() for l in text.replace("\r", "\n").split("\n") if l.strip()]
    if not lines or not lines[0].startswith(">"):
        raise ValueError("no FASTA header")
    if sum(1 for l in lines if l.startswith(">")) != 1:
        raise ValueError("more than one record")
    sequence = "".join(lines[1:]).upper().rstrip("*")
    if not sequence or any(ch not in PROTEIN_ALPHABET for ch in sequence):
        raise ValueError("not a protein sequence")
    return _one_line(lines[0][1:], 500), sequence


class Resolver:
    """Turns an accession into name, organism and sequence. Answers are kept for
    ``LOOKUP_TTL_DAYS`` so a popular accession is fetched once, and so a search can
    still be started while UniProt or NCBI are unreachable."""

    def __init__(self, store, fetch=http_get, clock=time.time):
        self.store = store
        self.fetch = fetch
        self.clock = clock

    def cached(self, accession):
        return self.store.get_lookup(accession, int(self.clock()) - LOOKUP_TTL_DAYS * 86400)

    def resolve(self, text):
        """Return ``(record, from_cache)``. Raises QueryError with a message for the
        visitor when the accession is unknown or the database does not answer."""
        source, accession = classify_accession(text)
        record = self.cached(accession)
        if record is not None:
            return record, True
        try:
            if source == "uniprot":
                record = self._uniprot(accession)
            else:
                record = self._ncbi(accession)
        except FetchError as exc:
            where = "UniProt" if source == "uniprot" else "NCBI"
            if exc.kind == "not_found":
                raise QueryError("accession_not_found",
                                 "%s has no protein with the accession %s." % (where, accession),
                                 "accession", status=404)
            raise QueryError("lookup_unavailable",
                             "%s did not answer. Try again in a minute, or paste the "
                             "sequence instead." % where, "accession", status=503)
        self.store.put_lookup(accession, record, int(self.clock()))
        return record, False

    # ── UniProt ─────────────────────────────────────────────────────────────
    def _uniprot(self, accession):
        base = "https://rest.uniprot.org/uniprotkb/" + quote(accession, safe="")
        try:
            data = json.loads(self.fetch(base + ".json", accept="application/json").decode("utf-8"))
        except ValueError:
            raise FetchError("unavailable", "UniProt sent something that is not JSON")
        if not isinstance(data, dict):
            raise FetchError("unavailable", "unexpected answer from UniProt")
        if str(data.get("entryType", "")).lower().startswith("inactive"):
            raise QueryError("accession_inactive",
                             "The accession %s is no longer in UniProt (it was removed or "
                             "merged into another entry)." % accession, "accession", status=404)
        sequence = (data.get("sequence") or {}).get("value") or ""
        if not sequence:
            raise FetchError("not_found", "entry without a sequence")
        try:
            header, fasta_sequence = _parse_fasta_record(self.fetch(base + ".fasta").decode("utf-8"))
        except ValueError as exc:
            raise FetchError("unavailable", "unexpected FASTA from UniProt: %s" % exc)
        if fasta_sequence != sequence.upper():
            raise FetchError("unavailable", "UniProt entry and FASTA disagree")
        description = data.get("proteinDescription") or {}
        name = ((description.get("recommendedName") or {}).get("fullName") or {}).get("value", "")
        if not name:
            for item in description.get("submissionNames") or []:
                name = (item.get("fullName") or {}).get("value", "")
                if name:
                    break
        gene = ""
        for item in data.get("genes") or []:
            gene = (item.get("geneName") or {}).get("value", "")
            if gene:
                break
        organism = data.get("organism") or {}
        taxon = organism.get("taxonId")
        return {
            "source": "uniprot",
            "accession": str(data.get("primaryAccession") or accession),
            "name": _one_line(name, 160),
            "gene": _one_line(gene, 40),
            "organism": _one_line(organism.get("scientificName", ""), 120),
            "organism_common": _one_line(organism.get("commonName", ""), 80),
            "taxon_id": taxon if isinstance(taxon, int) else None,
            "length": len(fasta_sequence),
            "sequence": fasta_sequence,
            "header": header,
            "reviewed": "reviewed" in str(data.get("entryType", "")).lower()
                        and "unreviewed" not in str(data.get("entryType", "")).lower(),
        }

    # ── NCBI ────────────────────────────────────────────────────────────────
    def _ncbi(self, accession):
        base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        acc = quote(accession, safe="")
        try:
            summary = json.loads(self.fetch(
                base + "esummary.fcgi?db=protein&retmode=json&id=" + acc,
                accept="application/json").decode("utf-8"))
        except ValueError:
            raise FetchError("unavailable", "NCBI sent something that is not JSON")
        result = (summary or {}).get("result") or {}
        uids = result.get("uids") or []
        if not uids or "error" in summary:
            raise FetchError("not_found", "no such protein record")
        entry = result.get(str(uids[0])) or {}
        if entry.get("error"):
            raise FetchError("not_found", str(entry.get("error")))
        try:
            header, sequence = _parse_fasta_record(self.fetch(
                base + "efetch.fcgi?db=protein&rettype=fasta&retmode=text&id=" + acc).decode("utf-8"))
        except ValueError as exc:
            raise FetchError("not_found", "no protein FASTA for this accession: %s" % exc)
        title = _one_line(entry.get("title", ""), 200)
        organism = _one_line(entry.get("organism", ""), 120)
        name = re.sub(r"\s*\[[^\[\]]*\]\s*$", "", title)      # drop the "[Genus species]" suffix
        taxon = entry.get("taxid")
        return {
            "source": "ncbi",
            "accession": _one_line(str(entry.get("accessionversion") or accession), 30),
            "name": name,
            "gene": "",
            "organism": organism,
            "organism_common": "",
            "taxon_id": taxon if isinstance(taxon, int) and taxon > 0 else None,
            "length": len(sequence),
            "sequence": sequence,
            "header": header,
            "reviewed": False,
        }
