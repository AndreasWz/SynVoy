"""Taxonomy of a species from NCBI, read once when a genome is added to the library or
fetched into the genome store.

The lineage (taxon ids from the top down to the species) is what lets the search form
list "closest relatives first". It is read from NCBI's Datasets service; the older
E-utilities are the second source. Both give the same list (compared on 2026-10-11 for
a yeast, a yeast strain, the honeybee and a nematode), and the second one answers
"HTTP 400" every now and then for a request that works a second later, so it is asked
up to three times.
"""

import json
import time
import xml.etree.ElementTree as ElementTree
from urllib.parse import quote

from .net import FetchError, http_get

_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
_DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/"
_TRIES = 3
_PAUSE = 0.8          # seconds between tries


def _passing(exc):
    """Is this the kind of failure that the same request does not show a moment later?"""
    return exc.kind == "unavailable" or (exc.kind == "not_found" and "HTTP 400" in exc.detail)


def _eutils(fetch, address, sleep, accept="*/*"):
    for attempt in range(1, _TRIES + 1):
        try:
            return fetch(_EUTILS + address, timeout=10, accept=accept)
        except FetchError as exc:
            if attempt == _TRIES or not _passing(exc):
                raise
            sleep(_PAUSE)


def _from_datasets(taxon_id, fetch):
    try:
        data = json.loads(fetch(_DATASETS + "%d" % taxon_id, timeout=10,
                                accept="application/json").decode("utf-8"))
        taxonomy = data["taxonomy_nodes"][0]["taxonomy"]
        lineage = [int(x) for x in taxonomy["lineage"] if int(x) != 1]     # 1 is the root of all
    except (ValueError, KeyError, IndexError, TypeError):
        raise FetchError("not_found", "NCBI Datasets has no taxon %s" % taxon_id)
    if not lineage or lineage[-1] != taxon_id:
        lineage.append(taxon_id)
    common = taxonomy.get("genbank_common_name") or taxonomy.get("common_name") or ""
    return taxon_id, lineage, str(common).strip()


def _from_eutils(taxon_id, fetch, sleep):
    data = _eutils(fetch, "efetch.fcgi?db=taxonomy&retmode=xml&id=%d" % taxon_id, sleep)
    try:
        root = ElementTree.fromstring(data)
    except ElementTree.ParseError:
        raise FetchError("unavailable", "NCBI sent something that is not XML")
    taxon = root.find("Taxon")
    if taxon is None:
        raise FetchError("not_found", "NCBI Taxonomy has no taxon %s" % taxon_id)
    lineage = []
    for entry in taxon.findall("LineageEx/Taxon"):
        try:
            lineage.append(int(entry.findtext("TaxId", "")))
        except ValueError:
            continue
    lineage.append(taxon_id)
    common = (taxon.findtext("OtherNames/GenbankCommonName")
              or taxon.findtext("OtherNames/CommonName") or "").strip()
    return taxon_id, lineage, common


def lookup_taxonomy(species, taxon_id=None, fetch=http_get, sleep=time.sleep):
    """Return ``(taxon_id, lineage, common_name)`` for a species name (or a known taxon
    id). ``lineage`` runs from the top down to the species itself. Raises FetchError."""
    if not taxon_id:
        term = quote('"%s"[Scientific Name]' % species, safe="")
        try:
            found = json.loads(_eutils(fetch, "esearch.fcgi?db=taxonomy&retmode=json&term=" + term,
                                       sleep, accept="application/json").decode("utf-8"))
        except ValueError:
            raise FetchError("unavailable", "NCBI sent something that is not JSON")
        ids = ((found or {}).get("esearchresult") or {}).get("idlist") or []
        if not ids:
            raise FetchError("not_found", "NCBI Taxonomy has no species named %r" % species)
        taxon_id = ids[0]
    taxon_id = int(taxon_id)
    try:
        return _from_datasets(taxon_id, fetch)
    except FetchError:
        return _from_eutils(taxon_id, fetch, sleep)
