"""Taxonomy of a species from NCBI, read once when a genome is added to the library.

The lineage (taxon ids from the root down to the species) is what lets the search form
list "closest relatives first". Nothing here runs while visitors use the site.
"""

import json
import xml.etree.ElementTree as ElementTree
from urllib.parse import quote

from .net import FetchError, http_get

_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def lookup_taxonomy(species, taxon_id=None, fetch=http_get):
    """Return ``(taxon_id, lineage, common_name)`` for a species name (or a known taxon
    id). ``lineage`` runs from the root to the species itself. Raises FetchError."""
    if not taxon_id:
        term = quote('"%s"[Scientific Name]' % species, safe="")
        try:
            found = json.loads(fetch(_EUTILS + "esearch.fcgi?db=taxonomy&retmode=json&term=" + term,
                                     accept="application/json").decode("utf-8"))
        except ValueError:
            raise FetchError("unavailable", "NCBI sent something that is not JSON")
        ids = ((found or {}).get("esearchresult") or {}).get("idlist") or []
        if not ids:
            raise FetchError("not_found", "NCBI Taxonomy has no species named %r" % species)
        taxon_id = int(ids[0])
    data = fetch(_EUTILS + "efetch.fcgi?db=taxonomy&retmode=xml&id=%d" % int(taxon_id))
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
    lineage.append(int(taxon_id))
    common = (taxon.findtext("OtherNames/GenbankCommonName")
              or taxon.findtext("OtherNames/CommonName") or "").strip()
    return int(taxon_id), lineage, common
