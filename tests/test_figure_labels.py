"""Review fix RF-3 (figures): a GOI call is drawn as resolved when its label is HIGH or MEDIUM,
whatever its class; an AMBIGUOUS call gets the needs-curation style; the best GOI call and the
overlap dedup rank by label, then by _model_strength. Four labels x two evidence types (a gene
model, whose class follows the label, and a tandem copy, whose class is always
tandem_goi_copy)."""
import html
import json
import os
import re
import sys
from types import SimpleNamespace

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import plot_synteny as ps  # noqa: E402

LABELS = ("HIGH", "MEDIUM", "AMBIGUOUS", "LOW")
EVIDENCE = ("exon_annotation", "tandem_copy")
# As the genome-end judge writes them (iterative_search_runner.GOI_CLASS_OF_LABEL).
CLASS_OF = {"HIGH": "confident_goi", "MEDIUM": "probable_goi",
            "AMBIGUOUS": "syntenic_candidate_unconfirmed", "LOW": "ambiguous_goi_family_member"}


def goi(label, evidence, start=250, end=330, eeff=1e-20, name="GOI_q"):
    return {"chrom": "chr1", "start": start, "end": end, "start_plot": start,
            "end_plot": end, "name": name, "home_gene_id": "GOI_q", "strand": "+",
            "identity": 60.0, "exon_coords": [], "role": "goi", "confidence": label,
            "goi_class": "tandem_goi_copy" if evidence == "tandem_copy" else CLASS_OF[label],
            "evidence_type": evidence, "eeff": eeff, "bits": 80.0, "block_esyn": 1e-4}


def tracks(target_gois):
    def gene(name, start, end, home_gene_id=None):
        return {"chrom": "chr1", "start": start, "end": end, "start_plot": start,
                "end_plot": end, "name": name, "home_gene_id": home_gene_id or name,
                "strand": "+", "identity": 75.0, "exon_coords": []}
    return [
        {"label": "Home species (chrH)", "is_home": True, "genome_id": "home",
         "goi_status": "resolved", "offset": 0, "breaks": [], "minus_strand_row": None,
         "genes": [gene("geneA", 100, 180), gene("GOI_q", 240, 320), gene("geneB", 580, 660)]},
        {"label": "Target species (target1)", "is_home": False, "genome_id": "target1",
         "goi_status": ps._track_goi_status({"genes": target_gois}), "offset": 0,
         "breaks": [], "minus_strand_row": None,
         "genes": [gene("orthA", 120, 200, home_gene_id="geneA")] + target_gois},
    ]


def target_goi_fills(target_gois):
    out = ps.render_synteny_html(tracks(target_gois), {"geneA": "#0072B2", "geneB": "#009E73"},
                                 {}, SimpleNamespace(plot_width=0, plot_height=0,
                                                     scale_bar_len=10000), ["test"], 0, 0, 1)
    fills = re.findall(r'data-is-goi="true"[^>]*>\s*<path d="[^"]+" fill="([^"]+)"', out)
    return fills[1:]          # the first one is the home GOI


@pytest.mark.parametrize("evidence", EVIDENCE)
@pytest.mark.parametrize("label", LABELS)
def test_resolved_style_and_tag_follow_the_label(label, evidence):
    g = goi(label, evidence)
    resolved = label in ("HIGH", "MEDIUM")
    assert ps._is_resolved_goi_target_gene(g) is resolved
    assert ps._track_goi_status({"genes": [g]}) == ("resolved" if resolved else "ambiguous")
    # Ribbon plot: red for a resolved call, the needs-curation pattern otherwise.
    assert target_goi_fills([g]) == [ps.GOI_COLOUR if resolved else "url(#ambiguousGoi)"]
    tag = json.loads(html.unescape(ps._build_tooltip_json(g, tracks([g])[1], {})))["goiTag"]
    assert tag == {"HIGH": "GOI", "MEDIUM": "GOI", "AMBIGUOUS": "GOI call, needs curation",
                   "LOW": "GOI-like, low confidence"}[label]
    # The anchor grid's copy fill: AMBIGUOUS palest with the most broken outline.
    fill, dash = ps._goi_copy_fill(60.0, label, False)
    if label == "AMBIGUOUS":
        assert dash == ' stroke-dasharray="1.5,2.5"'
    assert (dash == "") is (label == "HIGH")


def test_best_goi_is_by_label_then_strength_whatever_the_class():
    strong_ambiguous = goi("AMBIGUOUS", "exon_annotation", eeff=1e-60, name="GOI_a")
    weak_medium_copy = goi("MEDIUM", "tandem_copy", start=400, end=480, eeff=1e-4, name="GOI_m")
    assert max([strong_ambiguous, weak_medium_copy], key=ps._goi_priority_key) is weak_medium_copy
    # Same label: the stronger match wins, a tandem copy against a gene model (E-values
    # below 1e-10 tie in _model_strength, as in the judge's order).
    model = goi("HIGH", "exon_annotation", eeff=1e-7, name="GOI_h")
    copy = goi("HIGH", "tandem_copy", start=400, end=480, eeff=1e-9, name="GOI_t")
    assert max([model, copy], key=ps._goi_priority_key) is copy
    track = tracks([strong_ambiguous, weak_medium_copy])[1]
    assert ps._track_label_jobs(track, {})[0][0] is weak_medium_copy
    assert ps._goi_best([model, copy]) is copy


def _gff_row(label, evidence, start, end, eeff):
    ftype = "gene" if evidence == "tandem_copy" else "mRNA"
    cls = "tandem_goi_copy" if evidence == "tandem_copy" else CLASS_OF[label]
    return (f"chr1\tsynvoy\t{ftype}\t{start}\t{end}\t60.0\t+\t.\tID=GOI_q|g_{label}_{evidence};"
            f"Name=GOI_q;SynVoy_Parent=GOI_q;SynVoyRole=goi;EvidenceType={evidence};"
            f"Confidence={label};GOIClass={cls};Identity=60.0;QueryCoverage=0.9;"
            f"Eeff={eeff};BlockEsyn=1e-4;Bits=80\n")


@pytest.mark.parametrize("evidence", EVIDENCE)
def test_overlap_dedup_keeps_the_better_label(tmp_path, evidence):
    other = "exon_annotation" if evidence == "tandem_copy" else "tandem_copy"
    gff = tmp_path / "t.fna.gff"
    # The AMBIGUOUS call has the far stronger E; the overlapping MEDIUM call still wins.
    gff.write_text(_gff_row("AMBIGUOUS", other, 1000, 1900, "1e-60")
                   + _gff_row("MEDIUM", evidence, 1100, 1950, "1e-4"))
    genes = ps.parse_target_gff(str(gff), cap=False)
    assert [(g["confidence"], g["evidence_type"]) for g in genes] == [("MEDIUM", evidence)]
