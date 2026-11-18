"""The figures name the GOI after a home gene only when the query hits fall on that gene's CDS.

Without an annotated gene under the hits' CDS (a curated seed, an ORF, a gene in an intron) the
home GOI is the gene the span overlap finds (identify_goi_names), and that can be an lncRNA over
the hit. The home row's placeholder and the GOI's own id were already guarded (freeze fix e3);
the GOI column, the figure titles and the legends (_goi_display_label) and the ribbon plot's
home label still took that gene's name. Now they take the GOI's own name, as its target models
carry it, and the annotated gene names the GOI only on the CDS path.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import plot_synteny as ps  # noqa: E402

ROWS = {"gene-NB1": {"label": "NB1", "source": "gff_symbol", "full_name": "", "ncbi_gene_id": ""},
        "gene-lncX": {"label": "lncX", "source": "gff_symbol", "full_name": "long non-coding RNA X",
                      "ncbi_gene_id": ""},
        "gene-NB2": {"label": "NB2", "source": "gff_symbol", "full_name": "", "ncbi_gene_id": ""}}


@pytest.fixture
def state(monkeypatch):
    monkeypatch.setattr(ps, "_GENE_NAMES", dict(ROWS))
    monkeypatch.setattr(ps, "_GOI_NAMES", set())
    monkeypatch.setattr(ps, "_GOI_FROM_CDS", [False])
    monkeypatch.setattr(ps, "_GOI_QUERY_NAME", [None])
    monkeypatch.setattr(ps, "_GOI_DISPLAY_LABEL", [None])


def _home():
    return [{"chrom": "c", "start": 1000, "end": 2000, "name": "gene-NB1", "display_name": "NB1"},
            # an lncRNA over the hit span: no CDS, in the home BED as a neighbour
            {"chrom": "c", "start": 4800, "end": 6200, "name": "gene-lncX", "display_name": "lncX"},
            {"chrom": "c", "start": 9000, "end": 9900, "name": "gene-NB2", "display_name": "NB2"}]


HITS = [{"chrom": "c", "start": 5000, "end": 5600, "strand": "+"}]


def _tracks(home, target_names=("GOI_SEED7",)):
    genes = [{"name": f"{n}|Apis_florea_fna_b0_l1_exon_ann", "home_gene_id": n} for n in target_names]
    return [{"is_home": True, "genes": home}, {"is_home": False, "genes": genes}]


def test_a_gene_the_span_overlap_found_does_not_name_the_goi(state):
    home = _home()
    ps.set_goi_names(home, HITS, [])              # nothing under the hits' CDS: the span rule
    assert ps._GOI_NAMES == {"gene-lncX"} and ps._GOI_FROM_CDS[0] is False
    assert ps._goi_display_label(_tracks(home)) == "SEED7"


def test_the_gene_under_the_hits_cds_names_the_goi(state):
    home = _home()
    ps.set_goi_names(home, HITS, ["gene-NB1"])
    assert ps._GOI_FROM_CDS[0] is True
    assert ps._goi_display_label(_tracks(home)) == "NB1"


def test_without_a_model_name_the_goi_is_goi(state):
    home = _home()
    ps.set_goi_names(home, HITS, [])
    assert ps._goi_display_label(_tracks(home, target_names=())) == "GOI"


def test_the_ribbon_plot_names_the_goi_as_the_other_figures(state, monkeypatch):
    home = [dict(g) for g in _home()]
    ps.set_goi_names(home, HITS, [])
    tracks = _tracks(home)
    monkeypatch.setattr(ps, "_GOI_DISPLAY_LABEL", [ps._goi_display_label(tracks)])
    (gene, text), = ps._track_label_jobs(tracks[0], {})
    assert gene["name"] == "gene-lncX" and text == "★ SEED7"
    # on the CDS path the gene's own label, as before
    monkeypatch.setattr(ps, "_GOI_NAMES", {"gene-NB1"})
    monkeypatch.setattr(ps, "_GOI_FROM_CDS", [True])
    monkeypatch.setattr(ps, "_GOI_DISPLAY_LABEL", [ps._goi_display_label(tracks)])
    (gene, text), = ps._track_label_jobs(tracks[0], {})
    assert gene["name"] == "gene-NB1" and text == "★ NB1"
