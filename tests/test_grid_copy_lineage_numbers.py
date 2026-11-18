"""The anchor grid says what a GOI model's identity is measured against.

With --tandem_copy_parent seed (the default) a single-segment tandem copy found in one genome
is added to the queries of the next wave as a GOI of its own, GOI_copy_<n>. The models it finds
carry SynVoy_Parent=GOI_copy_<n>, and their Identity, QueryCoverage and Bits are to that copy,
not to the query. The grid printed the number after such a model and said "% identity to the
query". Now the legend and the draft figure legend say which models are measured against a
copy, the source table names the reference of every value (identity_to), and the full grid's
tooltip names the copy. The numbers themselves are unchanged.
"""
import os
import sys
from types import SimpleNamespace

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import plot_synteny as ps  # noqa: E402
import synvoy_grid as sg  # noqa: E402


def _g(name, start, end, **kw):
    g = {"chrom": "chr1", "start": start, "end": end, "start_plot": start, "end_plot": end,
         "name": name, "strand": "+", "exon_coords": []}
    g.update(kw)
    return g


def _goi(start, parent="GOI_T", conf="HIGH", ident=70.0, cov=0.95):
    return _g(f"{parent}|m{start}", start, start + 300, role="goi", home_gene_id=parent,
              confidence=conf, identity=ident, query_coverage=cov, evidence_type="exon_annotation",
              exon_coords=[(start, start + 100), (start + 200, start + 300)])


def _tracks(target_goi):
    names = [f"gene-F{k}" for k in range(6)]
    home = [_g(n, 1000 + 2000 * k, 2000 + 2000 * k) for k, n in enumerate(names)]
    home.insert(3, _g("GOI_T", 6600, 6900))
    target = [_g(f"t_{k}", 1000 + 2000 * k, 2000 + 2000 * k, home_gene_id=n, identity=80.0,
                 confidence="HIGH", query_coverage=0.95) for k, n in enumerate(names)]
    tracks = [{"label": "Home", "species": "Apis mellifera", "is_home": True, "genome_id": "home",
               "offset": 0, "breaks": [], "genes": home},
              {"label": "T0", "species": "Bombus species0", "is_home": False, "genome_id": "g0",
               "offset": 0, "breaks": [], "genes": target + target_goi}]
    return tracks, {n: ps.GENE_PALETTE[k] for k, n in enumerate(names)}


def _numbers(groups):
    return [text for title, rows in groups if title == "Numbers" for _, text in rows]


def test_a_model_found_from_a_seeded_copy_is_measured_against_that_copy():
    calls = [_goi(7000), _goi(9000, parent="GOI_copy_2", ident=93.3, cov=1.0)]
    svg, legend, source, info = ps.render_anchor_grid_print(*_tracks(calls), {}, SimpleNamespace())
    assert info["copy_lineage_models"] == 1
    goi = {r["copy_number"]: r for r in source if r["row"] == 2 and r["column"] == "T"}
    assert goi[2]["identity_to"] == "GOI_copy_2" and goi[2]["identity_pct"] == "93.3"
    assert goi[1]["identity_to"] == "query"
    assert ("1 of these models was found from a tandem copy of T that the search of an earlier "
            "genome added to the queries (GOI_copy_<n>): its numbers are the identity and coverage "
            "to that copy, not to the query (source data, column identity_to).") in legend
    res = sg.render(ps, *_tracks(calls), {}, SimpleNamespace(), sg.style_for())
    assert any("or to the tandem copy of the gene it was found from" in t
               for t in _numbers(res["legend_groups"]))
    # the number in the figure is the one in the GFF, as before
    assert ">93<" in svg
    full = ps.render_anchor_grid(*_tracks(calls), {}, {}, SimpleNamespace())
    assert "identity 93%, coverage 100% (to GOI_copy_2, a seeded tandem copy, not to the query)" in full


def test_models_of_the_query_say_nothing_more():
    calls = [_goi(7000), _goi(9000, conf="MEDIUM", ident=55.0)]
    svg, legend, source, info = ps.render_anchor_grid_print(*_tracks(calls), {}, SimpleNamespace())
    assert info["copy_lineage_models"] == 0
    assert {r["identity_to"] for r in source if r["row"] == 2 and r["column"] == "T"} == {"query"}
    assert "tandem copy" not in legend
    res = sg.render(ps, *_tracks(calls), {}, SimpleNamespace(), sg.style_for())
    assert _numbers(res["legend_groups"])[0] == ("after a GOI model: % identity to the query; "
                                                 "bold = high confidence")


def test_every_source_row_names_what_its_identity_is_to():
    calls = [_goi(7000, parent="GOI_copy_1")]
    tracks, colours = _tracks(calls)
    tracks.append(dict(tracks[1], genome_id="g1", label="T1", species="Bombus species1",
                       genes=[dict(g) for g in tracks[1]["genes"]]))
    tracks[1]["genes"] = [g for g in tracks[1]["genes"] if g.get("home_gene_id") != "gene-F5"]
    _svg, _legend, source, _info = ps.render_anchor_grid_print(tracks, colours, {}, SimpleNamespace())
    assert set(source[0]) == set(sg.SOURCE_COLUMNS) and "identity_to" in sg.SOURCE_COLUMNS_CLADES
    by = {(r["row"], r["column"]): r["identity_to"] for r in source}
    assert by[(2, "F0")] == "home gene"                  # a flanking ortholog
    assert by[(2, "F5")] == ""                           # not placed
    assert by[(1, "F0")] == "" and by[(1, "T")] == ""    # the home row
    assert by[(2, "T")] == "GOI_copy_1"
    assert sg.copy_lineage({"home_gene_id": "GOI_copy_12"}) == "GOI_copy_12"
    for parent in ("GOI_T", "GOI_T~rel", "gene-F1", "", None):
        assert sg.copy_lineage({"home_gene_id": parent}) == ""


def test_a_copy_lineage_model_drawn_without_a_number_needs_no_sentence():
    """Only the numbers drawn after a model are explained: a copy drawn as a small arrow carries
    its copy number, not an identity."""
    calls = [_goi(7000), _goi(9000, parent="GOI_copy_2", ident=20.0)]    # too low for a number
    _svg, legend, source, info = ps.render_anchor_grid_print(*_tracks(calls), {}, SimpleNamespace())
    assert info["copy_lineage_models"] == 0 and "tandem copy of T" not in legend
    assert [r["identity_to"] for r in source if r["row"] == 2 and r["column"] == "T"] == [
        "query", "GOI_copy_2"]
