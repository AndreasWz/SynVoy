"""Figure legends that describe what the drawing does.

- The tree page said its clades were cut "at its longest branches". The clades come from
  synvoy_tree.partition_clades, which splits the largest clade from the top and never reads a
  branch length.
- The anchor grid said "paler = lower identity" for every flanking arrow. Only a
  high-confidence ortholog is shaded by its identity (cell_style); a medium- or low-confidence
  one has a fixed pale fill whatever its identity, so that a weak call never looks confident.
"""
import os
import re
import sys
from types import SimpleNamespace

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import plot_synteny as ps  # noqa: E402
import synvoy_grid as sg  # noqa: E402
import synvoy_tree as stree  # noqa: E402


def test_the_tree_page_names_the_rule_that_makes_its_clades(tmp_path):
    nwk = tmp_path / "locus_1_tree.nwk"
    nwk.write_text("((GOI_X|A_fna_b0_l1:0.1,GOI_X|B_fna_b0_l1:0.1):0.2,"
                   "(GOI_X|C_fna_b0_l1:0.1,GOI_X|D_fna_b0_l1:0.1):0.2,GOI_X|E_fna_b0_l1:5.0);\n")
    out = tmp_path / "synteny_block_locus_1_tree.html"
    ps._render_tree_svg(str(nwk), {}, str(out), species_map={}, clade_count=3, gene_legend=False)
    html = out.read_text()
    legend = " ".join(re.findall(r"<text[^>]*>([^<]*)</text>", html.split('class="figure-legend"', 1)[1]))
    assert ("split from the root into 3 clades, the largest clade first (by topology, not branch "
            "length)") in legend
    assert "longest branches" not in html and "longest branches" not in legend
    # the rule the page names: the leaf on the 5.0 branch is not cut off first
    t = stree.midpoint_root(stree.parse_newick_tree(nwk.read_text()))
    groups = {}
    for leaf, cid in stree.partition_clades(t, target_k=3).items():
        groups.setdefault(cid, set()).add(leaf.split("|")[1][0])
    assert sorted(map(sorted, groups.values())) == [["A", "B"], ["C", "D"], ["E"]]
    src = open(os.path.join(os.path.dirname(__file__), "..", "bin", "plot_synteny.py")).read()
    assert "longest non-root branches" not in src and "K-longest-branches" not in src


def _tracks():
    names = [f"gene-F{k}" for k in range(4)]

    def g(name, start, **kw):
        d = {"chrom": "chr1", "start": start, "end": start + 900, "start_plot": start,
             "end_plot": start + 900, "name": name, "strand": "+", "exon_coords": []}
        d.update(kw)
        return d
    home = [g(n, 1000 + 2000 * k) for k, n in enumerate(names)]
    home.insert(2, g("GOI_T", 4600))
    target = [g(f"t_{k}", 1000 + 2000 * k, home_gene_id=n, identity=95.0 - 30 * (k % 2),
                confidence=("HIGH", "MEDIUM", "HIGH", "MEDIUM")[k], query_coverage=0.95)
              for k, n in enumerate(names)]
    target.append(g("GOI_T|m1", 4700, role="goi", home_gene_id="GOI_T", confidence="HIGH",
                    identity=70.0, evidence_type="exon_annotation",
                    exon_coords=[(4700, 4900), (5200, 5600)]))
    tracks = [{"label": "Home", "species": "Apis mellifera", "is_home": True, "genome_id": "home",
               "offset": 0, "breaks": [], "genes": home},
              {"label": "T0", "species": "Bombus species0", "is_home": False, "genome_id": "g0",
               "offset": 0, "breaks": [], "genes": target}]
    return tracks, {n: ps.GENE_PALETTE[k] for k, n in enumerate(names)}


@pytest.mark.parametrize("style", [sg.style_for, sg.clade_style, sg.clade_muted_style])
def test_the_paleness_legend_names_the_tier_it_applies_to(style):
    res = sg.render(ps, *_tracks(), {}, SimpleNamespace(), style())
    flank = dict(res["legend_groups"])["Flanking genes"]
    assert [t for _, t in flank if "paler" in t] == ["high confidence: paler = lower identity"]
    text = res["legend_text"]
    assert ("a high-confidence ortholog is paler the lower its protein identity to the home gene; "
            "a medium- or low-confidence one has a fixed pale fill.") in text
    assert "paler fill = lower protein identity" not in text


def test_with_numbers_the_fill_says_the_confidence():
    res = sg.render(ps, *_tracks(), {}, SimpleNamespace(), sg.style_for(numbers=True))
    assert not any("paler" in t for _, rows in res["legend_groups"] for _, t in rows)
    # AMBIGUOUS is drawn paler than LOW though it ranks above it, so the text
    # says what cell_style draws, not "paler the lower the confidence".
    assert ("Colour: the home gene the ortholog belongs to; the fill follows the label, high "
            "confidence darkest, then medium and low, ambiguous palest. Number: protein identity "
            "to the home gene in %") in res["legend_text"]
    assert "paler the lower the confidence" not in res["legend_text"]


def test_the_fill_the_legend_describes(monkeypatch):
    monkeypatch.setattr(sg, "_PS", ps)
    base = ps.GENE_PALETTE[0]
    fill = lambda ident, conf, numbers=False: sg.cell_style(base, ident, conf, False, numbers)[0]  # noqa: E731
    assert fill(95, "HIGH") != fill(40, "HIGH")                       # shaded by identity
    for conf in ("MEDIUM", "LOW", "AMBIGUOUS"):                        # a fixed pale fill
        assert fill(95, conf) == fill(40, conf)
    assert fill(95, "HIGH", True) == fill(40, "HIGH", True)           # numbers: by confidence only
    assert len({fill(80, c, True) for c in ("HIGH", "MEDIUM", "LOW", "AMBIGUOUS")}) == 4
    # the order the legend text states: HIGH darkest, then MEDIUM, LOW, AMBIGUOUS palest
    lightness = [sum(int(fill(80, c, True)[i:i + 2], 16) for i in (1, 3, 5))
                 for c in ("HIGH", "MEDIUM", "LOW", "AMBIGUOUS")]
    assert lightness == sorted(lightness) and len(set(lightness)) == 4
