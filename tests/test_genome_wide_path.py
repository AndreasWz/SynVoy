"""The genome-wide path for strong GOI hits outside every searched window (FINAL step 5, batch
of 2026-10-15; --goi_genome_wide_path, off by default): places of GOI hits outside the windows,
modelled when significant over the genome (E <= 1e-3) and specific, at most the home family
size + 1 per GOI and genome, best E first; windows judged in the genome place. All synthetic."""
import os
import sys
import types

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import generate_report as gr  # noqa: E402
import iterative_search_runner as isr  # noqa: E402
import synvoy_evalue  # noqa: E402


def _hit(q, chrom, s, e, strand="+", pident=30.0):
    return {"query": q, "chrom": chrom, "start": s, "end": e, "strand": strand, "pident": pident,
            "qstart": 1, "qend": (e - s) // 3, "evalue": 1e-5, "bits": 40.0}


def test_places_are_the_goi_hits_outside_every_window():
    hits = [_hit("GOI_q", "c1", 1000, 1300), _hit("GOI_q|exon_2", "c1", 4000, 4300),   # one place
            _hit("GOI_q", "c1", 90000, 90300),                                         # a second
            _hit("GOI_q|Sp_fna_b0_l1_exon_ann", "c2", 500, 800),                        # a model hit
            _hit("GOI_q", "c1", 200000, 200300),                                       # in a window
            _hit("gene-FLANK", "c1", 5000, 5300),                                      # flanking
            _hit("GOI_q~rev", "c3", 100, 400)]                                         # a decoy
    places = isr.genome_wide_places(hits, [("c1", 150000, 250000)], {"GOI_q"}, max_gap=10000)
    assert [(p["chrom"], p["start"], p["end"], p["hits"]) for p in places] == \
        [("c1", 1000, 4300, 2), ("c1", 90000, 90300, 1), ("c2", 500, 800, 1)]
    assert {p["parent"] for p in places} == {"GOI_q"}


def test_hits_of_another_gene_need_their_own_home_parent():
    hits = [_hit("GOI_copy_1|Sp_fna_b0_l1", "c1", 1000, 1300), _hit("GOI_a", "c1", 5000, 5300)]
    assert [p["parent"] for p in isr.genome_wide_places(hits, [], {"GOI_q"}, 100)] == ["GOI_q", "GOI_q"]
    two = isr.genome_wide_places(hits, [], {"GOI_a", "GOI_b"}, 100)
    assert [(p["parent"], p["start"]) for p in two] == [("GOI_a", 5000)]       # copy: no parent


def _place(parent, chrom, start, e, specific=True):
    return {"parent": parent, "chrom": chrom, "start": start, "end": start + 300, "evalue": e,
            "specific": specific}


def test_choice_significant_specific_best_e_first_up_to_the_cap():
    places = [_place("GOI_q", "c1", 100, 2e-10), _place("GOI_q", "c2", 50, 1e-5),
              _place("GOI_q", "c1", 900, 1e-5), _place("GOI_q", "c3", 0, 0.5),
              _place("GOI_q", "c4", 0, 1e-20, specific=False)]
    kept = isr.choose_genome_wide_places(places, lambda parent: 2)
    assert [(p["chrom"], p["start"]) for p in kept] == [("c1", 100), ("c1", 900)]   # tie: position
    assert [p["decision"] for p in places] == ["modelled", "cut_by_cap", "modelled",
                                               "not_significant", "not_specific"]
    assert [p["rank"] for p in places] == [1, 3, 2, "", ""] and [p["cap"] for p in places] == [2, 2, 2, "", ""]
    asked = []
    isr.choose_genome_wide_places([_place("GOI_q", "c1", 0, 0.5)], lambda parent: asked.append(parent) or 1)
    assert asked == []                                                  # no eligible place: no count


def test_the_cap_is_per_home_gene():
    places = [_place("GOI_a", "c1", 0, 1e-9), _place("GOI_a", "c2", 0, 1e-8), _place("GOI_b", "c3", 0, 1e-7)]
    kept = isr.choose_genome_wide_places(places, lambda parent: 1)
    assert [(p["parent"], p["chrom"]) for p in kept] == [("GOI_a", "c1"), ("GOI_b", "c3")]


def test_the_window_stays_outside_every_other_window():
    place = {"chrom": "c1", "start": 10000, "end": 11000}
    assert isr.genome_wide_window(place, 5000, 100000, []) == (5000, 16000)
    assert isr.genome_wide_window(place, 5000, 12000, []) == (5000, 12000)               # chromosome end
    wins = [("c1", 0, 8000), ("c1", 13000, 20000), ("c2", 0, 100000)]
    assert isr.genome_wide_window(place, 5000, 100000, wins) == (8000, 13000)


class FakeProteome:
    def __init__(self, scores, family, own="gene-q"):
        self._scores, self._family, self._own = scores, family, own
        self.records = [(rid, "M" * length) for rid, length in
                        {"gene-q": 100, "gene-a": 100, "gene-b": 100, "gene-c": 100}.items()]
        self.tested = []

    def own_gene(self, q_id, q_seq, gene_id=None):
        return gene_id or self._own

    def scores(self, seq):
        return self._scores

    def family_evalue(self, q_id, q_seq, p_id):
        self.tested.append(p_id)
        return self._family[p_id]


def test_home_family_size_tests_only_what_the_prefilter_lets_through():
    isr._FAMILY_SIZE_CACHE.clear()
    # raw SW scores: q itself, two relatives, one chance protein (KA E > 0.1 for 100 x 100)
    proteome = FakeProteome([(900.0, "gene-q"), (120.0, "gene-a"), (80.0, "gene-b"), (30.0, "gene-c")],
                            {"gene-q": 1e-90, "gene-a": 1e-20, "gene-b": 0.05, "gene-c": 1e-9})
    assert isr.home_family_size(proteome, "GOI_q", "M" * 100) == 1
    assert proteome.tested == ["gene-a", "gene-b"]                       # q and gene-c not tested
    proteome.tested.clear()
    assert isr.home_family_size(proteome, "GOI_q", "M" * 100, "gene-a") == 1   # own gene by id
    assert proteome.tested == ["gene-q", "gene-b"]


def test_genome_wide_windows_are_no_blocks():
    ctx = {"genome_len": 10_000_000, "places": {}, "block_esyn": {0: 1e-5}}
    args = types.SimpleNamespace(genome_ctx=ctx, home_rank={}, home_geometry={}, expansion=1.0)
    terms = isr.block_place_terms(args, 7, 100, 9000, 100, 9000, 3, [], genome_wide=True)
    assert terms["significant"] is False and terms["slot"] is None and ctx["block_esyn"] == {0: 1e-5}
    assert synvoy_evalue.place_of((200, 400), terms["significant"], terms["slot"], terms["core"],
                                  terms["window"], terms["genome_len"])[0] == "genome"


def test_hit_pieces_are_split_at_stops_on_both_strands():
    seq = "ATGGCT" + "TAA" + "GCTGCTGCTGCTGCT"
    assert isr._hit_pieces({"c": seq}, {"chrom": "c", "start": 0, "end": 24, "strand": "+"}) == ["AAAAA"]
    rc = isr.reverse_complement(seq)
    assert isr._hit_pieces({"c": rc}, {"chrom": "c", "start": 0, "end": 24, "strand": "-"}) == ["AAAAA"]


def test_the_pass_needs_a_significant_block_and_the_home_proteome(monkeypatch):
    args = types.SimpleNamespace(genome_ctx={"block_esyn": {0: 0.5}}, home_goi_ids={"GOI_q"})
    assert isr.genome_wide_goi_pass(args, [], {}, {}, "g", "u", 1, None, "", 5000) == ([], [])
    args.genome_ctx["block_esyn"][1] = 1e-4
    monkeypatch.setattr(isr, "_home_proteome", lambda a: (None, None))
    assert isr.genome_wide_goi_pass(args, [], {}, {}, "g", "u", 1, None, "", 5000) == ([], [])


def test_the_report_lists_the_places_and_warns_only_when_one_was_modelled(tmp_path):
    path = tmp_path / "locus_1_G.fna.genome_wide.tsv"
    path.write_text("\t".join(isr.GENOME_WIDE_COLUMNS) + "\n"
                    + "G.fna\tGOI_q\tc1\t100\t400\t2\t27.5\t2e-10\t60.1\tgene-p\t17.0\t1\t1\t1\t2\tmodelled\n"
                    + "G.fna\tGOI_q\tc2\t100\t400\t1\t40.0\t1e-05\t50.0\tgene-p\t12.0\t1\t1\t3\t2\tcut_by_cap\n")
    records = gr.collect_genome_wide_places([str(path)])
    block = gr.summarize_genome_wide_places(records)
    assert block["total"] == 2 and block["by_decision"] == {"modelled": 1, "cut_by_cap": 1}
    assert [r["evalue_genome"] for r in block["places"]] == ["2e-10", "1e-05"]
    assert "no decoy calibration" in block["note"]
