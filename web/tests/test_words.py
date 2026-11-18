"""Sentences of the result page (site/assets/js/words.js), run with Node: a strong hit
without synteny (an answer of the genome-wide path) is named as the report and the figures
name it, listed apart from the calls set aside as another gene, and never a weak match."""

import json
import shutil
import subprocess

import pytest

from conftest import WEB_DIR

WORDS = WEB_DIR / "site" / "assets" / "js" / "words.js"
RENDER = WEB_DIR / "site" / "assets" / "js" / "render.js"


def _node():
    node = shutil.which("node")
    if not node:
        return None
    try:
        version = subprocess.run([node, "--version"], capture_output=True, text=True,
                                 timeout=20).stdout.strip().lstrip("v")
        return node if int(version.split(".")[0]) >= 18 else None
    except (OSError, ValueError, subprocess.SubprocessError):
        return None


NODE = _node()
pytestmark = pytest.mark.skipif(NODE is None, reason="needs Node 18 or newer")


def say(*calls):
    """[(name, argument), ...] -> what words.js returns for each (a constant: its value)."""
    script = ("const w = await import(process.argv[1]);"
              "const calls = JSON.parse(process.argv[2]);"
              "console.log(JSON.stringify(calls.map(([f, a]) => "
              "typeof w[f] === 'function' ? w[f](a) : w[f])));")
    done = subprocess.run([NODE, "--input-type=module", "-e", script, WORDS.as_uri(),
                           json.dumps([list(c) for c in calls])],
                          capture_output=True, text=True, timeout=60)
    assert done.returncode == 0, done.stderr
    return json.loads(done.stdout)


STRONG = {"label": "set_aside", "reason": "strong_hit", "identity": 49.3, "coverage": 0.772,
          "chrom": "chr9", "start": 20001, "end": 21500}
NEIGHBOURS = {"found": 4, "together": 4, "places": 1, "total": 10, "locus": 1}


def test_a_genome_with_a_strong_hit_names_it_and_not_as_a_weak_match():
    alone, beside = say(
        ("notFoundText", {"neighbours": NEIGHBOURS, "leads": None, "calls": [STRONG]}),
        ("notFoundText", {"neighbours": NEIGHBOURS, "calls": [STRONG],
                          "leads": {"n": 2, "identity": 21.5, "coverage": 0.43}}))
    assert alone == ("Part of its neighbourhood is there (4 genes side by side); the gene was not "
                     "found in it. A strong hit without synteny, not counted, is listed under the table.")
    assert "weak" not in alone
    assert beside.startswith(alone.split(" A strong hit")[0] + " 2 weak matches of the gene")
    assert beside.endswith("not enough for a call. A strong hit without synteny, not counted, is "
                           "listed under the table.")
    two, nothing = say(("notFoundText", {"calls": [STRONG, dict(STRONG, start=1)]}),
                       ("notFoundText", {"neighbours": NEIGHBOURS, "calls": []}))
    assert two == "2 strong hits without synteny, not counted, are listed under the table."
    assert "strong" not in nothing


def test_strong_hits_have_their_own_list_under_the_table():
    out = say(("setAsideHeading", 1), ("setAsideHeading", 2), ("strongHitHeading", 1),
              ("strongHitHeading", 3), ("setAsideReason", STRONG),
              ("setAsideReason", {"reason": "related_gene", "related_gene": "BGN"}),
              ("setAsideReason", {"reason": "short_match"}), ("STRONG_HIT_MEANS", None))
    assert out[:4] == ["Set aside: 1 call that is not the gene", "Set aside: 2 calls that are not the gene",
                       "Not counted: 1 strong hit without synteny",
                       "Not counted: 3 strong hits without synteny"]
    assert out[4:7] == ["a strong match where no neighbourhood places the gene",
                        "matches the related gene BGN better than the query",
                        "matches only a short piece of the query"]
    assert "not counted" in out[7]
    render = RENDER.read_text()
    assert "words.strongHitHeading(strong.length)" in render
    assert "words.setAsideHeading(aside.length)" in render
    assert "call that is not the gene" not in render         # the heading's words live in words.js


def test_the_headline_and_the_details_name_the_strong_hits():
    summary = {"genomes": 3, "not_searched": 0, "found": 1, "clear": 1, "candidate": 0, "none": 2,
               "candidates": 1, "set_aside": 3, "strong_hits": 2}
    older = {k: v for k, v in summary.items() if k != "strong_hits"}     # a bundle made before
    with_hits, without, fact, old_fact = say(
        ("headline", {"summary": summary, "query": {"gene": "GEN1"}}),
        ("headline", {"summary": older, "query": {"gene": "GEN1"}}),
        ("asideFact", summary), ("asideFact", older))
    assert with_hits == without + " 2 strong hits without synteny, not counted, are listed under the table."
    assert "strong" not in without
    assert (fact, old_fact) == ("1 ambiguous, 3 set aside (2 strong hits without synteny among them)",
                                "1 ambiguous, 3 set aside")
