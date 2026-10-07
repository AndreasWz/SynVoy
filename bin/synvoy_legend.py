#!/usr/bin/env python3
"""Structured legend shared by every SynVoy figure (import-only helper).

A legend used to be a row of keys followed by one or two lines of text joined
with dots, which a reader has to scan from the start to find anything. Here a
legend is a set of titled groups laid out as columns:

    Layout                Flanking genes            GOI
    Columns  home order   [>] high confidence       [model] exons to scale
    Rows     genomes      [>] medium                [1] one copy, numbered
                          [>] low                   //  elsewhere in the genome

Each group answers one question (how is the figure laid out, what does a
flanking arrow mean, how is the GOI drawn, what do the numbers say), so the
place to look is found by its title.

    groups = [("Flanking genes", [(glyph, "high confidence"), ...]), ...]
    parts, height = legend_svg(groups, x, y, avail_w, text_width, fs=10)

A glyph is ``None`` (text only), a string (a key such as "×N", drawn bold) or
``(width, draw)`` with ``draw(x, y_mid) -> svg``, which paints the symbol with
its left edge at x and its vertical centre at y_mid.
"""
from __future__ import annotations

from html import escape
from typing import Callable, List, Optional, Sequence, Tuple, Union

Glyph = Union[None, str, Tuple[float, Callable[[float, float], str]]]
Entry = Tuple[Glyph, str]
Group = Tuple[str, Sequence[Entry]]

INK = "#1a1d26"
MUTED = "#4b5563"
RULE = "#d8dee8"


def _wrap(text: str, max_w: float, measure) -> List[str]:
    """Greedy word wrap of *text* into lines no wider than *max_w*."""
    lines, cur = [], ""
    for word in text.split():
        trial = f"{cur} {word}".strip()
        if cur and measure(trial) > max_w:
            lines.append(cur)
            cur = word
        else:
            cur = trial
    if cur:
        lines.append(cur)
    return lines or [""]


def glyph_width(glyph: Glyph, fs: float, text_width) -> float:
    if glyph is None:
        return 0.0
    if isinstance(glyph, str):
        return text_width(glyph, fs, True)
    return float(glyph[0])


def _boxes(groups, avail_w, text_width, fs, max_text_w):
    """(title, key_w, text_x, rows, width, height) per non-empty group."""
    lh, head_h, key_gap = fs * 1.5, fs * 1.9, fs * 0.55
    max_text_w = max_text_w or fs * 36
    boxes = []
    for title, entries in groups:
        entries = [e for e in entries if e]
        if not entries:
            continue
        key_w = max(glyph_width(g, fs, text_width) for g, _ in entries)
        text_x = key_w + (key_gap if key_w else 0.0)
        limit = max(fs * 8, min(max_text_w, avail_w - text_x))
        rows = [(g, _wrap(t, limit, lambda s: text_width(s, fs))) for g, t in entries]
        text_w = max(text_width(line, fs) for _, lines in rows for line in lines)
        width = max(text_x + text_w, text_width(title, fs, True))
        height = head_h + sum(len(lines) for _, lines in rows) * lh
        boxes.append((title, key_w, text_x, rows, width, height))
    return boxes


def natural_width(groups: Sequence[Group], text_width, fs: float = 10.0,
                  max_text_w: Optional[float] = None) -> float:
    """Width the legend takes with all its groups side by side."""
    boxes = _boxes(groups, float("inf"), text_width, fs, max_text_w)
    return sum(b[4] for b in boxes) + fs * 2.6 * max(0, len(boxes) - 1)


def legend_svg(groups: Sequence[Group], x: float, y: float, avail_w: float, text_width,
               fs: float = 10.0, ink: str = INK, muted: str = MUTED, rule: Optional[str] = RULE,
               max_text_w: Optional[float] = None, font_family: str = "") -> Tuple[List[str], float]:
    """SVG of a grouped legend with its top-left corner at (x, y).

    Groups are set side by side and continue on a new band when the width
    *avail_w* is used up. ``text_width(text, size, bold=False)`` measures text in
    the figure's font. Entry text wider than *max_text_w* (default 36 x the font
    size) wraps onto further lines. Returns ``(svg_parts, height)``; groups
    without entries are skipped, and no groups give height 0.
    """
    boxes = _boxes(groups, avail_w, text_width, fs, max_text_w)
    if not boxes:
        return [], 0.0
    lh = fs * 1.5                      # line pitch
    head_h = fs * 1.9                  # heading line + rule
    col_gap = fs * 2.6                 # between two groups
    band_gap = fs * 1.2                # between two bands of groups
    fam = f' font-family="{font_family}"' if font_family else ""

    parts: List[str] = []
    cx, cy, band_h = x, y, 0.0
    for title, key_w, text_x, rows, width, height in boxes:
        if cx > x and cx + width > x + avail_w + 0.01:
            cx, cy, band_h = x, cy + band_h + band_gap, 0.0
        parts.append(f'<text x="{cx:.2f}" y="{cy + fs:.2f}" font-size="{fs:g}" font-weight="700" '
                     f'fill="{ink}"{fam}>{escape(title)}</text>')
        if rule:
            parts.append(f'<line x1="{cx:.2f}" y1="{cy + fs * 1.42:.2f}" x2="{cx + width:.2f}" '
                         f'y2="{cy + fs * 1.42:.2f}" stroke="{rule}" stroke-width="{fs * 0.07:.2f}"/>')
        ly = cy + head_h
        for glyph, lines in rows:
            mid = ly + lh / 2
            base = mid + fs * 0.35
            if isinstance(glyph, str):
                parts.append(f'<text x="{cx:.2f}" y="{base:.2f}" font-size="{fs:g}" font-weight="700" '
                             f'fill="{ink}"{fam}>{escape(glyph)}</text>')
            elif glyph is not None:
                parts.append(glyph[1](cx + (key_w - glyph[0]) / 2, mid))
            for k, line in enumerate(lines):
                parts.append(f'<text x="{cx + text_x:.2f}" y="{base + k * lh:.2f}" font-size="{fs:g}" '
                             f'fill="{muted}"{fam}>{escape(line)}</text>')
            ly += len(lines) * lh
        band_h = max(band_h, height)
        cx += width + col_gap
    return parts, (cy + band_h) - y


def legend_text(groups: Sequence[Group]) -> str:
    """The legend as plain text, one group per line (for a caption draft)."""
    out = []
    for title, entries in groups:
        items = []
        for glyph, text in entries:
            items.append(f"{glyph} = {text}" if isinstance(glyph, str) else text)
        if items:
            out.append(f"{title}: " + "; ".join(items) + ".")
    return "\n".join(out)
