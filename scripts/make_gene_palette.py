#!/usr/bin/env python3
"""Derive the flanking-gene palette of the figures (plot_synteny.GENE_PALETTE).

Pale, distinct colours, so that the GOI (red, or the clade colours) stands out. The core is
Paul Tol's "light" scheme (Tol 2021, technical note SRON/EPS/TN/09-002 issue 3.2: "reasonably
distinct for both normal or colour-blind vision", made to fill labelled cells), without its
grey (taken by "no ortholog") and its salmon (CIEDE2000 20 from the GOI red). It is extended
from the same pale band (CIELAB lightness 66-86, chroma 20-38): the next colour is always the
candidate farthest (CIEDE2000) from every colour chosen so far and from the reserved ones (GOI
red and border, white, the figure greys). The seven Tol colours come first, in the same
farthest-first order, so that genes that sit next to each other contrast.

    python3 scripts/make_gene_palette.py [N]      # prints the palette and its distances
"""
import math
import sys

import numpy as np

TOL_LIGHT = ("#77aadd", "#99ddff", "#44bb99", "#bbcc33", "#aaaa00", "#eedd88", "#ffaabb")
RESERVED = ("#e31a1c", "#8b0000", "#ffffff", "#d9d9d9", "#9aa3b2", "#1a1a1a")

_M = np.array([[0.4124564, 0.3575761, 0.1804375], [0.2126729, 0.7151522, 0.0721750],
               [0.0193339, 0.1191920, 0.9503041]])
_WP = np.array([0.95047, 1.0, 1.08883])


def _f(t):
    d = 6 / 29
    return np.where(t > d ** 3, np.cbrt(t), t / (3 * d * d) + 4 / 29)


def _finv(t):
    d = 6 / 29
    return np.where(t > d, t ** 3, 3 * d * d * (t - 4 / 29))


def hex_to_lab(h):
    rgb = np.array([int(h[i:i + 2], 16) for i in (1, 3, 5)], float) / 255.0
    lin = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    fx, fy, fz = _f(_M @ lin / _WP)
    return (116 * fy - 16, 500 * (fx - fy), 200 * (fy - fz))


def _lab_to_rgb(lab):
    L, a, b = lab
    fy = (L + 16) / 116
    lin = np.linalg.inv(_M) @ (np.array([_finv(fy + a / 500), _finv(fy), _finv(fy - b / 200)]) * _WP)
    return np.where(lin <= 0.0031308, 12.92 * lin,
                    1.055 * np.power(np.clip(lin, 0, None), 1 / 2.4) - 0.055)


def de2000(l1, l2):
    """CIEDE2000 colour difference of two CIELAB colours."""
    L1, a1, b1 = l1
    L2, a2, b2 = l2
    Cb = (math.hypot(a1, b1) + math.hypot(a2, b2)) / 2
    G = 0.5 * (1 - math.sqrt(Cb ** 7 / (Cb ** 7 + 25 ** 7)))
    a1p, a2p = (1 + G) * a1, (1 + G) * a2
    C1p, C2p = math.hypot(a1p, b1), math.hypot(a2p, b2)
    h1p, h2p = math.degrees(math.atan2(b1, a1p)) % 360, math.degrees(math.atan2(b2, a2p)) % 360
    dh = h2p - h1p
    if C1p * C2p == 0:
        dh = 0
    elif dh > 180:
        dh -= 360
    elif dh < -180:
        dh += 360
    dHp = 2 * math.sqrt(C1p * C2p) * math.sin(math.radians(dh / 2))
    Lbp, Cbp = (L1 + L2) / 2, (C1p + C2p) / 2
    if C1p * C2p == 0:
        hbp = h1p + h2p
    elif abs(h1p - h2p) <= 180:
        hbp = (h1p + h2p) / 2
    else:
        hbp = (h1p + h2p + (360 if h1p + h2p < 360 else -360)) / 2
    T = (1 - 0.17 * math.cos(math.radians(hbp - 30)) + 0.24 * math.cos(math.radians(2 * hbp))
         + 0.32 * math.cos(math.radians(3 * hbp + 6)) - 0.20 * math.cos(math.radians(4 * hbp - 63)))
    Rc = 2 * math.sqrt(Cbp ** 7 / (Cbp ** 7 + 25 ** 7))
    Rt = -math.sin(math.radians(2 * 30 * math.exp(-(((hbp - 275) / 25) ** 2)))) * Rc
    Sl = 1 + 0.015 * (Lbp - 50) ** 2 / math.sqrt(20 + (Lbp - 50) ** 2)
    Sc, Sh = 1 + 0.045 * Cbp, 1 + 0.015 * Cbp * T
    dL, dC = (L2 - L1) / Sl, (C2p - C1p) / Sc
    return math.sqrt(dL ** 2 + dC ** 2 + (dHp / Sh) ** 2 + Rt * dC * (dHp / Sh))


def make_palette(n=30):
    red, border = hex_to_lab(RESERVED[0]), hex_to_lab(RESERVED[1])
    cands = []
    for L in (66, 70, 74, 78, 82, 86):
        for C in (20, 26, 32, 38):
            for h in range(0, 360, 4):
                rgb = _lab_to_rgb((L, C * math.cos(math.radians(h)), C * math.sin(math.radians(h))))
                if np.all(rgb >= 0) and np.all(rgb <= 1):
                    hx = "#%02x%02x%02x" % tuple(int(round(v * 255)) for v in rgb)
                    lab = hex_to_lab(hx)
                    if de2000(lab, red) > 30 and de2000(lab, border) > 22:
                        cands.append((hx, lab))
    seeds = [(h, hex_to_lab(h)) for h in TOL_LIGHT]
    chosen = [seeds[0][0]]
    taken = [hex_to_lab(h) for h in RESERVED] + [seeds[0][1]]
    pool = seeds[1:]
    while len(chosen) < n:
        if not pool:                    # the Tol colours are placed; now the extension
            pool = cands
        best = max(pool, key=lambda c: min(de2000(c[1], t) for t in taken))
        chosen.append(best[0])
        taken.append(best[1])
        pool = [c for c in pool if c[0] != best[0]]
    return chosen


if __name__ == "__main__":
    pal = make_palette(int(sys.argv[1]) if len(sys.argv) > 1 else 30)
    labs = [hex_to_lab(h) for h in pal]
    print(pal)
    print("smallest CIEDE2000 between two colours:",
          round(min(de2000(a, b) for i, a in enumerate(labs) for b in labs[i + 1:]), 1))
