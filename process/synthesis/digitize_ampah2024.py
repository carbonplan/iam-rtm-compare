"""ERW bars of Ampah et al. (2024), Nat. Commun. 15, 6342, read from the vector figures in the article PDF.

The Zenodo deposit (10.5281/zenodo.11254051) holds only GCAM-TJU inputs, so the figures are the only outputs:
  Fig. 3a, HIGH panel: Asia's gross CO2 removal by type, 2025-2050 (Gt CO2/yr)
  Fig. 4c: removal by type and country/region in 2050, HIGH scenario (Mt CO2/yr)
Only HIGH deploys ERW. Bar heights come from the PDF drawing commands, scaled by the y tick labels.
Needs pdfplumber (not in the repo env). Usage: python digitize_ampah2024.py s41467-024-50594-5.pdf out.csv
"""
import csv
import sys

import numpy as np
import pdfplumber

ERW = (0.93333, 0.77647, 0.69412)  # fill colour of the ERW legend swatch
REGIONS = ("Central Asia", "China", "India", "Indonesia", "Japan", "Middle East", "Pakistan", "South Asia",
           "South Korea", "Southeast Asia", "Taiwan")  # Fig. 4c x axis, left to right ("Taiwan (China)")
YEARS = range(2025, 2051, 5)


def yscale(page, ticks, x_max, top, bottom):
    """Value per point from the panel's y tick labels (left of x_max), fit through the label centres."""
    pos = {}
    for w in page.extract_words():
        if w["text"] in ticks and w["x1"] < x_max and top < w["top"] < bottom:
            assert w["text"] not in pos, w
            pos[w["text"]] = (w["top"] + w["bottom"]) / 2
    assert set(pos) == set(ticks), pos
    slope = np.polyfit([pos[t] for t in ticks], [float(t) for t in ticks], 1)[0]
    return -slope


def erw_rects(page, x0, x1, top, bottom):
    return [r for r in page.rects
            if r["non_stroking_color"] == ERW and x0 < r["x0"] < x1 and top < r["top"] < bottom]


def bar_centres(page, x0, x1, top, bottom, n):
    """x centres of the n stacked bars in the panel (any fill colour except greys and white)."""
    xs = sorted({round((r["x0"] + r["x1"]) / 2, 1) for r in page.rects
                 if x0 < r["x0"] < x1 and top < r["top"] < bottom and 2 < r["width"] < 12
                 and len(set(r["non_stroking_color"] or ())) > 1})
    groups = np.split(np.array(xs), np.where(np.diff(xs) > 3)[0] + 1)
    assert len(groups) == n, xs
    return [g.mean() for g in groups]


def nearest(x, centres):
    i = int(np.argmin([abs(x - c) for c in centres]))
    assert abs(x - centres[i]) < 3, (x, centres)
    return i


pdf = pdfplumber.open(sys.argv[1])
rows = []

# Fig. 3a, HIGH panel (page 5): Asia total by year
page = pdf.pages[4]
box = (125, 205, 70, 260)                                    # x0, x1, top, bottom of the HIGH panel
per_pt = yscale(page, ("30", "15", "0", "-15"), 125, *box[2:])         # Gt per point
centres = bar_centres(page, *box, n=len(YEARS))
asia = dict.fromkeys(YEARS, 0.0)
for r in erw_rects(page, *box):
    asia[YEARS[nearest((r["x0"] + r["x1"]) / 2, centres)]] += r["height"] * per_pt
rows += [("Fig. 3a", "HIGH", "Asia", y, v) for y, v in asia.items()]

# Fig. 4c (page 6): 2050 by country/region
page = pdf.pages[5]
box = (120, 250, 285, 415)
per_pt = yscale(page, ("6000", "4000", "2000", "0"), 125, *box[2:]) / 1000   # Gt per point
centres = bar_centres(page, *box, n=len(REGIONS))
region = dict.fromkeys(REGIONS, 0.0)
for r in erw_rects(page, *box):
    region[REGIONS[nearest((r["x0"] + r["x1"]) / 2, centres)]] += r["height"] * per_pt
rows += [("Fig. 4c", "HIGH", k, 2050, v) for k, v in region.items()]

with open(sys.argv[2], "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["figure", "scenario", "region", "year", "erw_cdr_GtCO2_per_yr"])
    w.writerows((fig, s, reg, y, "%.4g" % v) for fig, s, reg, y, v in rows)
print("Asia (Fig. 3a):", {y: round(float(v), 3) for y, v in asia.items()})
print("2050 by region (Fig. 4c):", {k: round(float(v), 3) for k, v in region.items()}, "sum %.3f" % sum(region.values()))
