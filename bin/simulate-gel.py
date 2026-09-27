#!/usr/bin/env python3
"""
simulate-gel.py — Simulate a photorealistic agarose gel electrophoresis image (PNG)
from FASTQ read-length distributions.

Features:
  - Intercalating dye mass weighting (fluorescence proportional to base pairs).
  - Adaptive smoothing: sharp bands for concentrated amplicon peaks alongside silky,
    continuous smears for dispersed shearing/genomic libraries.
  - Authentic physical and optical artifacts:
      * Meniscus curvature ("smile" effect) within lanes.
      * Optical bloom / saturation on bright bands.
      * Agarose slab boundaries and transilluminator illumination vignette.
      * Comb wells with realistic pocket depth and loaded sample residual dye.
      * Diffusion-broadened molecular weight reference ladders.
  - Multiple visualization themes:
      * 'sybr': SYBR Green on dark background.
      * 'ethidium': Ethidium bromide (EtBr) fluorescence (warm orange-red).
      * 'chemidoc': Monochrome GelDoc / ChemiDoc transillumination.
  - Automatic multi-gel splitting & vertical stacking when samples exceed comb capacity.

Usage:
    python simulate-gel.py sample1.fastq.gz sample2.fastq.gz -o gel.png
    python simulate-gel.py *.fastq.gz --theme ethidium --ladder 1kb_plus
    python simulate-gel.py *.fastq.gz --lanes 8 -o multi_gel.png

Dependencies: numpy, matplotlib.
"""

import argparse
import base64
import gzip
import io
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches


# ---------------------------------------------------------------------------
# Ladders: name -> list of (size_bp, relative_intensity 0-1)
# Intensities approximate real reference-lane band brightness patterns.
# ---------------------------------------------------------------------------
LADDERS = {
    "1kb": [
        (10000, 0.55), (8000, 0.45), (6000, 0.45), (5000, 0.50),
        (4000, 0.55), (3000, 1.00), (2000, 0.60), (1500, 0.55),
        (1000, 0.85), (500, 0.65),
    ],
    "1kb_plus": [
        (12000, 0.4), (11000, 0.4), (10000, 0.4), (9000, 0.4),
        (8000, 0.4), (7000, 0.4), (6000, 0.4), (5000, 0.4),
        (4000, 0.4), (3000, 0.4), (2000, 1.0), (1650, 0.4),
        (1000, 0.9), (850, 0.4), (650, 0.4), (500, 0.85),
        (400, 0.4), (300, 0.4), (200, 0.4), (100, 0.4),
    ],
    "100bp": [
        (1500, 0.6), (1000, 0.6), (900, 0.5), (800, 0.5), (700, 0.5),
        (600, 0.5), (500, 1.0), (400, 0.5), (300, 0.5), (200, 0.5),
        (100, 0.7),
    ],
}

# ---------------------------------------------------------------------------
# Visual themes
# ---------------------------------------------------------------------------
THEMES = {
    "sybr": {
        "cmap": mcolors.LinearSegmentedColormap.from_list(
            "gel_sybr",
            [
                (0.00, "#040b06"),
                (0.06, "#081c10"),
                (0.18, "#0f451e"),
                (0.38, "#1fa33f"),
                (0.65, "#5ee673"),
                (0.85, "#bcf7c3"),
                (1.00, "#ffffff"),
            ]
        ),
        "bg_color": "#030805",
        "text_color": "#e0e8e2",
        "well_color": "#051008",
        "border_color": "#0d2b17",
        "accent_color": "#5ee673",
    },
    "ethidium": {
        "cmap": mcolors.LinearSegmentedColormap.from_list(
            "gel_ethidium",
            [
                (0.00, "#080406"),
                (0.06, "#14080c"),
                (0.18, "#420f18"),
                (0.38, "#b02e20"),
                (0.65, "#f07535"),
                (0.85, "#fad48c"),
                (1.00, "#ffffff"),
            ]
        ),
        "bg_color": "#060304",
        "text_color": "#ece2e4",
        "well_color": "#0e0508",
        "border_color": "#2a0d14",
        "accent_color": "#f07535",
    },
    "chemidoc": {
        "cmap": mcolors.LinearSegmentedColormap.from_list(
            "gel_chemidoc",
            [
                (0.00, "#060606"),
                (0.06, "#121212"),
                (0.20, "#333333"),
                (0.45, "#777777"),
                (0.70, "#bbbbbb"),
                (0.88, "#e6e6e6"),
                (1.00, "#ffffff"),
            ]
        ),
        "bg_color": "#050505",
        "text_color": "#e0e0e0",
        "well_color": "#0a0a0a",
        "border_color": "#1f1f1f",
        "accent_color": "#e6e6e6",
    },
}


def read_fastq_lengths(path):
    """Return a list of read lengths from a FASTQ file (gzip-aware)."""
    opener = gzip.open if path.endswith(".gz") else open
    lengths = []
    with opener(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:  # sequence line
                lengths.append(len(line.strip()))
    if not lengths:
        raise ValueError(f"No reads found in {path} (not a valid FASTQ?)")
    return lengths


def clean_sample_name(path):
    """Strip directory path and common FASTQ extensions from sample filename."""
    name = os.path.basename(path)
    for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(ext):
            name = name[:-len(ext)]
            break
    return name


def smooth1d(arr, sigma):
    """1D Gaussian convolution with normalized kernel."""
    if sigma <= 0:
        return arr
    radius = max(1, int(3.5 * sigma))
    x = np.arange(-radius, radius + 1)
    kernel = np.exp(-(x ** 2) / (2 * sigma ** 2))
    kernel /= kernel.sum()
    return np.convolve(arr, kernel, mode="same")


def bp_to_y(bp, min_bp, max_bp, y_top, y_bot):
    """Map a base-pair size to a vertical pixel coordinate on a log scale."""
    log_min, log_max = np.log10(min_bp), np.log10(max_bp)
    frac = (np.log10(max_bp) - np.log10(bp)) / (log_max - log_min)
    return np.clip(y_top + frac * (y_bot - y_top), y_top, y_bot)


def profile_from_lengths(lengths, min_bp, max_bp, y_top, y_bot, n_pixels,
                         band_sigma=3.2, smear_sigma=20.0):
    """Build an intensity profile from read lengths with adaptive smoothing:
    sharp bands for concentrated peaks, broad continuous smear for dispersed fragments.

    Fragments larger than max_bp are too big to migrate into the gel matrix, so
    rather than discarding them they're reported back as `well_frac`: the share
    of the sample's total dye-weighted mass that stays trapped in the loading
    well. The caller renders that as a retained band in the pocket itself.
    """
    lengths = np.asarray(lengths, dtype=float)
    total_mass = lengths.sum() if lengths.size else 0.0
    oversized_mass = lengths[lengths > max_bp].sum() if lengths.size else 0.0
    well_frac = (oversized_mass / total_mass) if total_mass > 0 else 0.0

    lengths = lengths[(lengths >= min_bp) & (lengths <= max_bp)]
    if len(lengths) == 0:
        return np.zeros(n_pixels), well_frac

    n_reads = len(lengths)
    ys = bp_to_y(lengths, min_bp, max_bp, y_top, y_bot)
    y_clipped = np.clip(ys.astype(int), 0, n_pixels - 1)
    mass = np.bincount(y_clipped, weights=lengths, minlength=n_pixels)[:n_pixels]
    counts = np.bincount(y_clipped, minlength=n_pixels)[:n_pixels].astype(float)

    # Pilot count density across gel
    pilot_count = smooth1d(counts, 6.0)
    frac_per_px = pilot_count / max(1, n_reads)

    # Adaptive sharp vs broad weighting
    # Peaks with high fraction of library reads render as sharp bands
    w_sharp = np.clip((frac_per_px - 0.006) / (0.025 - 0.006), 0.0, 1.0)

    sharp = smooth1d(mass, band_sigma)
    broad = smooth1d(mass, smear_sigma)
    profile = w_sharp * sharp + (1.0 - w_sharp) * broad

    if profile.max() > 0:
        profile = profile / profile.max()
        # Soft photographic gamma: preserves subtle library smear
        profile = profile ** 0.55
        # Conserve total dye signal between the well and the gel: only the
        # (1 - well_frac) share of mass actually migrated, so that's the most
        # the in-gel bands should ever read as, however bright they'd look
        # in isolation. Without this, a lane that's 90% oversized DNA would
        # still show a full-brightness band from its leftover 10%.
        profile = profile * (1.0 - well_frac)

    return profile, well_frac


def ladder_profile(ladder_bands, min_bp, max_bp, y_top, y_bot, n_pixels, ladder_sigma=2.6):
    """Build a molecular weight ladder profile with slight diffusion down the gel."""
    profile = np.zeros(n_pixels)
    for bp, intensity in ladder_bands:
        if bp < min_bp or bp > max_bp:
            continue
        y = bp_to_y(bp, min_bp, max_bp, y_top, y_bot)
        # Slight diffusion increase for smaller fragments farther down the gel
        band_w = ladder_sigma * (1.0 + 0.3 * (y - y_top) / max(1, y_bot - y_top))
        radius = int(round(3.5 * band_w))
        y_int = int(round(y))
        y_range = np.arange(max(0, y_int - radius), min(n_pixels, y_int + radius + 1))
        kernel = np.exp(-((y_range - y) ** 2) / (2 * band_w ** 2))
        profile[y_range] += intensity * kernel

    if profile.max() > 0:
        profile = profile / profile.max()
    return profile


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>__DOC_TITLE__</title>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
  :root {
    --bg-color: __BG_COLOR__;
    --text-color: __TEXT_COLOR__;
    --accent-color: __ACCENT_COLOR__;
    --border-color: __BORDER_COLOR__;
    --well-color: __WELL_COLOR__;
  }
  * {
    box-sizing: border-box;
    margin: 0;
    padding: 0;
  }
  body {
    background-color: var(--bg-color);
    color: var(--text-color);
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    min-height: 100vh;
    display: flex;
    flex-direction: column;
    align-items: center;
    padding: 24px 16px;
    overflow-x: auto;
  }
  header {
    width: 100%;
    max-width: 1100px;
    margin-bottom: 20px;
    display: flex;
    flex-wrap: wrap;
    justify-content: space-between;
    align-items: center;
    gap: 12px;
    padding-bottom: 12px;
    border-bottom: 1px solid rgba(255, 255, 255, 0.1);
  }
  .title-group h1 {
    font-size: 1.2rem;
    font-weight: 600;
    letter-spacing: -0.01em;
  }
  .title-group p {
    font-size: 0.72rem;
    color: rgba(255, 255, 255, 0.6);
    margin-top: 4px;
  }
  .controls {
    display: flex;
    align-items: center;
    gap: 10px;
  }
  .btn {
    background: var(--well-color);
    border: 1px solid var(--border-color);
    color: var(--text-color);
    padding: 6px 14px;
    font-size: 0.7rem;
    border-radius: 6px;
    cursor: pointer;
    transition: all 0.15s ease;
    display: inline-flex;
    align-items: center;
    gap: 6px;
  }
  .btn:hover {
    border-color: var(--accent-color);
    background: rgba(255, 255, 255, 0.08);
  }
  .badge-chip {
    background: rgba(255, 255, 255, 0.05);
    border: 1px solid var(--border-color);
    border-radius: 4px;
    padding: 2px 8px;
    font-size: 0.64rem;
    color: var(--accent-color);
  }
  .panels-container {
    display: flex;
    flex-direction: column;
    gap: 28px;
    align-items: center;
  }
  .gel-card {
    position: relative;
    user-select: none;
  }
  svg.gel-svg {
    display: block;
    max-width: 100%;
    max-height: 800px;
    width: auto;
    height: auto;
  }
  .crosshair-line {
    stroke: var(--accent-color);
    stroke-width: 1.5;
    stroke-dasharray: 4, 3;
    pointer-events: none;
    opacity: 0.95;
  }
  .crosshair-glow {
    stroke: var(--accent-color);
    stroke-width: 5;
    opacity: 0.3;
    pointer-events: none;
  }
  .pinned-line {
    stroke: #ffda6b;
    stroke-width: 1.5;
    stroke-dasharray: 2, 2;
    cursor: pointer;
  }
  .pinned-glow {
    stroke: #ffda6b;
    stroke-width: 4;
    opacity: 0.3;
  }
  .lane-rect {
    transition: opacity 0.15s ease;
  }
  .legend-entry {
    transition: opacity 0.15s ease;
  }
  .gel-tooltip {
    position: fixed;
    pointer-events: none;
    background: rgba(14, 18, 16, 0.95);
    backdrop-filter: blur(8px);
    -webkit-backdrop-filter: blur(8px);
    border: 1px solid var(--border-color);
    border-radius: 8px;
    padding: 8px 12px;
    font-size: 0.7rem;
    color: var(--text-color);
    box-shadow: 0 8px 24px rgba(0, 0, 0, 0.6), 0 0 10px rgba(94, 230, 115, 0.15);
    z-index: 1000;
    opacity: 0;
    transition: opacity 0.12s ease;
    white-space: nowrap;
  }
  .gel-tooltip.above {
    transform: translate(-50%, calc(-100% - 15px));
  }
  .gel-tooltip.below {
    transform: translate(-50%, 18px);
  }
  .gel-tooltip::after {
    content: "";
    position: absolute;
    left: 50%;
    transform: translateX(-50%);
    border-width: 6px;
    border-style: solid;
  }
  .gel-tooltip.above::after {
    top: 100%;
    border-color: rgba(14, 18, 16, 0.95) transparent transparent transparent;
  }
  .gel-tooltip.below::after {
    bottom: 100%;
    border-color: transparent transparent rgba(14, 18, 16, 0.95) transparent;
  }
  .gel-tooltip.visible {
    opacity: 1;
  }
  .tt-size {
    font-size: 0.9rem;
    font-weight: 700;
    color: var(--accent-color);
    margin-bottom: 2px;
  }
  .tt-detail {
    font-size: 0.66rem;
    color: rgba(255, 255, 255, 0.75);
    display: flex;
    gap: 8px;
    align-items: center;
  }
  .tt-lane {
    font-weight: 600;
    color: #ffffff;
  }
  .tt-closest {
    font-size: 0.63rem;
    color: rgba(255, 255, 255, 0.5);
    margin-top: 3px;
  }
  footer {
    margin-top: 30px;
    font-size: 0.64rem;
    color: rgba(255, 255, 255, 0.4);
    text-align: center;
  }
</style>
</head>
<body>

<header>
  <div class="title-group">
    <h1>__DOC_TITLE__</h1>
    <p>Interactive Agarose Gel Viewer &bull; Hover to measure fragment sizes across lanes &bull; Click to pin markers</p>
  </div>
  <div class="controls">
    <button class="btn" id="btn-clear-pins">&#10005; Clear Pins (<span id="pins-count">0</span>)</button>
    <div class="badge-chip">Theme: __THEME_NAME__</div>
  </div>
</header>

<div class="panels-container" id="panels-root"></div>

<div class="gel-tooltip" id="global-tooltip">
  <div class="tt-size" id="tt-size">--</div>
  <div class="tt-detail">
    <span class="tt-lane" id="tt-lane">Gel</span>
    <span id="tt-bp">--</span>
  </div>
  <div class="tt-closest" id="tt-closest">--</div>
</div>

<footer>
  Generated by <code>simulate-gel.py</code> &bull; Calibrated across __MIN_BP__&ndash;__MAX_BP__ bp scale
</footer>

<script>
const panels = __PANELS_JSON__;
const root = document.getElementById("panels-root");
const tooltip = document.getElementById("global-tooltip");
const ttSize = document.getElementById("tt-size");
const ttLane = document.getElementById("tt-lane");
const ttBp = document.getElementById("tt-bp");
const ttClosest = document.getElementById("tt-closest");
const pinsCountEl = document.getElementById("pins-count");
let totalPins = 0;

function yToBp(y, yTop, yBot, minBp, maxBp) {
  if (y <= yTop) return maxBp;
  if (y >= yBot) return minBp;
  const frac = (y - yTop) / (yBot - yTop);
  const logMin = Math.log10(minBp);
  const logMax = Math.log10(maxBp);
  const logBp = logMax - frac * (logMax - logMin);
  return Math.round(Math.pow(10, logBp));
}

function formatBp(bp) {
  if (bp >= 1000) {
    const kb = (bp / 1000);
    const digits = kb >= 10 ? 1 : 2;
    return `${kb.toFixed(digits)} kb`;
  }
  return `${bp} bp`;
}

function findClosestLadder(bp, ladderMarkers) {
  if (!ladderMarkers || ladderMarkers.length === 0) return null;
  let closest = ladderMarkers[0];
  let minDiff = Math.abs(bp - closest.bp);
  for (let i = 1; i < ladderMarkers.length; i++) {
    const diff = Math.abs(bp - ladderMarkers[i].bp);
    if (diff < minDiff) {
      minDiff = diff;
      closest = ladderMarkers[i];
    }
  }
  return closest;
}

panels.forEach((p, pIdx) => {
  const card = document.createElement("div");
  card.className = "gel-card";
  card.id = `panel-${pIdx}`;
  root.appendChild(card);

  const totalH = p.gel_h + 100;
  const svg = d3.select(card)
    .append("svg")
    .attr("class", "gel-svg")
    .attr("viewBox", `0 0 ${p.total_w} ${totalH}`)
    .attr("width", p.total_w)
    .attr("height", totalH);

  // Background
  svg.append("rect")
    .attr("width", p.total_w)
    .attr("height", totalH)
    .attr("fill", p.theme.bg_color);

  // Title if present
  if (p.panel_title) {
    svg.append("text")
      .attr("x", p.gel_w / 2.0)
      .attr("y", 28)
      .attr("text-anchor", "middle")
      .attr("fill", p.theme.text_color)
      .attr("font-size", 17)
      .attr("font-weight", "bold")
      .attr("font-family", "sans-serif")
      .text(p.panel_title);
  }

  // Gel Group offset for top title / margin
  const gY = 40;
  const g = svg.append("g").attr("transform", `translate(0, ${gY})`);

  // Embedded Photorealistic Gel Slab Image
  g.append("image")
    .attr("x", 0)
    .attr("y", 0)
    .attr("width", p.gel_w)
    .attr("height", p.gel_h)
    .attr("href", "data:image/png;base64," + p.gel_image_b64)
    .attr("preserveAspectRatio", "none");

  // Lane Hover Rectangles (subtle highlight behind bands)
  const laneHighlights = g.append("g").attr("class", "lane-highlights");
  p.lanes.forEach(lane => {
    laneHighlights.append("rect")
      .attr("class", `lane-rect lane-rect-${lane.lane_num}`)
      .attr("x", lane.x_start)
      .attr("y", 0)
      .attr("width", lane.x_end - lane.x_start)
      .attr("height", p.gel_h)
      .attr("fill", "white")
      .attr("opacity", 0)
      .attr("pointer-events", "none");
  });

  // Ladder ticks and labels
  const ladderG = g.append("g").attr("class", "ladder-markers");
  const ladderXc = p.lanes[0].xc;
  p.ladder_markers.forEach(lm => {
    ladderG.append("line")
      .attr("x1", ladderXc - p.lane_width / 2.0 - 14)
      .attr("x2", ladderXc - p.lane_width / 2.0 - 4)
      .attr("y1", lm.y)
      .attr("y2", lm.y)
      .attr("stroke", p.theme.text_color)
      .attr("stroke-width", 0.9)
      .attr("opacity", 0.8);

    ladderG.append("text")
      .attr("x", ladderXc - p.lane_width / 2.0 - 18)
      .attr("y", lm.y)
      .attr("text-anchor", "end")
      .attr("dominant-baseline", "central")
      .attr("fill", p.theme.text_color)
      .attr("font-size", 13)
      .attr("font-weight", "600")
      .attr("font-family", "sans-serif")
      .text(lm.label);
  });

  // Lane Numbers above wells and below gel
  const numsG = g.append("g").attr("class", "lane-numbers");
  p.lanes.forEach(lane => {
    numsG.append("text")
      .attr("x", lane.xc)
      .attr("y", p.y_well - 15)
      .attr("text-anchor", "middle")
      .attr("fill", p.theme.text_color)
      .attr("font-size", 15)
      .attr("font-weight", "bold")
      .attr("font-family", "sans-serif")
      .attr("opacity", 0.95)
      .text(lane.lane_num);

    numsG.append("text")
      .attr("x", lane.xc)
      .attr("y", p.gel_h + 24)
      .attr("text-anchor", "middle")
      .attr("dominant-baseline", "hanging")
      .attr("fill", p.theme.text_color)
      .attr("font-size", 15)
      .attr("font-weight", "bold")
      .attr("font-family", "sans-serif")
      .attr("opacity", 0.9)
      .text(lane.lane_num);
  });

  // Right-side Legend Card
  const legG = g.append("g").attr("class", "legend-card");
  const cardX = p.gel_w + 25;
  const cardW = p.legend_w - 45;
  const cardY = p.y_well - 5;
  const lineSpacing = 52;
  const cardH = 62 + p.legend_entries.length * lineSpacing + 18;

  legG.append("rect")
    .attr("x", cardX)
    .attr("y", cardY)
    .attr("width", cardW)
    .attr("height", cardH)
    .attr("rx", 8)
    .attr("ry", 8)
    .attr("fill", p.theme.well_color)
    .attr("stroke", p.theme.border_color)
    .attr("stroke-width", 1.2)
    .attr("opacity", 0.92);

  legG.append("text")
    .attr("x", cardX + 20)
    .attr("y", cardY + 28)
    .attr("dominant-baseline", "central")
    .attr("fill", p.theme.text_color)
    .attr("font-size", 12)
    .attr("font-weight", "bold")
    .attr("font-family", "sans-serif")
    .attr("opacity", 0.7)
    .text("LANES");

  legG.append("line")
    .attr("x1", cardX + 20)
    .attr("x2", cardX + cardW - 8)
    .attr("y1", cardY + 46)
    .attr("y2", cardY + 46)
    .attr("stroke", p.theme.border_color)
    .attr("stroke-width", 1.1)
    .attr("opacity", 0.85);

  p.legend_entries.forEach((entry, idx) => {
    const yPos = cardY + 72 + idx * lineSpacing;
    const entryG = legG.append("g")
      .attr("class", `legend-entry entry-${entry.lane_num}`)
      .style("cursor", "pointer")
      .on("mouseenter", () => highlightLane(entry.lane_num, true))
      .on("mouseleave", () => highlightLane(entry.lane_num, false));

    entryG.append("text")
      .attr("x", cardX + 36)
      .attr("y", yPos)
      .attr("text-anchor", "end")
      .attr("dominant-baseline", "central")
      .attr("fill", p.theme.text_color)
      .attr("font-size", 14)
      .attr("font-weight", "bold")
      .attr("font-family", "sans-serif")
      .text(entry.lane_num);

    entryG.append("text")
      .attr("x", cardX + 48)
      .attr("y", yPos)
      .attr("text-anchor", "middle")
      .attr("dominant-baseline", "central")
      .attr("fill", p.theme.text_color)
      .attr("font-size", 14)
      .attr("opacity", 0.75)
      .text("-");

    entryG.append("text")
      .attr("x", cardX + 62)
      .attr("y", yPos)
      .attr("text-anchor", "start")
      .attr("dominant-baseline", "central")
      .attr("fill", p.theme.text_color)
      .attr("font-size", 14)
      .attr("font-family", "sans-serif")
      .attr("opacity", 0.98)
      .text(entry.name);
  });

  // Pinned Lines Container
  const pinnedG = g.append("g").attr("class", "pinned-lines-group");

  // Crosshair Elements (horizontal line across the whole gel)
  const crosshairG = g.append("g").attr("class", "crosshair-group").style("display", "none");

  const crossGlow = crosshairG.append("line")
    .attr("class", "crosshair-glow")
    .attr("x1", 0)
    .attr("x2", p.gel_w);

  const crossLine = crosshairG.append("line")
    .attr("class", "crosshair-line")
    .attr("x1", 0)
    .attr("x2", p.gel_w);

  // Badge on the left of crosshair showing size (positioned above the line,
  // with a small gap so the box doesn't touch the crosshair itself)
  const badgeG = crosshairG.append("g").attr("class", "crosshair-badge");
  badgeG.append("rect")
    .attr("x", 4)
    .attr("y", -38)
    .attr("width", 90)
    .attr("height", 28)
    .attr("rx", 5)
    .attr("ry", 5)
    .attr("fill", p.theme.well_color)
    .attr("stroke", p.theme.accent_color)
    .attr("stroke-width", 1.2);

  const badgeText = badgeG.append("text")
    .attr("x", 49)
    .attr("y", -24)
    .attr("text-anchor", "middle")
    .attr("dominant-baseline", "central")
    .attr("fill", p.theme.accent_color)
    .attr("font-size", 13)
    .attr("font-weight", "bold")
    .attr("font-family", "sans-serif");

  // Interactive Overlay Rect covering the gel slab
  const overlay = g.append("rect")
    .attr("class", "interactive-overlay")
    .attr("x", 0)
    .attr("y", 0)
    .attr("width", p.gel_w)
    .attr("height", p.gel_h)
    .attr("fill", "transparent")
    .style("cursor", "crosshair");

  function highlightLane(laneNum, active) {
    d3.select(card).selectAll(`.lane-rect-${laneNum}`)
      .transition().duration(120)
      .attr("opacity", active ? 0.08 : 0);

    d3.select(card).selectAll(`.entry-${laneNum} text`)
      .transition().duration(120)
      .attr("fill", active ? p.theme.accent_color : p.theme.text_color);
  }

  overlay
    .on("mouseenter", () => {
      crosshairG.style("display", null);
      tooltip.classList.add("visible", "above");
    })
    .on("mouseleave", () => {
      crosshairG.style("display", "none");
      tooltip.classList.remove("visible");
      p.lanes.forEach(l => highlightLane(l.lane_num, false));
    })
    .on("mousemove", function(event) {
      const [mx, my] = d3.pointer(event, this);
      const bp = yToBp(my, p.y_top, p.y_bot, p.min_bp, p.max_bp);
      const formatted = formatBp(bp);

      crossGlow.attr("y1", my).attr("y2", my);
      crossLine.attr("y1", my).attr("y2", my);
      badgeG.attr("transform", `translate(0, ${my})`);
      badgeText.text(formatted);

      // Detect active lane
      let activeLane = null;
      for (const lane of p.lanes) {
        if (mx >= lane.x_start && mx <= lane.x_end) {
          activeLane = lane;
          break;
        }
      }

      p.lanes.forEach(l => highlightLane(l.lane_num, activeLane && activeLane.lane_num === l.lane_num));

      // Closest reference ladder band
      const closest = findClosestLadder(bp, p.ladder_markers);
      const diff = closest ? bp - closest.bp : 0;
      const diffStr = diff === 0 ? "Exact ladder match" :
                      diff > 0 ? `+${diff.toLocaleString()} bp above ${closest.label}` :
                      `${Math.abs(diff).toLocaleString()} bp below ${closest.label}`;

      ttSize.textContent = formatted;
      ttBp.textContent = `(${bp.toLocaleString()} bp)`;
      ttLane.textContent = activeLane ? (activeLane.is_ladder ? "Ladder Lane" : `Lane ${activeLane.lane_num}: ${activeLane.name}`) : "Gel Slab";
      ttClosest.textContent = `Ladder reference: ${closest ? closest.label : "--"} (${diffStr})`;

      // Flip the tooltip below the cursor when there isn't room above it,
      // so the arrow direction (and its border-color, which only .above/.below
      // define) always matches where the tooltip actually renders.
      const showBelow = event.clientY < 140;
      tooltip.classList.toggle("above", !showBelow);
      tooltip.classList.toggle("below", showBelow);

      tooltip.style.left = `${event.clientX}px`;
      tooltip.style.top = `${event.clientY}px`;
    })
    .on("click", function(event) {
      const [mx, my] = d3.pointer(event, this);
      const bp = yToBp(my, p.y_top, p.y_bot, p.min_bp, p.max_bp);
      const formatted = formatBp(bp);

      totalPins++;
      pinsCountEl.textContent = totalPins;

      const pinGroup = pinnedG.append("g")
        .attr("class", "pinned-marker")
        .attr("transform", `translate(0, ${my})`)
        .style("cursor", "pointer")
        .attr("title", "Click to remove pin");

      pinGroup.append("line")
        .attr("class", "pinned-glow")
        .attr("x1", 0)
        .attr("x2", p.gel_w)
        .attr("y1", 0)
        .attr("y2", 0);

      pinGroup.append("line")
        .attr("class", "pinned-line")
        .attr("x1", 0)
        .attr("x2", p.gel_w)
        .attr("y1", 0)
        .attr("y2", 0);

      pinGroup.append("rect")
        .attr("x", p.gel_w - 110)
        .attr("y", -14)
        .attr("width", 110)
        .attr("height", 28)
        .attr("rx", 4)
        .attr("ry", 4)
        .attr("fill", "#2a2205")
        .attr("stroke", "#ffda6b")
        .attr("stroke-width", 1);

      pinGroup.append("text")
        .attr("x", p.gel_w - 55)
        .attr("y", 0)
        .attr("text-anchor", "middle")
        .attr("dominant-baseline", "central")
        .attr("fill", "#ffda6b")
        .attr("font-size", 13)
        .attr("font-weight", "bold")
        .attr("font-family", "sans-serif")
        .text(`📌 ${formatted}`);

      pinGroup.on("click", function(ev) {
        ev.stopPropagation();
        d3.select(this).remove();
        totalPins = Math.max(0, totalPins - 1);
        pinsCountEl.textContent = totalPins;
      });
    });
});

document.getElementById("btn-clear-pins").addEventListener("click", () => {
  d3.selectAll(".pinned-marker").remove();
  totalPins = 0;
  pinsCountEl.textContent = "0";
});
</script>
</body>
</html>"""


def synthesize_gel(ladder_name, ladder_bands, samples_chunk, samples_per_gel,
                   min_bp, max_bp, gel_h, band_sigma, smear_sigma, ladder_sigma, smile_amp, lane_gap=22):
    """Compute physical gel image array and layout geometry."""
    lane_width = 76
    margin_left = 90
    margin_right = 35

    n_lanes = samples_per_gel + 1
    gel_w = margin_left + n_lanes * lane_width + (n_lanes - 1) * lane_gap + margin_right

    # Scale vertical proportions relative to standard 880px height
    y_well = int(round(gel_h * (65.0 / 880.0)))
    well_h = int(round(gel_h * (16.0 / 880.0)))
    well_w = int(lane_width * 0.76)
    well_gap = max(1, int(round(gel_h * (8.0 / 880.0))))
    y_top = y_well + well_h + well_gap
    y_bot = gel_h - int(round(gel_h * (45.0 / 880.0)))

    # Background illumination vignette + agarose matrix speckle
    yy, xx = np.mgrid[0:gel_h, 0:gel_w]
    cx, cy = gel_w / 2.0, gel_h * 0.45
    rad = ((xx - cx) / (gel_w / 2.0)) ** 2 + ((yy - cy) / (gel_h / 2.0)) ** 2
    bg = 0.045 * (1.0 - 0.15 * rad)
    rng = np.random.default_rng(42)
    noise = rng.normal(0, 0.003, (gel_h, gel_w))
    gel_img = np.clip(bg + noise, 0.01, 1.0)

    # Gel slab borders & chamfer edge highlight
    slab_pad = 12
    mask_outside = (xx < slab_pad) | (xx > gel_w - slab_pad) | (yy < slab_pad) | (yy > gel_h - slab_pad)
    gel_img[mask_outside] *= 0.35

    edge_x = (np.abs(xx - slab_pad) <= 1) | (np.abs(xx - (gel_w - slab_pad)) <= 1)
    edge_y = (np.abs(yy - slab_pad) <= 1) | (np.abs(yy - (gel_h - slab_pad)) <= 1)
    gel_img[(edge_x | edge_y) & ~mask_outside] += 0.03

    # Assemble lanes
    ladder_prof = ladder_profile(ladder_bands, min_bp, max_bp, y_top, y_bot, gel_h, ladder_sigma)
    all_lanes = [(f"{ladder_name} ladder", ladder_prof, False, False, 0.0)]

    for path, lengths in samples_chunk:
        prof, well_frac = profile_from_lengths(lengths, min_bp, max_bp, y_top, y_bot, gel_h, band_sigma, smear_sigma)
        all_lanes.append((clean_sample_name(path), prof, True, False, well_frac))

    # Pad to fixed comb capacity if needed
    while len(all_lanes) - 1 < samples_per_gel:
        all_lanes.append(("", np.zeros(gel_h), False, True, 0.0))

    lane_centers = []
    y_indices = np.arange(gel_h, dtype=float)

    for i, (name, prof, is_sample, is_empty, well_frac) in enumerate(all_lanes):
        x_start = margin_left + i * (lane_width + lane_gap)
        x_center = x_start + lane_width / 2.0
        lane_centers.append(x_center)

        # Well pocket rendering
        w_x0 = int(round(x_center - well_w / 2.0))
        w_x1 = w_x0 + well_w
        w_y0 = y_well
        w_y1 = y_well + well_h

        gel_img[w_y0:w_y1, w_x0:w_x1] = 0.015 + rng.normal(0, 0.0015, (well_h, well_w))
        gel_img[w_y0:w_y0 + 1, w_x0:w_x1] = 0.05
        gel_img[w_y0:w_y1, w_x0:w_x0 + 1] = 0.04
        gel_img[w_y0:w_y1, w_x1 - 1:w_x1] = 0.04

        # Well loading line (residual dye for loaded samples)
        if is_sample and np.max(prof) > 0.1:
            gel_img[w_y1 - 1:w_y1 + 1, w_x0 + 2:w_x1 - 2] += 0.06

        # Fragments bigger than max_bp never entered the gel — they stay
        # trapped in the well pocket as a retained band. Brightness scales
        # with how much of the sample's mass was too large to migrate; shape
        # follows the same horizontal meniscus and grain as an in-gel band,
        # concentrated toward the bottom of the well (closest to the gel
        # interface it couldn't cross) rather than filling it as a flat block.
        if is_sample and well_frac > 0.0:
            peak_glow = 0.92 * (np.clip(well_frac, 0.0, 1.0) ** 0.55)

            x_local = np.arange(w_x0, w_x1) - x_center
            h_prof = np.exp(-(np.abs(x_local) / (well_w / 2.0 * 0.88)) ** 6)

            y_local = np.arange(w_y0, w_y1) - w_y0
            v_prof = (0.35 + 0.65 * y_local / max(1, well_h - 1)) ** 1.4

            band = peak_glow * np.outer(v_prof, h_prof)
            band += rng.normal(0, 0.012 * peak_glow, band.shape)
            gel_img[w_y0:w_y1, w_x0:w_x1] += band

            # Faint diffuse leading edge: the smallest of the oversized
            # fragments start creeping toward the gel before stalling,
            # softening the boundary between well and slab.
            bleed_h = max(1, min(y_top - w_y1 - 1, well_h // 3))
            if bleed_h > 0:
                fade = np.linspace(1.0, 0.0, bleed_h) ** 1.5
                bleed = (peak_glow * 0.3 * h_prof)[None, :] * fade[:, None]
                gel_img[w_y1:w_y1 + bleed_h, w_x0:w_x1] += bleed

        if is_empty or np.max(prof) == 0:
            continue

        # 2D lane rendering with flat core and meniscus smile curvature
        for x in range(x_start, x_start + lane_width):
            u = (x - x_center) / (lane_width / 2.0)
            if abs(u) > 1.0:
                continue
            h_weight = np.exp(-(abs(u) / 0.82) ** 8)
            y_shifted = y_indices + smile_amp * (u ** 2)
            lane_col = np.interp(y_shifted, y_indices, prof)
            gel_img[:, x] += lane_col * h_weight * 0.92

    # Optical bloom & sensor saturation
    bloom_source = np.maximum(0, gel_img - 0.55) ** 1.6
    for r in range(gel_h):
        bloom_source[r, :] = smooth1d(bloom_source[r, :], 2.5)
    for c in range(gel_w):
        bloom_source[:, c] = smooth1d(bloom_source[:, c], 2.5)
    gel_img += bloom_source * 0.35
    gel_img = np.clip(gel_img, 0.0, 1.0)

    return (gel_img, gel_w, lane_centers, all_lanes,
            y_well, well_h, well_w, y_top, y_bot,
            margin_left, lane_width)


def render_panel(ax, ladder_name, ladder_bands, samples_chunk, samples_per_gel,
                 min_bp, max_bp, gel_h, band_sigma, smear_sigma, ladder_sigma,
                 smile_amp, theme, legend_w, panel_title=None, lane_gap=22):
    """Render a single agarose gel slab into the specified Matplotlib Axes."""
    (gel_img, gel_w, lane_centers, all_lanes,
     y_well, well_h, well_w, y_top, y_bot,
     margin_left, lane_width) = synthesize_gel(
        ladder_name, ladder_bands, samples_chunk, samples_per_gel,
        min_bp, max_bp, gel_h, band_sigma, smear_sigma, ladder_sigma, smile_amp, lane_gap
    )

    # Plot onto axes
    ax.set_facecolor(theme["bg_color"])
    ax.imshow(gel_img, cmap=theme["cmap"], aspect="auto", vmin=0, vmax=1,
              extent=[0, gel_w, gel_h, 0])

    # Molecular weight ladder markers (left side of ladder)
    ladder_xc = lane_centers[0]
    for bp, _ in ladder_bands:
        if bp < min_bp or bp > max_bp:
            continue
        y = bp_to_y(bp, min_bp, max_bp, y_top, y_bot)
        label = f"{bp/1000:g} kb" if bp >= 1000 else f"{bp} bp"
        ax.plot([ladder_xc - lane_width / 2.0 - 14, ladder_xc - lane_width / 2.0 - 4], [y, y],
                color=theme["text_color"], lw=0.9, alpha=0.8)
        ax.text(ladder_xc - lane_width / 2.0 - 18, y, label,
                color=theme["text_color"], fontsize=13.0, ha="right", va="center",
                fontfamily="sans-serif", fontweight="semibold")

    # Lane numbers above wells and below gel slab
    for i, xc in enumerate(lane_centers):
        lane_str = str(i + 1)
        ax.text(xc, y_well - 15, lane_str, color=theme["text_color"],
                fontsize=15.0, fontweight="bold", ha="center", va="bottom",
                fontfamily="sans-serif", alpha=0.95)
        ax.text(xc, gel_h + 24, lane_str, color=theme["text_color"],
                fontsize=15.0, fontweight="bold", ha="center", va="top",
                fontfamily="sans-serif", alpha=0.9)

    # Legend on the right: "number - sample name"
    legend_entries = []
    for i, (name, _, _, is_empty, _) in enumerate(all_lanes):
        if not is_empty and name:
            legend_entries.append((i + 1, name))

    if legend_entries:
        card_x = gel_w + 25
        card_w = legend_w - 45
        y_start = y_well - 5
        line_spacing = 52
        card_h = 62 + len(legend_entries) * line_spacing + 18

        # Card container with subtle border matching theme
        rect = mpatches.FancyBboxPatch(
            (card_x, y_start), card_w, card_h,
            boxstyle="round,pad=8,rounding_size=8",
            linewidth=1.2,
            edgecolor=theme["border_color"],
            facecolor=theme["well_color"],
            alpha=0.92,
            zorder=3
        )
        ax.add_patch(rect)

        # Header
        ax.text(card_x + 20, y_start + 28, "LANES", color=theme["text_color"],
                fontsize=12.0, fontweight="bold", fontfamily="sans-serif",
                alpha=0.7, va="center", zorder=4)

        # Divider line
        ax.plot([card_x + 20, card_x + card_w - 8], [y_start + 46, y_start + 46],
                color=theme["border_color"], lw=1.1, alpha=0.85, zorder=4)

        # Entries: number - sample name
        for idx, (lane_num, sname) in enumerate(legend_entries):
            y_pos = y_start + 72 + idx * line_spacing
            ax.text(card_x + 36, y_pos, str(lane_num), color=theme["text_color"],
                    fontsize=14.0, fontweight="bold", fontfamily="sans-serif",
                    ha="right", va="center", zorder=4)
            ax.text(card_x + 48, y_pos, "-", color=theme["text_color"],
                    fontsize=14.0, fontweight="normal", fontfamily="sans-serif",
                    ha="center", va="center", alpha=0.75, zorder=4)
            ax.text(card_x + 62, y_pos, sname, color=theme["text_color"],
                    fontsize=14.0, fontweight="normal", fontfamily="sans-serif",
                    ha="left", va="center", zorder=4, alpha=0.98)

    total_w = gel_w + legend_w
    ax.set_xlim(0, total_w)
    ax.set_ylim(gel_h + 58, -36)
    ax.axis("off")

    if panel_title:
        ax.text(gel_w / 2.0, -18, panel_title, color=theme["text_color"],
                fontsize=17.0, ha="center", va="bottom",
                fontfamily="sans-serif", fontweight="bold")

    return total_w


def build_panel_d3_data(ladder_name, ladder_bands, samples_chunk, samples_per_gel,
                        min_bp, max_bp, gel_h, band_sigma, smear_sigma, ladder_sigma,
                        smile_amp, theme, theme_name, legend_w, panel_title=None, lane_gap=22):
    """Build data packet and base64 gel image for D3 HTML visualization."""
    (gel_img, gel_w, lane_centers, all_lanes,
     y_well, well_h, well_w, y_top, y_bot,
     margin_left, lane_width) = synthesize_gel(
        ladder_name, ladder_bands, samples_chunk, samples_per_gel,
        min_bp, max_bp, gel_h, band_sigma, smear_sigma, ladder_sigma, smile_amp, lane_gap
    )

    # Encode gel image to PNG base64
    buf = io.BytesIO()
    plt.imsave(buf, gel_img, cmap=theme["cmap"], format="png", vmin=0, vmax=1)
    buf.seek(0)
    img_b64 = base64.b64encode(buf.read()).decode("ascii")

    # Ladder markers
    ladder_markers = []
    for bp, intensity in ladder_bands:
        if min_bp <= bp <= max_bp:
            y = bp_to_y(bp, min_bp, max_bp, y_top, y_bot)
            label = f"{bp/1000:g} kb" if bp >= 1000 else f"{bp} bp"
            ladder_markers.append({
                "bp": int(bp),
                "y": float(y),
                "label": label,
                "intensity": float(intensity)
            })

    # Lanes metadata
    lanes_data = []
    legend_entries = []
    for i, (name, _, is_sample, is_empty, _) in enumerate(all_lanes):
        xc = lane_centers[i]
        lanes_data.append({
            "lane_num": i + 1,
            "name": name,
            "xc": float(xc),
            "x_start": float(xc - lane_width / 2.0),
            "x_end": float(xc + lane_width / 2.0),
            "is_ladder": (i == 0),
            "is_empty": bool(is_empty)
        })
        if not is_empty and name:
            legend_entries.append({"lane_num": i + 1, "name": name})

    return {
        "panel_title": panel_title,
        "gel_image_b64": img_b64,
        "gel_w": int(gel_w),
        "gel_h": int(gel_h),
        "legend_w": int(legend_w),
        "total_w": int(gel_w + legend_w),
        "lane_width": int(lane_width),
        "y_well": int(y_well),
        "well_h": int(well_h),
        "y_top": float(y_top),
        "y_bot": float(y_bot),
        "min_bp": float(min_bp),
        "max_bp": float(max_bp),
        "ladder_markers": ladder_markers,
        "lanes": lanes_data,
        "legend_entries": legend_entries,
        "theme": {
            "name": theme_name,
            "bg_color": theme["bg_color"],
            "text_color": theme["text_color"],
            "well_color": theme["well_color"],
            "border_color": theme["border_color"],
            "accent_color": theme["accent_color"],
        }
    }


def render_html(panels_data, out_path, title=None):
    """Render interactive D3 HTML document."""
    theme = panels_data[0]["theme"]
    doc_title = title or ""
    panels_json = json.dumps(panels_data)

    content = HTML_TEMPLATE
    replacements = {
        "__DOC_TITLE__": doc_title,
        "__PANELS_JSON__": panels_json,
        "__BG_COLOR__": theme["bg_color"],
        "__TEXT_COLOR__": theme["text_color"],
        "__ACCENT_COLOR__": theme["accent_color"],
        "__BORDER_COLOR__": theme["border_color"],
        "__WELL_COLOR__": theme["well_color"],
        "__THEME_NAME__": theme["name"],
        "__MIN_BP__": f"{panels_data[0]['min_bp']:g}",
        "__MAX_BP__": f"{panels_data[0]['max_bp']:g}",
    }
    for key, val in replacements.items():
        content = content.replace(key, str(val))

    with open(out_path, "w", encoding="utf-8") as f:
        f.write(content)


def render(lengths_by_file, ladder_name, out_path, min_bp, max_bp,
           n_pixels, band_sigma, smear_sigma, ladder_sigma, smile_amp,
           theme_name, title, lanes, pad_single=False, out_format="png", invert=False, lane_gap=22):
    """Render one or multiple vertically stacked gels as static PNG or interactive HTML."""
    theme = THEMES.get(theme_name, THEMES["sybr"]).copy()
    if invert:
        theme["cmap"] = theme["cmap"].reversed()
        theme["bg_color"] = "#ffffff"
        theme["text_color"] = "#111111"
        theme["well_color"] = "#f5f5f5"
        theme["border_color"] = "#cccccc"
        theme["accent_color"] = "#d9534f" if theme_name == "chemidoc" else "#2b8a3e"
    ladder_bands = LADDERS[ladder_name]
    samples_per_gel = max(1, lanes - 1)

    chunks = [lengths_by_file[i:i + samples_per_gel]
              for i in range(0, len(lengths_by_file), samples_per_gel)] or [[]]
    n_panels = len(chunks)
    effective_lanes = samples_per_gel if (n_panels > 1 or pad_single) else len(chunks[0])

    lane_width = 76
    margin_left = 90
    margin_right = 35

    # Compute legend width dynamically based on longest sample name
    all_names = [f"{ladder_name} ladder"]
    for path, _ in lengths_by_file:
        all_names.append(clean_sample_name(path))
    max_len = max([len(n) for n in all_names], default=10)
    legend_w = max(330, int(max_len * 12.5) + 130)

    if out_format == "html":
        panels_data = []
        for idx, chunk in enumerate(chunks):
            if n_panels > 1:
                start = idx * samples_per_gel + 1
                end = idx * samples_per_gel + len(chunk)
                panel_title = f"Gel {idx + 1}/{n_panels} — samples {start}-{end}"
                if title:
                    panel_title = f"{title}\n{panel_title}"
            else:
                panel_title = title

            pdata = build_panel_d3_data(
                ladder_name=ladder_name,
                ladder_bands=ladder_bands,
                samples_chunk=chunk,
                samples_per_gel=effective_lanes,
                min_bp=min_bp,
                max_bp=max_bp,
                gel_h=n_pixels,
                band_sigma=band_sigma,
                smear_sigma=smear_sigma,
                ladder_sigma=ladder_sigma,
                smile_amp=smile_amp,
                theme=theme,
                theme_name=theme_name,
                legend_w=legend_w,
                panel_title=panel_title,
                lane_gap=lane_gap
            )
            panels_data.append(pdata)

        render_html(panels_data, out_path, title=title)
        return

    # PNG rendering
    n_total_lanes = effective_lanes + 1
    gel_w = margin_left + n_total_lanes * lane_width + (n_total_lanes - 1) * lane_gap + margin_right
    total_w = gel_w + legend_w

    fig_w = total_w / 90.0
    panel_h = (n_pixels + 120) / 90.0
    fig, axes = plt.subplots(nrows=n_panels, ncols=1, figsize=(fig_w, panel_h * n_panels),
                             dpi=120, squeeze=False)
    fig.patch.set_facecolor(theme["bg_color"])

    for idx, chunk in enumerate(chunks):
        if n_panels > 1:
            start = idx * samples_per_gel + 1
            end = idx * samples_per_gel + len(chunk)
            panel_title = f"Gel {idx + 1}/{n_panels} — samples {start}-{end}"
            if title:
                panel_title = f"{title}\n{panel_title}"
        else:
            panel_title = title

        render_panel(
            ax=axes[idx, 0],
            ladder_name=ladder_name,
            ladder_bands=ladder_bands,
            samples_chunk=chunk,
            samples_per_gel=effective_lanes,
            min_bp=min_bp,
            max_bp=max_bp,
            gel_h=n_pixels,
            band_sigma=band_sigma,
            smear_sigma=smear_sigma,
            ladder_sigma=ladder_sigma,
            smile_amp=smile_amp,
            theme=theme,
            legend_w=legend_w,
            panel_title=panel_title,
            lane_gap=lane_gap
        )

    plt.tight_layout()
    plt.savefig(out_path, facecolor=theme["bg_color"], edgecolor="none", bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Simulate a photorealistic agarose gel image from FASTQ read-length distributions."
    )
    parser.add_argument("fastq", nargs="+", help="Input FASTQ file(s) (.fastq or .fastq.gz), one lane each")
    parser.add_argument("-o", "--output", default="simulated_gel.png", help="Output path (default: simulated_gel.png)")
    parser.add_argument("--format", choices=["png", "html"], default=None,
                        help="Output format: 'png' for static image, 'html' for interactive D3 plot. "
                             "If not specified, inferred from the extension of --output (default: png).")
    parser.add_argument("--theme", choices=sorted(THEMES), default="sybr",
                        help="Gel visualization theme (sybr, ethidium, chemidoc; default: sybr)")
    parser.add_argument("--ladder", choices=sorted(LADDERS), default="1kb",
                        help="Ladder to use (default: 1kb)")
    parser.add_argument("--lanes", type=int, default=12,
                        help="Max lanes per comb, including the ladder (default: 12). "
                             "If more samples than that are given, they're split across "
                             "multiple gels stacked vertically in one image, each with its own ladder.")
    parser.add_argument("--pad", action="store_true",
                        help="Pad single gel to full comb width with empty wells (default: compact width)")
    parser.add_argument("--min-bp", type=float, default=100, help="Smallest fragment size shown (default: 100)")
    parser.add_argument("--max-bp", type=float, default=12000, help="Largest fragment size shown (default: 12000)")
    parser.add_argument("--pixels", type=int, default=680, help="Vertical resolution of each gel (default: 680)")
    parser.add_argument("--band-sigma", type=float, default=3.2,
                        help="Band sharpness sigma for concentrated amplicon peaks (default: 3.2)")
    parser.add_argument("--smear-sigma", type=float, default=20.0,
                        help="Smear smoothness sigma for dispersed libraries (default: 20.0)")
    parser.add_argument("--ladder-sigma", type=float, default=2.6,
                        help="Ladder band sharpness sigma (default: 2.6)")
    parser.add_argument("--smile", type=float, default=3.2,
                        help="Meniscus smile curvature amplitude in pixels (default: 3.2)")
    parser.add_argument("--lane-gap", type=int, default=15,
                        help="Spacing in pixels between sample lanes (default: 15)")
    parser.add_argument("--invert", action="store_true",
                        help="Invert colors (white background, dark DNA bands)")
    parser.add_argument("--title", default=None, help="Optional title printed above the gel(s)")
    args = parser.parse_args()

    if args.lanes < 2:
        sys.exit("Error: --lanes must be at least 2 (1 ladder + 1 sample)")

    ext = os.path.splitext(args.output)[1].lower()
    if args.format is not None:
        out_format = args.format
    elif ext in (".html", ".htm"):
        out_format = "html"
    else:
        out_format = "png"

    if out_format == "html" and args.output == "simulated_gel.png":
        args.output = "simulated_gel.html"

    lengths_by_file = []
    for path in args.fastq:
        if not os.path.exists(path):
            sys.exit(f"Error: file not found: {path}")
        try:
            lengths = read_fastq_lengths(path)
        except Exception as e:
            sys.exit(f"Error reading {path}: {e}")
        lengths_by_file.append((path, lengths))
        print(f"{path}: {len(lengths)} reads, "
              f"length range {min(lengths)}-{max(lengths)} bp, "
              f"median {int(np.median(lengths))} bp")

    n_gels = -(-len(lengths_by_file) // max(1, args.lanes - 1))  # ceil div
    if n_gels > 1:
        print(f"{len(lengths_by_file)} samples > {args.lanes - 1} per comb "
              f"-> splitting into {n_gels} gels stacked vertically")

    render(
        lengths_by_file=lengths_by_file,
        ladder_name=args.ladder,
        out_path=args.output,
        min_bp=args.min_bp,
        max_bp=args.max_bp,
        n_pixels=args.pixels,
        band_sigma=args.band_sigma,
        smear_sigma=args.smear_sigma,
        ladder_sigma=args.ladder_sigma,
        smile_amp=args.smile,
        theme_name=args.theme,
        title=args.title,
        lanes=args.lanes,
        pad_single=args.pad,
        out_format=out_format,
        invert=args.invert,
        lane_gap=args.lane_gap
    )
    if out_format == "html":
        print(f"Saved interactive D3 gel visualization to {args.output}")
    else:
        print(f"Saved gel image to {args.output}")


if __name__ == "__main__":
    main()
