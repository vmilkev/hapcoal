#!/usr/bin/env python3
"""Plot hapcoal results: LD decay curves and genetic maps.

Usage:
    python3 scripts/plot_ld.py [OUT_DIR] [--chrom C]      (default OUT_DIR = out/chr22)

Writes into OUT_DIR:
    ld_decay.{png|svg}        observed_*.ld_decay.tsv and sim_*.ld_decay.tsv together
    <run>.maps.{png|svg}      for every <run>.map.true.tsv: true vs realised vs
                              LD-inferred genetic map (cumulative cM and rate)

Uses matplotlib when installed (PNG); otherwise writes dependency-free SVG.
"""
import bisect
import csv
import glob
import math
import os
import sys

PALETTE = ["#111111", "#4269d0", "#efb118", "#ff725c", "#6cc5b0", "#3ca951", "#ff8ab7", "#a463f2", "#97bbf5"]


# ------------------------------------------------------------------ data ----
def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def read_decay(path):
    xs, ys = [], []
    for row in read_tsv(path):
        if row["chrom"] == "all":
            continue
        xs.append(float(row["bin_mid_bp"]) / 1000.0)
        ys.append(float(row["mean_r2"]))
    return xs, ys


def pretty_label(name):
    name = name.replace("smcp", "SMC'").replace("smc", "SMC").replace("macs", "MaCS")
    return name.replace("_ld_fit", " / LD-fitted demography").replace("_tennessen2012", " / Tennessen 2012")


# --------------------------------------------------------------- figures ----
# A figure is a list of panels; a panel is a dict with title, xlabel, ylabel,
# logx, and series = [(x, y, label, colour, dashed, width)].
def render(panels, path_noext, title):
    try:
        return render_matplotlib(panels, path_noext, title)
    except ImportError:
        return render_svg(panels, path_noext, title)


def render_matplotlib(panels, path_noext, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(len(panels), 1, figsize=(9, 4.6 * len(panels)), dpi=130, squeeze=False)
    for ax, p in zip(axes[:, 0], panels):
        for x, y, label, col, dashed, width in p["series"]:
            ax.plot(x, y, color=col, lw=width, ls="--" if dashed else "-", label=label)
        if p.get("logx"):
            ax.set_xscale("log")
        ax.set_xlabel(p["xlabel"])
        ax.set_ylabel(p["ylabel"])
        ax.set_title(p["title"], fontsize=11)
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle(title, fontweight="bold")
    fig.tight_layout()
    path = path_noext + ".png"
    fig.savefig(path)
    plt.close(fig)
    return path


def nice_ticks(lo, hi, n=6):
    if hi <= lo:
        hi = lo + 1
    raw = (hi - lo) / n
    mag = 10 ** math.floor(math.log10(raw))
    step = min((m * mag for m in (1, 2, 2.5, 5, 10) if m * mag >= raw), default=10 * mag)
    t = math.ceil(lo / step) * step
    out = []
    while t <= hi + 1e-9 * step:
        out.append(round(t, 10))
        t += step
    return out


def render_svg(panels, path_noext, title):
    W, PH, L, R, T, B = 900, 430, 75, 30, 50, 55
    H = T + PH * len(panels) + 10
    s = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" font-family="sans-serif" font-size="12">',
         f'<rect width="{W}" height="{H}" fill="white"/>',
         f'<text x="{L}" y="26" font-size="15" font-weight="bold">{title}</text>']
    for k, p in enumerate(panels):
        top = T + k * PH
        y0, y1 = top + 28, top + PH - B          # plot area (y)
        xs = [v for x, *_ in p["series"] for v in x]
        ys = [v for _, y, *_ in p["series"] for v in y]
        logx = p.get("logx", False)
        if logx:
            xlo, xhi = math.log10(max(min(v for v in xs if v > 0), 1e-12)), math.log10(max(xs))
        else:
            xlo, xhi = min(xs), max(xs)
        ylo, yhi = min(0.0, min(ys)), max(ys) * 1.05
        fx = lambda v: L + ((math.log10(max(v, 1e-12)) if logx else v) - xlo) / ((xhi - xlo) or 1) * (W - L - R)
        fy = lambda v: y1 - (v - ylo) / ((yhi - ylo) or 1) * (y1 - y0)
        s.append(f'<text x="{L}" y="{top + 16}" font-size="13" font-weight="bold">{p["title"]}</text>')
        xt = [v for v in (0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000) if 10 ** xlo <= v <= 10 ** xhi] if logx \
            else nice_ticks(xlo, xhi)
        for v in xt:
            s.append(f'<line x1="{fx(v):.1f}" y1="{y0}" x2="{fx(v):.1f}" y2="{y1}" stroke="#e3e3e3"/>')
            s.append(f'<text x="{fx(v):.1f}" y="{y1 + 16}" text-anchor="middle">{v:g}</text>')
        for v in nice_ticks(ylo, yhi):
            s.append(f'<line x1="{L}" y1="{fy(v):.1f}" x2="{W - R}" y2="{fy(v):.1f}" stroke="#e3e3e3"/>')
            s.append(f'<text x="{L - 8}" y="{fy(v) + 4:.1f}" text-anchor="end">{v:g}</text>')
        s.append(f'<text x="{(L + W - R) / 2}" y="{y1 + 36}" text-anchor="middle">{p["xlabel"]}</text>')
        s.append(f'<text transform="translate(18,{(y0 + y1) / 2}) rotate(-90)" text-anchor="middle">{p["ylabel"]}</text>')
        for x, y, label, col, dashed, width in p["series"]:
            pts = " ".join(f"{fx(a):.1f},{fy(b):.1f}" for a, b in zip(x, y))
            dash = ' stroke-dasharray="6 4"' if dashed else ""
            s.append(f'<polyline fill="none" stroke="{col}" stroke-width="{width}"{dash} points="{pts}"/>')
        # legend on a white box so it stays readable over the curves
        lx0 = p.get("legend_x", 0.62) * W
        wmax = max(len(lab) for _, _, lab, *_ in p["series"]) * 6.6 + 44
        s.append(f'<rect x="{lx0 - 6}" y="{y0 + 2}" width="{wmax:.0f}" height="{18 * len(p["series"]) + 6}" '
                 f'fill="white" fill-opacity="0.9" stroke="#ddd"/>')
        for i, (x, y, label, col, dashed, width) in enumerate(p["series"]):
            dash = ' stroke-dasharray="6 4"' if dashed else ""
            lx, ly = lx0, y0 + 14 + 18 * i
            s.append(f'<line x1="{lx}" y1="{ly}" x2="{lx + 26}" y2="{ly}" stroke="{col}" stroke-width="2.5"{dash}/>')
            s.append(f'<text x="{lx + 32}" y="{ly + 4}">{label}</text>')
    s.append("</svg>")
    path = path_noext + ".svg"
    with open(path, "w") as f:
        f.write("\n".join(s))
    return path


# ------------------------------------------------------------- LD decay ----
def plot_decay(out):
    series = []
    for p in sorted(glob.glob(os.path.join(out, "observed_*.ld_decay.tsv"))):
        series.append(("observed " + os.path.basename(p).split(".")[0][len("observed_"):], p))
    for p in sorted(glob.glob(os.path.join(out, "sim_*.ld_decay.tsv"))):
        series.append((pretty_label(os.path.basename(p)[len("sim_"):-len(".ld_decay.tsv")]), p))
    if not series:
        return None
    lines = []
    for i, (label, p) in enumerate(series):
        x, y = read_decay(p)
        obs = label.startswith("observed")
        lines.append((x, y, label, PALETTE[i % len(PALETTE)], not obs and "Tennessen" in label, 2.8 if obs else 1.5))
    panel = dict(title="mean r² of SNP pairs (MAF ≥ 0.05) by distance", xlabel="distance between SNPs (kb, log scale)",
                 ylabel="mean r²", logx=True, series=lines, legend_x=0.60)
    return render([panel], os.path.join(out, "ld_decay"), "LD decay: observed vs simulated")


# ---------------------------------------------------------- genetic maps ----
def plot_maps(true_path, chrom=None):
    base = true_path[:-len(".map.true.tsv")]
    run = os.path.basename(base)
    true = read_tsv(true_path)
    chrom = chrom or true[0]["chrom"]
    true = [r for r in true if r["chrom"] == chrom]
    tpos = [float(r["pos_bp"]) for r in true]
    tcm = [float(r["pos_cM"]) for r in true]

    def true_cm(x):  # interpolate the true map
        i = bisect.bisect_right(tpos, x) - 1
        if i < 0:
            return tcm[0]
        if i >= len(tpos) - 1:
            return tcm[-1]
        return tcm[i] + (x - tpos[i]) / (tpos[i + 1] - tpos[i]) * (tcm[i + 1] - tcm[i])

    def windows(path):
        if not os.path.exists(path):
            return []
        return [(float(r["start_bp"]), float(r["end_bp"]), float(r["rate_cM_per_Mb"]))
                for r in read_tsv(path) if r["chrom"] == chrom]

    def cumulative(win, scale=1.0):
        x, y, c = [win[0][0] / 1e6], [0.0], 0.0
        for a, b, rate in win:
            c += rate * scale * (b - a) / 1e6
            x.append(b / 1e6)
            y.append(c)
        return x, y

    def steps(win, scale=1.0):
        x, y = [], []
        for a, b, rate in win:
            x += [a / 1e6, b / 1e6]
            y += [rate * scale, rate * scale]
        return x, y

    real = windows(base + ".map.realised.tsv")
    ldw = windows(base + ".map.ld.tsv")
    total = tcm[-1] - tcm[0]
    cum_s, rate_s = [], []

    cum_s.append(([p / 1e6 for p in tpos], [c - tcm[0] for c in tcm], "true map (simulation input)", PALETTE[0], False, 2.6))
    if real:
        x, y = cumulative(real)
        cum_s.append((x, y, "realised (simulated recombination events)", PALETTE[1], False, 1.6))
        # true rate averaged over the same windows, for a like-for-like comparison
        tw = [(a, b, (true_cm(b) - true_cm(a)) / ((b - a) / 1e6)) for a, b, _ in real]
        rate_s.append((*steps(tw), "true map (window average)", PALETTE[0], False, 2.2))
        rate_s.append((*steps(real), "realised", PALETTE[1], False, 1.4))
    if ldw:
        est = sum(r * (b - a) / 1e6 for a, b, r in ldw)
        k = total / est if est > 0 else 1.0
        x, y = cumulative(ldw, k)
        lab = f"LD-inferred from haplotypes (x{k:.2f}, scaled to map length)"
        cum_s.append((x, y, lab, PALETTE[3], True, 1.8))
        rate_s.append((*steps(ldw, k), lab, PALETTE[3], True, 1.8))

    panels = [dict(title=f"chr{chrom}: cumulative genetic position", xlabel="position (Mb)", ylabel="cM",
                   series=cum_s, legend_x=0.09)]
    if rate_s:
        panels.append(dict(title=f"chr{chrom}: recombination rate per window", xlabel="position (Mb)",
                           ylabel="cM/Mb", series=rate_s, legend_x=0.09))
    return render(panels, base + ".maps", f"Genetic maps: {pretty_label(run.replace('sim_', ''))}")


def main():
    args = sys.argv[1:]
    chrom = None
    if "--chrom" in args:
        i = args.index("--chrom")
        chrom = args[i + 1]
        del args[i:i + 2]
    out = args[0] if args else "out/chr22"
    written = []
    p = plot_decay(out)
    if p:
        written.append(p)
    for t in sorted(glob.glob(os.path.join(out, "*.map.true.tsv"))):
        written.append(plot_maps(t, chrom))
    if not written:
        sys.exit(f"nothing to plot in {out} (no *.ld_decay.tsv or *.map.true.tsv)")
    for w in written:
        print("written", w)


if __name__ == "__main__":
    main()
