"""
Common plot style of the TB2026 notes, as in the final-sample report (final_panels.py / final_profiles.py):

  data                     blue solid points with their error bars (and the bin width as a horizontal bar)
  simulation, no digi      orange, dotted histogram
  simulation, digitised    red, hatched histogram; several digitised models (geometry variants, campaigns) take a
                           gradient of reds, light to dark, each with its own hatch and line style

Binned quantities (per layer, per radial bin, per time bin) are drawn as histograms; ratios and trends with several
curves (against energy, k, Q, distance) draw the simulations as lines and the data as points with error bars.
"""
import numpy as np

DATA = "#29a3dc"
NODIGI = "#eda100"
DIGI = "#c0392b"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
# Digitised models, light to dark around DIGI; hatches and line styles cycle with them.
REDS = ["#f08c76", "#e0583c", "#c0392b", "#8e2116", "#5e1008"]
HATCHES = ["////", "\\\\\\\\", "xxx", "...", "++"]
LINESTYLES = ["-", "--", "-.", ":", (0, (5, 1, 1, 1))]
MARKERS = ["o", "s", "^", "D", "v"]

RC = {"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
      "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
      "axes.spines.right": False, "legend.frameon": False, "lines.linewidth": 1.8, "hatch.linewidth": 0.8}


def reds(n):
    """n digitised-model styles, light to dark: dicts with color, hatch, ls, marker."""
    idx = np.linspace(0, len(REDS) - 1, n).round().astype(int) if n > 1 else [2]
    return [{"color": REDS[i], "hatch": HATCHES[k % len(HATCHES)], "ls": LINESTYLES[k % len(LINESTYLES)],
             "marker": MARKERS[k % len(MARKERS)]} for k, i in enumerate(idx)]


def edges_of(x):
    """Bin edges around bin centres x (halfway between neighbours)."""
    x = np.asarray(x, float)
    mid = 0.5 * (x[1:] + x[:-1])
    return np.r_[x[0] - (mid[0] - x[0]), mid, x[-1] + (x[-1] - mid[-1])]


def data_points(ax, x, y, yerr=None, edges=None, label="data", color=DATA, **kw):
    """Data: solid points with error bars; with edges, a horizontal bar across the bin."""
    xerr = kw.pop("xerr", None) if edges is None else [x - edges[:-1], edges[1:] - x]
    ax.errorbar(x, y, yerr=yerr, xerr=xerr, fmt=kw.pop("fmt", "o"), color=color, ms=kw.pop("ms", 5.5), mew=0,
                elinewidth=1.3, capsize=0, label=label, zorder=5, **kw)


def sim_hist(ax, values, edges, style, label=None, fill=True):
    """A digitised (or any) simulation as a histogram: hatched fill in its colour plus the outline."""
    if fill:
        ax.stairs(values, edges, fill=True, facecolor="none", edgecolor=style["color"], hatch=style["hatch"], lw=0,
                  alpha=0.55)
    ax.stairs(values, edges, color=style["color"], lw=1.7, ls=style["ls"], label=label)


def nodigi_hist(ax, values, edges, label=None):
    ax.stairs(values, edges, color=NODIGI, lw=2.0, ls=":", label=label)


def sim_line(ax, x, y, style, label=None, markers=True, **kw):
    """A simulation in a ratio / trend comparison: a line in its colour and style."""
    ax.plot(x, y, color=style["color"], ls=style["ls"], marker=style["marker"] if markers else None, ms=4,
            label=label, **kw)


PAD_MM = 5.53   # SiW-ECAL pad pitch


def pad_band(ax, label=True):
    """Shade +-half a pad around 0 on a sim - data width axis. Not an uncertainty (R90 over thousands of events is
    known far better than a pad): it marks the detector's granularity, differences inside it are below one cell."""
    ax.axhspan(-PAD_MM / 2, PAD_MM / 2, color=GRID, alpha=0.6, lw=0, zorder=0,
               label=f"$\\pm$ half a pad ({PAD_MM / 2:.2f} mm)" if label else None)
