#!/usr/bin/env python3
"""
Plot the actual GRASS-format table's low-density region (e, p, h, Gamma1
vs n0), comparing the current Gamma=1.58425 polytropic tail against the
superseded Gamma~0.9495 tail (proxied by DD2_hot_equal_pinned_fine.dat, a
resample taken before the Gamma fix). Unlike plot_appendixA_extended.py,
this reads the actual GRASS table, so it includes the synthetic cold-
polytrope extension that the Appendix-A pinned-columns diagnostic omits.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPTDIR = os.path.dirname(os.path.realpath(__file__))
NEW_PATH = os.path.join(SCRIPTDIR, "DD2_hot_equal_pinned.dat")
OLD_PATH = os.path.join(SCRIPTDIR, "DD2_hot_equal_pinned_fine.dat")
OUT_PNG = os.path.join(SCRIPTDIR, "DD2_lowdensity_gamma_comparison.png")


def load_grass(path):
    with open(path) as f:
        n = int(f.readline())
    d = np.loadtxt(path, skiprows=1)
    e, p, h, n0 = d[:, 0], d[:, 1], d[:, 2], d[:, 3]
    return {"e": e, "p": p, "h": h, "n0": n0}


def gamma1(e, p):
    le, lp = np.log(e), np.log(p)
    return np.gradient(lp, le)


def main():
    new = load_grass(NEW_PATH)
    old = load_grass(OLD_PATH)

    # native CompOSE floor -- the stitch point for both tables
    n0_edge_new = new["n0"][np.searchsorted(new["e"], 1653.0)]

    panels = [
        ("e", "energy density  e  [g/cm$^3$]", "log"),
        ("p", "pressure  p  [dyn/cm$^2$]", "log"),
        ("h", "GRASS enthalpy  h", "linear"),
        ("gamma1", r"$\Gamma_1 = d\ln P/d\ln e$", "linear"),
    ]

    fig, axes = plt.subplots(len(panels), 1, figsize=(8, 11), sharex=True)

    g_new = gamma1(new["e"], new["p"])
    g_old = gamma1(old["e"], old["p"])

    for ax, (key, label, yscale) in zip(axes, panels):
        if key == "gamma1":
            ax.plot(new["n0"], g_new, "-", color="tab:blue", lw=1.6,
                     label=r"current ($\Gamma=1.58425$ tail)")
            ax.plot(old["n0"], g_old, "--", color="tab:gray", lw=1.4,
                     label=r"superseded ($\Gamma\approx0.9495$ tail)")
            ax.set_ylim(0, 2.0)
        else:
            ax.plot(new["n0"], new[key], "-", color="tab:blue", lw=1.6,
                     label=r"current ($\Gamma=1.58425$ tail)")
            ax.plot(old["n0"], old[key], "--", color="tab:gray", lw=1.4,
                     label=r"superseded ($\Gamma\approx0.9495$ tail)")
        ax.axvline(n0_edge_new, color="0.3", ls=":", lw=1.2)
        ax.set_xscale("log")
        if yscale == "log":
            ax.set_yscale("log")
        ax.set_ylabel(label)
        ax.grid(True, which="both", alpha=0.3)

    axes[0].annotate("DD2's own native\nCompOSE floor\n(stitch point)",
                      xy=(n0_edge_new, new["e"][np.searchsorted(new["n0"], n0_edge_new)]),
                      xytext=(0.02, 0.75), textcoords="axes fraction",
                      fontsize=8, color="0.3",
                      arrowprops=dict(arrowstyle="->", color="0.3", lw=0.8))
    axes[0].legend(loc="lower right", fontsize=9)
    axes[-1].axhline(1.58425, color="tab:blue", ls=":", lw=0.8, alpha=0.6)
    axes[-1].axhline(0.9495, color="tab:gray", ls=":", lw=0.8, alpha=0.6)

    fig.suptitle(
        "DD2_hot_equal_pinned: low-density polytropic tail, before vs after\n"
        "vertical line = stitch point (DD2's own native table floor); "
        "left of it is pure synthetic polytrope",
        fontsize=11,
    )
    axes[-1].set_xlabel(r"$n_0$  [cm$^{-3}$]  (GRASS table's own column, CGS convention)")

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(OUT_PNG, dpi=150)
    print(f"Wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
