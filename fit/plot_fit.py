"""
Figures for the Yb:YAG physical fit.

    python -m fit.plot_fit --sample B
"""

import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fit.run_fit import load_sample, DEFAULT_DATA, SAMPLES
from fit.slab import SYM_PAIRS, reflection_mm, transmittance
from fit import dispersion as dsp

LABEL = {"s12": r"$(m_{12}+m_{21})/2$", "s33": r"$(m_{33}+m_{44})/2$",
         "s34": r"$(m_{34}-m_{43})/2$"}


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--sample", default="B", choices=sorted(SAMPLES))
    p.add_argument("--data-dir", default=DEFAULT_DATA)
    p.add_argument("--out", default="fit/out")
    a = p.parse_args(argv)

    js = json.load(open(os.path.join(a.out, f"fit_{a.sample}_model.json")))
    m, d_nm = js["model"], js["thickness_mm"] * 1e6
    data = load_sample(a.sample, a.data_dir, (260.0, 1690.0))
    wls = data["wls"]
    N = dsp.nk(m, wls)
    asc = dsp.scatter_alpha(m, wls)

    # ---- fig 1: n, and the loss budget -----------------------------------
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    ax[0].plot(wls, N.real, lw=1.6, label="fit")
    ax[0].plot(wls, dsp.yag_reference_n(wls), "--", lw=1.2,
               label="Zelmon single-crystal YAG")
    ax[0].set_xlabel("wavelength (nm)")
    ax[0].set_ylabel("n")
    ax[0].legend(frameon=False)
    ax[0].set_title(f"{a.sample}: refractive index")

    ax[1].plot(wls, N.real - dsp.yag_reference_n(wls), lw=1.2, color="C3")
    ax[1].axhline(0, color="k", lw=0.6)
    ax[1].set_xlabel("wavelength (nm)")
    ax[1].set_ylabel(r"$n_{\rm fit}-n_{\rm YAG}$")
    ax[1].set_title("residual vs reference YAG")

    ad_abs = 4 * np.pi * N.imag / wls * d_nm
    ad_sc = asc * d_nm
    ax[2].plot(wls, ad_abs + ad_sc, lw=1.6, label="total")
    ax[2].plot(wls, ad_sc, lw=1.2, label=r"scattering (non-KK)")
    ax[2].plot(wls, ad_abs, lw=1.2, label=r"Yb$^{3+}$ absorption (KK)")
    ax[2].set_xlabel("wavelength (nm)")
    ax[2].set_ylabel(r"per-pass optical depth $\alpha d$")
    ax[2].legend(frameon=False)
    ax[2].set_title("loss budget (thickness-independent)")
    fig.tight_layout()
    f1 = os.path.join(a.out, f"fig_{a.sample}_nk.png")
    fig.savefig(f1, dpi=140)

    # ---- fig 2: measured vs modelled observables -------------------------
    ncol = 1 + (3 if data["S"] is not None else 0)
    fig, ax = plt.subplots(1, ncol, figsize=(4.6 * ncol, 4.0), squeeze=False)
    ax = ax[0]
    if data["S"] is not None:
        MM = reflection_mm(N, wls, data["angles"], d_nm, asc)
        show = [i for i, v in enumerate(data["angles"]) if v in (45.0, 60.0, 70.0)]
        for c, (key, _a, _b, _s, i, j) in enumerate(SYM_PAIRS):
            for kk, ia in enumerate(show):
                ax[c].plot(wls, data["S"][key][:, ia], lw=0.8, color=f"C{kk}",
                           label=f"{data['angles'][ia]:.0f}$^\\circ$ meas")
                ax[c].plot(wls, MM[:, ia, i, j], "--", lw=1.0, color="k")
            ax[c].set_xlabel("wavelength (nm)")
            ax[c].set_ylabel(LABEL[key])
            ax[c].legend(frameon=False, fontsize=8)
            ax[c].set_title(f"{LABEL[key]} (dashed = model)")
    if data["T"] is not None:
        Tm = transmittance(N, wls, [0.0], d_nm, asc)[:, 0]
        ax[-1].plot(wls, data["T"], lw=0.9, label="measured")
        ax[-1].plot(wls, Tm, "--", lw=1.1, color="k", label="model")
        ax[-1].set_xlabel("wavelength (nm)")
        ax[-1].set_ylabel("transmittance")
        ax[-1].legend(frameon=False)
        ax[-1].set_title("normal-incidence transmittance")
    fig.tight_layout()
    f2 = os.path.join(a.out, f"fig_{a.sample}_fit_vs_data.png")
    fig.savefig(f2, dpi=140)
    print(f"wrote {f1}\nwrote {f2}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
