"""
Global physical fit of the Yb:YAG RC2 Mueller-matrix data.

What makes this a *physical* fit rather than the earlier point-by-point inversion:
n and k are not free at every wavelength.  They come from one KK-consistent
oscillator model with ~8-17 parameters fitted simultaneously to every wavelength
and every angle of incidence at once.  That is what breaks the absorption/scatter
degeneracy the point-by-point fit could not touch:

  * the 10-AOI reflection Mueller matrix fixes n(lambda) absolutely, with no
    reference to the slab thickness;
  * the normal-incidence transmittance fixes the total per-pass loss (alpha_abs
    + alpha_sc) * d;
  * KK consistency then forbids the absorptive part from taking any shape it
    likes -- a k that rises towards the NIR must show up in n -- so the flat,
    non-KK remainder is pinned as scattering.

Thickness caveat, unchanged and unavoidable: a fringe-free thick slab measures
only the product alpha*d.  d is an input (--thickness, default 1 mm); n is
independent of it, while k, the oscillator amplitudes and the scatter
coefficients all scale as 1/d.  The reported alpha*d columns are the
thickness-independent quantities.

Usage
-----
    python -m fit.run_fit --sample B --thickness 1.0
    python -m fit.run_fit --sample A          # after B: reuses B's dispersion
    python -m fit.run_fit --list              # incl. the screened-out files
"""

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fit.mm_io import read_dat, grid
from fit.slab import ISO_ELEMS, SYM_PAIRS, reflection_mm, transmittance
from fit import dispersion as dsp

# The .dat exports are not tracked in git (they are large binaries); point at the
# measurement folder.  Note the trailing space in the directory name.
DEFAULT_DATA = "/home/aidas/Desktop/FTMC/YbYAG/YbYag txt "

# Sample -> measurement files.  `refl` files were verified to be reflection by
# the absence of TransmissionMeas in the header, NOT by filename: YbYag_5a_MMt /
# YbYag_5a_UMMt are named "...t" but are reflection scans at AOI 30-75 deg.
#
# Data quality was screened before choosing these (see METHOD.md).  Both A and B
# show the Yb(3+) 940 nm pump band at alpha*d ~ 0.54 against a ~0.36 baseline, so
# both are Yb-doped: "pure" in this repository means *uncoated*, not undoped.
# Every sample therefore carries the Yb oscillators.
#
# The 5% and 10% acquisitions are unusable; see UNUSABLE below for the evidence.
SAMPLES = {
    "B": dict(
        refl="YbYag_B_UMMr [2026-06-30,114116].dat",
        trans="YbYag_B_UMMt [2026-06-18,132244].dat",
        note="uncoated Yb:YAG -- full R(10 AOI) + T fit, the reference case",
    ),
    "A": dict(
        refl=None,
        trans="YbYag_A_UMMt  [2026-06-18,133022].dat",
        note="uncoated Yb:YAG, transmission only -> loss fit at fixed dispersion",
    ),
}

# Screened out, kept here so the reason is on the record rather than rediscovered.
UNUSABLE = {
    "10pct": "both YbYag_10%_UMMt files have negative and >1 transmittance; no reflection scan exists",
    "5pct": "no usable data at all -- transmission files corrupt, and BOTH 5a reflection "
            "scans read mm12 ~ 0 / mm33 ~ -1 (Psi=45 deg, Delta=180 deg, i.e. no "
            "polarisation change) at AOI 30-50 deg, the signature of a failed alignment. "
            "Re-measure before this sample can be fitted.",
    "YbYag_5%_UMMt [2026-06-23,123319].dat": "all-NaN",
    "YbYag_5%_UMMt_A [2026-06-23,123646].dat": "all-NaN",
    "YbYag_B [2026-06-23,123926].dat": "all-NaN",
    "YbYag_B_UMMt [2026-06-16,175112].dat": "median T=0.835, at the lossless-slab limit -- looks like a straight-through baseline, not a sample scan",
}

# Instrumental artefacts to exclude.  651-661 nm is the grating turret change;
# 965-980 nm is the Si/InGaAs detector crossover, which unfortunately sits on the
# Yb zero-phonon line -- it is masked by default and can be kept with --keep-zpl.
GRATING_MASK = (651.0, 661.0)
DETECTOR_MASK = (965.0, 980.0)


def _mask_ranges(wls, ranges):
    bad = np.zeros(len(wls), bool)
    for lo, hi in ranges:
        bad |= (wls >= lo) & (wls <= hi)
    return bad


def load_sample(name, data_dir, wl_range, keep_zpl=False):
    """Load and grid one sample's reflection + transmission data."""
    cfg = SAMPLES[name]
    out = {"name": name, "cfg": cfg}

    if cfg["trans"]:
        tdf, tmeta = read_dat(os.path.join(data_dir, cfg["trans"]))
        if not tmeta.get("is_transmission"):
            raise ValueError(f"{cfg['trans']} is not a transmission measurement")
    else:
        tdf, tmeta = None, {}

    if cfg["refl"]:
        rdf, rmeta = read_dat(os.path.join(data_dir, cfg["refl"]))
        if rmeta.get("is_transmission"):
            raise ValueError(f"{cfg['refl']} is a transmission measurement")
        wls = np.array(sorted(rdf.wl.unique()))
        angles = np.array(sorted(rdf.angle.unique()))
    else:
        rdf, rmeta = None, {}
        wls = np.array(sorted(tdf.wl.unique()))
        angles = np.array([])

    keep = (wls >= wl_range[0]) & (wls <= wl_range[1])
    ranges = [GRATING_MASK] + ([] if keep_zpl else [DETECTOR_MASK])
    keep &= ~_mask_ranges(wls, ranges)
    wls = wls[keep]

    out["wls"], out["angles"] = wls, angles
    if rdf is not None:
        raw = {n: grid(rdf, n, wls, angles) for n, _, _ in ISO_ELEMS}
        rawerr = {n: grid(rdf, n + "_err", wls, angles) for n, _, _ in ISO_ELEMS}
        out["raw"], out["depol"] = raw, grid(rdf, "depol", wls, angles)
        # symmetric = signal, antisymmetric = systematic (model predicts 0)
        S, A, E = {}, {}, {}
        for key, a, b, sgn, _, _ in SYM_PAIRS:
            S[key] = 0.5 * (raw[a] + sgn * raw[b])
            A[key] = 0.5 * (raw[a] - sgn * raw[b])
            rep = 0.5 * np.hypot(np.nan_to_num(rawerr[a], nan=1e-4),
                                 np.nan_to_num(rawerr[b], nan=1e-4))
            floor = np.sqrt(np.nanmean(A[key] ** 2))  # instrument systematic
            E[key] = np.hypot(rep, max(floor, 1e-4))
        out["S"], out["A"], out["Serr"] = S, A, E
    else:
        out["S"] = out["Serr"] = None

    if tdf is not None:
        t0 = tdf[tdf.angle == tdf.angle.min()].sort_values("wl")
        out["T"] = np.interp(wls, t0.wl, t0.mm11)
        terr = t0.mm11_err if "mm11_err" in t0 else pd.Series(np.full(len(t0), 1e-3))
        out["Terr"] = np.interp(wls, t0.wl, terr.fillna(1e-3))
        # a physical slab cannot transmit more than the lossless Fresnel limit,
        # nor less than zero -- anything outside is a failed acquisition
        nref = dsp.yag_reference_n(wls)
        Rf = ((1 - nref) / (1 + nref)) ** 2
        out["T"] = np.where((out["T"] > 0) & (out["T"] < (1 - Rf) / (1 + Rf) * 1.02),
                            out["T"], np.nan)
    else:
        out["T"] = out["Terr"] = None
    out["meta"] = {"refl": rmeta, "trans": tmeta}
    return out


def build_model(uv_edge=False, with_yb=True, fit_dispersion=True, with_zpl=False):
    """Seed model + list of free parameter names.

    `fit_dispersion=False` freezes n (used for the transmission-only sample, where
    reflection is absent and n would otherwise be degenerate with the loss).
    """
    m = dsp.yag_seed(with_yb=False)
    m["scatter"] = [3.6, 0.5]          # ~0.36 per-pass optical depth at 1 um
    free = ["scatter_C", "scatter_p"]
    if fit_dispersion:
        free = ["pole0_amp", "pole0_Ec", "pole1_amp"] + free
    if uv_edge:
        m["gaussians"].append([1.0, 6.5, 1.0])
        if fit_dispersion:
            free += ["gauss0_amp", "gauss0_En", "gauss0_Br"]
    if with_yb:
        i0 = len(m["gaussians"])
        # Yb3+ 2F7/2 -> 2F5/2: pump band, emission wing, and (only when the
        # detector-crossover mask is lifted) the 969 nm zero-phonon line.
        # amp ~5e-6 gives the ~0.18 excess per-pass optical depth measured at the
        # 940 nm band for a 1 mm slab.
        # Line list read off the empirical alpha*d spectrum (T inverted at known
        # n): a main band at 941 nm with FWHM ~18 nm (0.025 eV), a resolved
        # shoulder near 915 nm, and the 1030 nm emission wing.  A single broad
        # Gaussian straddling 915+941 reproduces the dip ~40% too shallow.
        lines = [(915.0, 0.020), (941.0, 0.025), (1030.0, 0.030)]
        if with_zpl:
            lines.insert(2, (969.0, 0.008))
        for cen, br in lines:
            m["gaussians"].append([5e-6, dsp.HC / cen, br])
        # The whole Yb manifold lies within 900-1100 nm; hold the oscillators
        # there so none of them can degenerate into a scattering surrogate.
        for j in range(i0, len(m["gaussians"])):
            free += [f"gauss{j}_amp", f"gauss{j}_En", f"gauss{j}_Br"]
            m["bounds"][f"gauss{j}_En"] = (dsp.HC / 1100.0, dsp.HC / 900.0)
            m["bounds"][f"gauss{j}_Br"] = (2e-3, 0.08)
            m["bounds"][f"gauss{j}_amp"] = (0.0, 1e-3)
    return m, free


def make_residuals(data, model, free, d_nm, d_rough, w_trans, t_sys=0.005):
    """Return (resid_fn, x0, bounds, names).  Blocks are balanced by 1/sqrt(N)."""
    wls = data["wls"]
    x0, bounds, names = dsp.pack(model, free)

    has_R = data["S"] is not None
    if has_R:
        rmask = {k: np.isfinite(data["S"][k]) for k, *_ in SYM_PAIRS}
        rerr = {k: np.fmax(data["Serr"][k], 1e-4) for k, *_ in SYM_PAIRS}
        nR = sum(int(m.sum()) for m in rmask.values())
    # The RC2 reports a transmittance sigma of ~6e-4 (0.06%), far below the real
    # photometric accuracy of the measurement.  Left at face value it makes the
    # T block outweigh the whole 10-AOI reflection data set by ~100x and n is
    # sacrificed to chase a systematic.  t_sys is the honest floor.
    has_T = data["T"] is not None
    if has_T:
        tmask = np.isfinite(data["T"]) & (data["T"] > 1e-3)
        terr = np.hypot(np.nan_to_num(data["Terr"], nan=1e-3), t_sys)
        nT = int(tmask.sum())
    else:
        nT = 0

    wR = 1.0 / np.sqrt(nR) if has_R else 0.0
    wT = w_trans / np.sqrt(max(nT, 1))

    def resid(x):
        m = dsp.unpack(model, names, x)
        N = dsp.nk(m, wls)
        asc = dsp.scatter_alpha(m, wls)
        parts = []
        if has_R:
            MM = reflection_mm(N, wls, data["angles"], d_nm, asc, d_rough)
            for key, _a, _b, _s, i, j in SYM_PAIRS:
                r = (MM[:, :, i, j] - data["S"][key]) / rerr[key]
                parts.append(wR * r[rmask[key]])
        if has_T:
            Tm = transmittance(N, wls, [0.0], d_nm, asc, d_rough)[:, 0]
            parts.append(wT * ((Tm - data["T"]) / terr)[tmask])
        return np.concatenate(parts)

    return resid, x0, bounds, names


def report(data, m, d_nm, d_rough, out_dir):
    """Write fitted optical constants and a comparison against reference YAG."""
    wls = data["wls"]
    N = dsp.nk(m, wls)
    asc = dsp.scatter_alpha(m, wls)
    alpha_abs = 4 * np.pi * N.imag / wls          # 1/nm
    df = pd.DataFrame({
        "wl_nm": wls,
        "n": N.real,
        "k": N.imag,
        "alpha_abs_percm": alpha_abs * 1e7,
        "alpha_scat_percm": asc * 1e7,
        "alphad_abs": alpha_abs * d_nm,           # thickness-independent
        "alphad_scat": asc * d_nm,
        "yag_ref_n": dsp.yag_reference_n(wls),
    })
    df["T_model"] = transmittance(N, wls, [0.0], d_nm, asc, d_rough)[:, 0]
    df["T_meas"] = data["T"] if data["T"] is not None else np.nan
    if data["S"] is not None:
        MM = reflection_mm(N, wls, data["angles"], d_nm, asc, d_rough)
        ia = int(np.argmin(np.abs(data["angles"] - 60.0)))  # representative AOI
        for key, _a, _b, _s, i, j in SYM_PAIRS:
            df[f"{key}_model"] = MM[:, ia, i, j]
            df[f"{key}_meas"] = data["S"][key][:, ia]
    path = os.path.join(out_dir, f"fit_{data['name']}_nk.csv")
    df.to_csv(path, index=False)
    return df, path


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample", default="B", choices=sorted(SAMPLES))
    p.add_argument("--data-dir", default=DEFAULT_DATA)
    p.add_argument("--thickness", type=float, default=1.0, help="slab thickness, mm")
    p.add_argument("--roughness", type=float, default=0.0, help="surface roughness, nm")
    p.add_argument("--wl-min", type=float, default=260.0)
    p.add_argument("--wl-max", type=float, default=1690.0)
    p.add_argument("--uv-edge", action="store_true",
                   help="add a band-edge Gaussian (use with --wl-min 210)")
    p.add_argument("--keep-zpl", action="store_true",
                   help="keep 965-980 nm (Yb zero-phonon line sits on the detector crossover)")
    p.add_argument("--t-sys", type=float, default=0.005,
                   help="systematic floor on transmittance (default 0.5%%)")
    p.add_argument("--w-trans", type=float, default=1.0,
                   help="relative weight of the transmission block")
    p.add_argument("--out", default="fit/out")
    p.add_argument("--list", action="store_true")
    a = p.parse_args(argv)

    if a.list:
        for k, v in SAMPLES.items():
            print(f"{k:7s} refl={str(v['refl']):48s} trans={v['trans']}  # {v['note']}")
        print("\nscreened out:")
        for k, v in UNUSABLE.items():
            print(f"  {k}: {v}")
        return 0

    os.makedirs(a.out, exist_ok=True)
    d_nm = a.thickness * 1e6

    data = load_sample(a.sample, a.data_dir, (a.wl_min, a.wl_max), a.keep_zpl)
    has_R = data["S"] is not None
    model, free = build_model(uv_edge=a.uv_edge, with_yb=True,
                              fit_dispersion=has_R, with_zpl=a.keep_zpl)
    if not has_R:
        # transmission-only: n is degenerate with the loss, so freeze the
        # dispersion at the reflection-fitted result from sample B
        ref = os.path.join(a.out, "fit_B_model.json")
        if os.path.exists(ref):
            base = json.load(open(ref))["model"]
            model["poles"] = [list(q) for q in base["poles"]]
            model["eps_inf"] = base["eps_inf"]
            print(f"[n frozen at the sample-B reflection fit: {ref}]")
        else:
            print("[n frozen at the Zelmon YAG reference; run --sample B first"
                  " to use the measured dispersion]")
    resid, x0, bounds, names = make_residuals(
        data, model, free, d_nm, a.roughness, a.w_trans, a.t_sys)

    print(f"=== sample {a.sample} ({SAMPLES[a.sample]['note']}) ===")
    print(f"wavelengths {len(data['wls'])} in [{data['wls'].min():.0f}, "
          f"{data['wls'].max():.0f}] nm; AOI {list(data['angles'].astype(int))}; "
          f"transmission: {'yes' if data['T'] is not None else 'no (files unusable)'}")
    print(f"thickness {a.thickness} mm (k, amplitudes and scatter scale as 1/d)")
    print(f"free parameters ({len(names)}): {', '.join(names)}")
    r0 = resid(x0)
    print(f"residuals: {r0.size}, start cost {0.5*np.sum(r0**2):.5g}")

    sol = least_squares(resid, x0, bounds=bounds, method="trf",
                        x_scale="jac", max_nfev=4000, verbose=0)
    m = dsp.unpack(model, names, sol.x)
    print(f"converged: {sol.success} ({sol.message.strip()})")
    print(f"final cost {sol.cost:.5g}, nfev={sol.nfev}")

    # 1-sigma parameter uncertainties from the Jacobian at the solution, scaled
    # by the residual variance (the weights are relative, not absolute sigmas).
    try:
        J = sol.jac
        dof = max(J.shape[0] - J.shape[1], 1)
        cov = np.linalg.inv(J.T @ J) * (2 * sol.cost / dof)
        perr = np.sqrt(np.abs(np.diag(cov)))
    except np.linalg.LinAlgError:
        perr = np.full(len(names), np.nan)
    print("\nfitted parameters (1-sigma):")
    for nm_, v, e in zip(names, sol.x, perr):
        print(f"  {nm_:14s} {v:12.6g} +/- {e:.3g}")

    df, path = report(data, m, d_nm, a.roughness, a.out)
    with open(os.path.join(a.out, f"fit_{a.sample}_model.json"), "w") as fh:
        json.dump({"sample": a.sample, "thickness_mm": a.thickness,
                   "roughness_nm": a.roughness, "model": m,
                   "free": names, "cost": float(sol.cost),
                   "stderr": dict(zip(names, map(float, perr)))}, fh, indent=2)

    good = np.isfinite(df.n) & (df.wl_nm > 350) & (df.wl_nm < 1650)
    dn = df.n[good] - df.yag_ref_n[good]
    print(f"\nn vs Zelmon single-crystal YAG: RMS {np.sqrt(np.mean(dn**2)):.5f}, "
          f"max |dev| {np.max(np.abs(dn)):.5f}")
    print("\n  lambda      n      YAG ref    alpha*d(abs)  alpha*d(scat)")
    for w in (400, 500, 633, 800, 940, 1000, 1064, 1300, 1550):
        if w < df.wl_nm.min() or w > df.wl_nm.max():
            continue
        i = int(np.argmin(np.abs(df.wl_nm.values - w)))
        print(f"  {w:5d} nm  {df.n[i]:.4f}   {df.yag_ref_n[i]:.4f}     "
              f"{df.alphad_abs[i]:.4f}        {df.alphad_scat[i]:.4f}")
    if data["T"] is not None:
        resT = (df.T_model - df.T_meas).values
        print(f"\ntransmittance residual: mean {np.nanmean(resT):+.4f}, "
              f"RMS {np.sqrt(np.nanmean(resT**2)):.4f}")
    print("\nscatter: alpha_sc = %.3f (1000nm/lam)^%.2f 1/cm" % tuple(m["scatter"]))
    for g in m["gaussians"]:
        print(f"  Gaussian: A={g[0]:.4g}  E0={g[1]:.4f} eV ({dsp.HC/g[1]:.1f} nm)"
              f"  Br={g[2]:.4f} eV")
    print(f"\nwrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
