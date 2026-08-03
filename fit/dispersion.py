"""
Kramers-Kronig-consistent dielectric-function model for (Yb:)YAG.

Everything that lives in eps is KK-consistent by construction:

  eps(E) = eps_inf
         + sum_p  A_p / (Ec_p^2 - E^2)             Sellmeier poles (lossless below Ec)
         + sum_g  Gaussian(A_g, En_g, Br_g)        KK pair via the Dawson function

The Gaussian is the WVASE/CompleteEASE "Gaussian" oscillator -- the right shape for
inhomogeneously broadened Yb(3+) 4f-4f crystal-field lines (a Lorentzian tail would
leak absorption across the whole transparent window):

  sigma  = Br / (2 sqrt(ln 2))
  eps2(E) = A [ exp(-((E-En)/sigma)^2) - exp(-((E+En)/sigma)^2) ]
  eps1(E) = (2A/sqrt(pi)) [ D((E+En)/sigma) - D((E-En)/sigma) ]      D = Dawson F

The second (negative-energy) term enforces the odd symmetry of eps2 that the KK
integral demands; dropping it biases eps1 at low E.

SCATTERING IS DELIBERATELY *NOT* IN eps.  A ceramic/defect haze removes light from
the specular beam by redirecting it, not by absorbing it, so it is not KK-consistent
and forcing it into k distorts both n and the Yb oscillator strengths.  It is carried
separately as an intensity loss coefficient (see `scatter_alpha`) applied to the
propagation through the slab only.
"""

import re

import numpy as np
from scipy.special import dawsn

HC = 1239.841984  # eV*nm


def ev(lam_nm):
    """Photon energy (eV) from vacuum wavelength (nm)."""
    return HC / np.asarray(lam_nm, dtype=float)


def sellmeier_pole(E, amp, Ec):
    """Lossless pole: eps1 = amp / (Ec^2 - E^2).  amp in eV^2."""
    return amp / (Ec ** 2 - E ** 2)


def gaussian_osc(E, amp, En, Br):
    """KK-consistent Gaussian oscillator, returns complex eps contribution."""
    sig = Br / (2.0 * np.sqrt(np.log(2.0)))
    e2 = amp * (np.exp(-(((E - En) / sig) ** 2)) - np.exp(-(((E + En) / sig) ** 2)))
    e1 = (2.0 * amp / np.sqrt(np.pi)) * (dawsn((E + En) / sig) - dawsn((E - En) / sig))
    return e1 + 1j * e2


# --------------------------------------------------------------------------
# Parameter handling
#
# A model is a plain dict:
#   {"eps_inf": float,
#    "poles":     [(amp, Ec), ...],
#    "gaussians": [(amp, En, Br), ...],
#    "scatter":   (C, p)}   alpha_sc = C (1000nm/lam)^p
# `PARAM_SPEC` flattens it to the vector scipy fits, with bounds.
# --------------------------------------------------------------------------

def pack(model, free):
    """Flatten the free parameters of `model` into a vector + bounds."""
    x, lo, hi, names = [], [], [], []
    for name, val, b0, b1 in _iter_params(model):
        if name in free:
            x.append(val)
            lo.append(b0)
            hi.append(b1)
            names.append(name)
    return np.array(x), (np.array(lo), np.array(hi)), names


def unpack(model, names, x):
    """Return a copy of `model` with the named parameters replaced by `x`."""
    m = {
        "eps_inf": model["eps_inf"],
        "poles": [list(p) for p in model["poles"]],
        "gaussians": [list(g) for g in model["gaussians"]],
        "scatter": list(model["scatter"]),
        "bounds": dict(model.get("bounds", {})),
        "yb_assign": list(model.get("yb_assign", [])),
    }
    for name, v in zip(names, x):
        if name == "eps_inf":
            m["eps_inf"] = v
        elif name.startswith("pole"):
            i, fld = _idx(name)
            m["poles"][i][{"amp": 0, "Ec": 1}[fld]] = v
        elif name.startswith("gauss"):
            i, fld = _idx(name)
            m["gaussians"][i][{"amp": 0, "En": 1, "Br": 2}[fld]] = v
        elif name == "scatter_C":
            m["scatter"][0] = v
        elif name == "scatter_p":
            m["scatter"][1] = v
        else:
            raise KeyError(name)
    return m


def _idx(name):
    """Split e.g. 'gauss12_amp' into (12, 'amp').

    The index must be parsed as a whole number, not as the trailing character:
    with ten or more oscillators `int(head[-1])` maps gauss10..gauss13 onto
    0..3, silently aliasing four parameters onto four others.  That makes the
    Jacobian singular (every uncertainty comes back NaN) and lets amplitudes be
    driven to absurd values, while the reported per-line table -- which reads
    model["gaussians"] directly -- still looks plausible.
    """
    m = re.fullmatch(r"([A-Za-z]+)(\d+)_(\w+)", name)
    if m is None:
        raise KeyError(name)
    return int(m.group(2)), m.group(3)


def _iter_params(model):
    """Yield (name, value, lo, hi) for every parameter of the model.

    `model["bounds"]` may override any entry as {name: (lo, hi)}.  This matters
    for the Yb oscillators: left with generic bounds one of them runs off to
    ~1570 nm with a 1.3 eV width and becomes a smooth background absorption --
    a KK-consistent stand-in for scattering, which is exactly the degeneracy this
    model is built to avoid.  Confining the centres to the Yb(3+)
    2F7/2 -> 2F5/2 manifold (900-1100 nm) and the widths to <= 0.08 eV keeps every
    oscillator a real transition.
    """
    ov = model.get("bounds", {})

    def b(name, default):
        return ov.get(name, default)

    yield ("eps_inf", model["eps_inf"], *b("eps_inf", (1.0, 5.0)))
    for i, (amp, Ec) in enumerate(model["poles"]):
        yield (f"pole{i}_amp", amp, *b(f"pole{i}_amp", (0.0, 1e4)))
        yield (f"pole{i}_Ec", Ec, *b(f"pole{i}_Ec", (0.02, 40.0)))
    for i, (amp, En, Br) in enumerate(model["gaussians"]):
        yield (f"gauss{i}_amp", amp, *b(f"gauss{i}_amp", (0.0, 10.0)))
        yield (f"gauss{i}_En", En, *b(f"gauss{i}_En", (0.05, 12.0)))
        yield (f"gauss{i}_Br", Br, *b(f"gauss{i}_Br", (1e-3, 5.0)))
    C, pw = model["scatter"]
    yield ("scatter_C", C, *b("scatter_C", (0.0, 1e3)))
    yield ("scatter_p", pw, *b("scatter_p", (0.0, 6.0)))


def eps(model, lam_nm):
    """Complex dielectric function of the model on a wavelength grid (nm)."""
    E = ev(lam_nm)
    out = np.full(E.shape, complex(model["eps_inf"]), dtype=complex)
    for amp, Ec in model["poles"]:
        out += sellmeier_pole(E, amp, Ec)
    for amp, En, Br in model["gaussians"]:
        out += gaussian_osc(E, amp, En, Br)
    return out


def nk(model, lam_nm):
    """Complex refractive index N = n + i k."""
    return np.sqrt(eps(model, lam_nm))


SCATTER_REF_NM = 1000.0


def scatter_alpha(model, lam_nm):
    """Non-KK haze loss coefficient alpha_sc (1/nm), as a free-exponent power law.

        alpha_sc(lam) = C * (1000 nm / lam)^p          C in 1/cm, converted to 1/nm

    The exponent is fitted rather than fixed at the Rayleigh value.  Rayleigh
    (p = 4) applies only to scatterers much smaller than the wavelength; the
    measured loss here falls by ~20% from 300 to 1650 nm, whereas p = 4 would
    vary by a factor of ~900 across that span, so a fixed lam^-4 term is simply
    driven to zero and leaves the NIR transmittance ~2% under-predicted.  Grain
    boundaries and residual pores in a ceramic are comparable to or larger than
    the wavelength, i.e. the Mie/geometric regime, where p is small (0-2).  The
    fitted p is therefore a diagnostic of scatterer size, not a nuisance knob.
    """
    C, p = model["scatter"]
    lam = np.asarray(lam_nm, dtype=float)
    return C * (SCATTER_REF_NM / lam) ** p * 1e-7  # 1/cm -> 1/nm


# --------------------------------------------------------------------------
# Seed: Zelmon single-crystal YAG (Appl. Opt. 37, 4933 (1998)), converted from
#   n^2 = 1 + 2.28200 l^2/(l^2-0.01185) + 3.27644 l^2/(l^2-282.734)   [l in um]
# to energy-space poles via  B l^2/(l^2-C) = (B Ec^2)/(Ec^2 - E^2),  Ec = hc/sqrt(C).
# --------------------------------------------------------------------------

def yag_seed(with_yb=False):
    uv_Ec = HC / 1000.0 / np.sqrt(0.01185)   # 11.388 eV
    ir_Ec = HC / 1000.0 / np.sqrt(282.734)   # 0.07373 eV
    model = {
        "eps_inf": 1.0,
        "poles": [(2.28200 * uv_Ec ** 2, uv_Ec), (3.27644 * ir_Ec ** 2, ir_Ec)],
        "gaussians": [],
        "scatter": [0.0, 0.0],
        "bounds": {},
    }
    if with_yb:
        # Yb(3+) 2F7/2 -> 2F5/2 manifold: pump band, zero-phonon line, emission wing.
        model["gaussians"] = [
            [2e-4, HC / 940.0, 0.030],
            [2e-4, HC / 969.0, 0.010],
            [1e-4, HC / 1030.0, 0.030],
        ]
    return model


def yag_reference_n(lam_nm):
    """Zelmon single-crystal YAG n(lambda) -- independent check, not fitted."""
    l = np.asarray(lam_nm, dtype=float) / 1000.0
    n2 = 1 + 2.28200 * l ** 2 / (l ** 2 - 0.01185) + 3.27644 * l ** 2 / (l ** 2 - 282.734)
    return np.sqrt(n2)


# --------------------------------------------------------------------------
# Self-test: verify the analytic Gaussian eps1 really is the KK transform of eps2.
# --------------------------------------------------------------------------

def _kk_selftest(verbose=True):
    """Numerically KK-transform eps2 of a Gaussian and compare to the analytic eps1."""
    amp, En, Br = 0.5, 1.32, 0.05
    Egrid = np.linspace(1e-4, 12.0, 400001)
    g = gaussian_osc(Egrid, amp, En, Br)
    e2 = g.imag

    # eps1(E) = (2/pi) P int_0^L x eps2(x)/(x^2 - E^2) dx.
    # The pole at x=E is removed by singularity subtraction, which is exact and
    # stays well behaved *on* resonance (a bare masked grid does not):
    #   P int x f/(x^2-a^2) = int [x f(x) - a f(a)]/(x^2-a^2) + (f(a)/2) ln((L-a)/(L+a))
    # using P int_0^L dx/(x^2-a^2) = (1/2a) ln((L-a)/(L+a)).
    trapz = getattr(np, "trapezoid", None) or np.trapz  # numpy <2 compat
    L = Egrid[-1]
    probes = np.array([0.8, 1.0, 1.2, 1.32, 1.5, 2.0, 3.0])
    err = []
    for E0 in probes:
        f0 = np.interp(E0, Egrid, e2)
        d = Egrid ** 2 - E0 ** 2
        reg = np.where(np.abs(d) > 1e-12, (Egrid * e2 - E0 * f0) / d, 0.0)
        integ = trapz(reg, Egrid) + 0.5 * f0 * np.log((L - E0) / (L + E0))
        num = (2.0 / np.pi) * integ
        ana = np.interp(E0, Egrid, g.real)
        err.append(abs(num - ana))
        if verbose:
            print(f"  E={E0:5.2f} eV  KK-numeric={num:+.6e}  analytic={ana:+.6e}  d={err[-1]:.2e}")
    return max(err)


def _param_roundtrip_selftest():
    """pack/unpack must survive >=10 oscillators (regression: gauss1X aliasing)."""
    m = yag_seed()
    m["gaussians"] = [[1e-5 * (i + 1), 1.0 + 0.01 * i, 0.01] for i in range(14)]
    free = [f"gauss{i}_amp" for i in range(14)]
    x, _, names = pack(m, free)
    x2 = x * 3.0
    m2 = unpack(m, names, x2)
    got = np.array([g[0] for g in m2["gaussians"]])
    want = np.array([1e-5 * (i + 1) * 3.0 for i in range(14)])
    return float(np.max(np.abs(got - want)))


if __name__ == "__main__":
    err = _param_roundtrip_selftest()
    print(f"parameter pack/unpack with 14 oscillators: max error {err:.2e}"
          f"  {'OK' if err < 1e-18 else 'FAIL'}")
    print()
    print("KK consistency of the Gaussian oscillator (numeric vs analytic eps1):")
    worst = _kk_selftest()
    print(f"worst |difference| = {worst:.2e}")

    lam = np.array([400.0, 500.0, 633.0, 1064.0, 1550.0])
    m = yag_seed()
    print("\nZelmon-seed n(lambda) vs reference:")
    for l, a, b in zip(lam, nk(m, lam).real, yag_reference_n(lam)):
        print(f"  {l:7.1f} nm  seed={a:.5f}  Zelmon={b:.5f}")
