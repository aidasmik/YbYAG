"""
Forward optical model: a thick, transparent, double-side-polished slab in air.

This is the physics CompleteEASE applies to a transparent substrate with
"Model Backside Reflections" switched on.  The slab is millimetre-thick and
fringe-free at the RC2's resolution, so the front and back beams add
*incoherently* -- Mueller matrices are summed, not Jones amplitudes.  That sum is
what produces the measured depolarisation and the angle-fanning of the pseudo
index, and it is why a naive single-interface fit reports n too low.

An optional surface-roughness layer (Bruggeman 50/50 void:material) sits on the
front face and is treated *coherently* (it is nm-thick), giving effective front
r/t coefficients that then enter the incoherent slab sum.  d_rough = 0 reduces
exactly to the bare interface.

Sign/phase convention: N = n + i k with exp(-i omega t), so a forward-propagating
wave carries exp(+i 2 pi N cos(theta) d / lambda) and Im(N) > 0 attenuates.

All routines are vectorised over a wavelength axis and an angle axis and return
arrays shaped [n_wl, n_angle].
"""

import numpy as np

# Reflection Mueller elements that are non-zero for an isotropic sample.
# The isotropic reflection Mueller matrix is block diagonal:
#   [[1, -cos2Psi, 0, 0], [-cos2Psi, 1, 0, 0],
#    [0, 0, sin2Psi cosD, sin2Psi sinD], [0, 0, -sin2Psi sinD, sin2Psi cosD]]
ISO_ELEMS = [("mm12", 0, 1), ("mm21", 1, 0), ("mm33", 2, 2),
             ("mm44", 3, 3), ("mm34", 2, 3), ("mm43", 3, 2)]

# Isotropy forces mm12 == mm21, mm33 == mm44 and mm34 == -mm43.  Each pair is
# therefore split into a symmetric combination, which carries all the physical
# signal, and an antisymmetric one, which the model predicts to be exactly zero
# and which consequently measures the instrument's systematic error directly.
# Fitting the symmetric part and taking the error floor from the antisymmetric
# part is what keeps a ~0.7% calibration systematic from biasing n: the RC2
# reports sigma ~1e-4, two orders of magnitude below the real disagreement, and
# ~10% of the measured mm12 values fall (impossibly) below -1.
SYM_PAIRS = [("s12", "mm12", "mm21", +1, 0, 1),
             ("s33", "mm33", "mm44", +1, 2, 2),
             ("s34", "mm34", "mm43", -1, 2, 3)]


def _cos_theta(N0, N1, cos0):
    """cos(theta_1) inside medium N1 from Snell, on the physical branch."""
    sin1 = N0 * np.sqrt(np.clip(1.0 - cos0 ** 2, 0.0, None)) / N1
    c1 = np.sqrt(1.0 - sin1 ** 2 + 0j)
    # keep the decaying branch (Im >= 0) so evanescent/absorbing waves decay
    return np.where(c1.imag < 0, -c1, c1)


def _fresnel(Na, Nb, ca, cb):
    """Amplitude r, t for a -> b, (rs, rp, ts, tp)."""
    rs = (Na * ca - Nb * cb) / (Na * ca + Nb * cb)
    rp = (Nb * ca - Na * cb) / (Nb * ca + Na * cb)
    ts = 2 * Na * ca / (Na * ca + Nb * cb)
    tp = 2 * Na * ca / (Nb * ca + Na * cb)
    return rs, rp, ts, tp


def bruggeman_void(eps_mat, f_void=0.5):
    """Bruggeman EMA of material + void.  Solves f1(e1-e)/(e1+2e)+f2(e2-e)/(e2+2e)=0."""
    e1, e2 = eps_mat, np.ones_like(eps_mat)
    f1, f2 = 1.0 - f_void, f_void
    # quadratic  2e^2 + b e - c = 0  with the standard two-component reduction
    b = (2 * f1 - f2) * e2 + (2 * f2 - f1) * e1
    disc = np.sqrt(b ** 2 + 8 * e1 * e2 + 0j)
    e = (b + disc) / 4.0
    return np.where(e.imag < 0, (b - disc) / 4.0, e)


def _mueller_diag(ap, asx):
    """Mueller matrix of a diagonal Jones matrix diag(ap, as), shape [...,4,4]."""
    Ipp = np.abs(ap) ** 2
    Iss = np.abs(asx) ** 2
    cr = ap * np.conj(asx)
    M = np.zeros(ap.shape + (4, 4), dtype=float)
    M[..., 0, 0] = 0.5 * (Ipp + Iss)
    M[..., 0, 1] = 0.5 * (Ipp - Iss)
    M[..., 1, 0] = M[..., 0, 1]
    M[..., 1, 1] = M[..., 0, 0]
    M[..., 2, 2] = cr.real
    M[..., 2, 3] = cr.imag
    M[..., 3, 2] = -cr.imag
    M[..., 3, 3] = cr.real
    return M


def _front_coeffs(N0, N1, lam, cos0, d_rough):
    """Effective front-face r/t including a coherent roughness layer.

    Returns (r_in, t_in, r_out, t_out, c1) where *_in are for ambient->substrate
    and *_out for substrate->ambient, each a (s, p) tuple.
    """
    c1 = _cos_theta(N0, N1, cos0)
    if d_rough <= 0:
        rs, rp, ts, tp = _fresnel(N0, N1, cos0, c1)
        rs_o, rp_o, ts_o, tp_o = _fresnel(N1, N0, c1, cos0)
        return (rs, rp), (ts, tp), (rs_o, rp_o), (ts_o, tp_o), c1

    Nr = np.sqrt(bruggeman_void(N1 ** 2))
    cr = _cos_theta(N0, Nr, cos0)
    b = 2j * np.pi * Nr * cr * d_rough / lam
    eb = np.exp(b)

    r01 = _fresnel(N0, Nr, cos0, cr)
    r12 = _fresnel(Nr, N1, cr, c1)
    r10 = _fresnel(Nr, N0, cr, cos0)
    r21 = _fresnel(N1, Nr, c1, cr)

    def comp(a, bb, ta, tb):
        den = 1.0 + a * bb * eb ** 2
        return (a + bb * eb ** 2) / den, (ta * tb * eb) / den

    rs_in, ts_in = comp(r01[0], r12[0], r01[2], r12[2])
    rp_in, tp_in = comp(r01[1], r12[1], r01[3], r12[3])
    rs_out, ts_out = comp(r21[0], r10[0], r21[2], r10[2])
    rp_out, tp_out = comp(r21[1], r10[1], r21[3], r10[3])
    return (rs_in, rp_in), (ts_in, tp_in), (rs_out, rp_out), (ts_out, tp_out), c1


def reflection_mm(N, lam, angles_deg, d_nm, alpha_sc=0.0, d_rough=0.0, n_back=3):
    """Normalised reflection Mueller matrix of the slab.

    Parameters
    ----------
    N : complex array [n_wl]      substrate index
    lam : array [n_wl]            wavelength, nm
    angles_deg : array [n_ang]    angles of incidence
    d_nm : float                  slab thickness, nm
    alpha_sc : array [n_wl]       non-KK scatter loss, 1/nm
    n_back : int                  number of backside round trips retained

    Returns array [n_wl, n_ang, 4, 4] normalised so M[...,0,0] == 1.
    """
    N1 = np.asarray(N)[:, None]
    lm = np.asarray(lam, dtype=float)[:, None]
    a_sc = np.broadcast_to(np.asarray(alpha_sc, dtype=float).reshape(-1, 1), N1.shape)
    N0 = 1.0
    cos0 = np.cos(np.radians(np.asarray(angles_deg, dtype=float)))[None, :]

    (rs, rp), (ts, tp), (rs_o, rp_o), (ts_o, tp_o), c1 = \
        _front_coeffs(N0, N1, lm, cos0, d_rough)

    # back face is bare substrate/air
    rs_b, rp_b, ts_b, tp_b = _fresnel(N1, N0, c1, cos0)

    # single-pass amplitude through the slab (absorption from Im N) times the
    # non-KK scatter loss, which is an intensity loss -> amplitude gets half of it
    phase = np.exp(1j * 2 * np.pi * N1 * c1 * d_nm / lm)
    haze = np.exp(-0.5 * a_sc * d_nm / np.real(c1))
    a1 = phase * haze

    M = _mueller_diag(rp, rs)                      # front-face beam
    for m in range(1, n_back + 1):
        # m-th emergent beam: in, (m-1) internal round trips, back bounce, out
        ap = tp[..., ] * tp_o * rp_b * (rp_b * rp_o) ** (m - 1) * a1 ** (2 * m)
        asx = ts * ts_o * rs_b * (rs_b * rs_o) ** (m - 1) * a1 ** (2 * m)
        M = M + _mueller_diag(ap, asx)
    return M / M[..., 0, 0][..., None, None]


def transmittance(N, lam, angles_deg, d_nm, alpha_sc=0.0, d_rough=0.0):
    """Total intensity transmittance of the slab (incoherent multiple beams).

    Closed form of the geometric series:  T = T_f T_b A / (1 - R_i R_b A^2),
    evaluated separately for s and p and averaged (unpolarised / normal incidence
    they coincide).  Returns [n_wl, n_ang].
    """
    N1 = np.asarray(N)[:, None]
    lm = np.asarray(lam, dtype=float)[:, None]
    a_sc = np.broadcast_to(np.asarray(alpha_sc, dtype=float).reshape(-1, 1), N1.shape)
    N0 = 1.0
    cos0 = np.cos(np.radians(np.asarray(angles_deg, dtype=float)))[None, :]

    (rs, rp), (ts, tp), (rs_o, rp_o), (ts_o, tp_o), c1 = \
        _front_coeffs(N0, N1, lm, cos0, d_rough)
    rs_b, rp_b, ts_b, tp_b = _fresnel(N1, N0, c1, cos0)

    # intensity transmittance of an interface carries the (N cos) flux factor
    flux_in = (N1 * c1).real / (N0 * cos0).real
    flux_out = (N0 * cos0).real / (N1 * c1).real

    A = np.exp(-4 * np.pi * N1.imag * d_nm / lm / np.real(c1)) * \
        np.exp(-a_sc * d_nm / np.real(c1))

    out = []
    for t_in, t_out, r_i, r_b in ((ts, ts_b, rs_o, rs_b), (tp, tp_b, rp_o, rp_b)):
        Tf = np.abs(t_in) ** 2 * flux_in
        Tb = np.abs(t_out) ** 2 * flux_out
        Ri = np.abs(r_i) ** 2
        Rb = np.abs(r_b) ** 2
        out.append(Tf * Tb * A / (1.0 - Ri * Rb * A ** 2))
    return 0.5 * (out[0] + out[1])


# --------------------------------------------------------------------------
# Self-test: energy conservation and limiting cases.
# --------------------------------------------------------------------------

if __name__ == "__main__":
    lam = np.array([500.0, 1000.0, 1500.0])
    n = np.array([1.842, 1.816, 1.807])
    N = n.astype(complex)
    d = 1.0e6  # 1 mm

    # lossless slab at normal incidence: T + R must equal 1
    T = transmittance(N, lam, [0.0], d)[:, 0]
    M = reflection_mm(N, lam, [0.0], d, n_back=60)
    # normalised MM loses the absolute R, so recompute R from the analytic series
    R0 = ((1 - n) / (1 + n)) ** 2
    R = R0 * (1 + (1 - R0) ** 2 / (1 - R0 ** 2))
    print("Lossless slab, normal incidence:")
    for i, l in enumerate(lam):
        print(f"  {l:6.0f} nm  n={n[i]:.4f}  T={T[i]:.6f}  R={R[i]:.6f}  T+R={T[i]+R[i]:.6f}")

    # isotropic reflection MM must satisfy mm12 == mm21 and mm33 == mm44
    ang = np.array([30.0, 55.0, 75.0])
    MM = reflection_mm(N, lam, ang, d)
    print("\nIsotropic symmetry of the reflection MM (should be ~0):")
    print(f"  max|mm12-mm21| = {np.max(np.abs(MM[:,:,0,1]-MM[:,:,1,0])):.2e}")
    print(f"  max|mm33-mm44| = {np.max(np.abs(MM[:,:,2,2]-MM[:,:,3,3])):.2e}")
    print(f"  max|mm34+mm43| = {np.max(np.abs(MM[:,:,2,3]+MM[:,:,3,2])):.2e}")

    # a strongly absorbing slab must approach the single-interface result
    Nabs = n + 1j * 0.05
    MMa = reflection_mm(Nabs, lam, ang, d)
    c0 = np.cos(np.radians(ang))[None, :]
    c1 = _cos_theta(1.0, Nabs[:, None], c0)
    rs, rp, _, _ = _fresnel(1.0, Nabs[:, None], c0, c1)
    rho = rp / rs
    psi = np.arctan(np.abs(rho))
    dl = np.angle(rho)
    print("\nOpaque limit vs single-interface Fresnel:")
    print(f"  max|mm12 + cos2Psi| = {np.max(np.abs(MMa[:,:,0,1] + np.cos(2*psi))):.2e}")
    print(f"  max|mm33 - sin2Psi cosD| = {np.max(np.abs(MMa[:,:,2,2] - np.sin(2*psi)*np.cos(dl))):.2e}")

    # zero-thickness roughness must reproduce the bare interface exactly
    m0 = reflection_mm(N, lam, ang, d, d_rough=0.0)
    m1 = reflection_mm(N, lam, ang, d, d_rough=1e-12)
    print(f"\nd_rough -> 0 continuity: max diff = {np.max(np.abs(m0-m1)):.2e}")
