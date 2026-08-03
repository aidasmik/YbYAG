# Physical fit of the Yb:YAG RC2 Mueller-matrix data

A global, Kramers-Kronig-consistent oscillator fit to the CompleteEASE exports in
`YbYag txt/`, replacing the earlier wavelength-by-wavelength inversion.

```
python -m fit.run_fit --sample B --thickness 1.0   # full R(10 AOI) + T fit
python -m fit.run_fit --sample A                   # reuses B's dispersion
python -m fit.run_fit --list                       # incl. screened-out files
python -m fit.plot_fit --sample B
python fit/dispersion.py                           # KK self-test
python fit/slab.py                                 # forward-model self-tests
```

## Why "physical"

The previous fit (`Desktop/Remote/YbYag/analysis/fit_slab.py`) let `n` and `k`
float independently at each of ~1000 wavelengths. It reproduced the data but
could not answer the question it raised: transmittance came out 2-3% high in the
NIR, meaning the loss needed to explain `T` exceeded what `k` could carry without
spoiling the reflection fit. With `n,k` free per point there is nothing to stop
absorption and scattering trading against each other.

Here `n` and `k` come from **one model with 14 parameters** fitted to every
wavelength and every angle simultaneously, so the trade is closed:

* the 10-AOI reflection Mueller matrix fixes `n(lambda)` absolutely, independent
  of the slab thickness;
* the normal-incidence transmittance fixes the total per-pass loss
  `(alpha_abs + alpha_sc) * d`;
* KK consistency forbids the absorptive part from taking an arbitrary shape --
  any `k` rising towards the NIR must show up in `n` -- so the smooth,
  non-KK remainder is pinned as scattering rather than absorption.

## Model

`fit/dispersion.py`. Everything inside `eps` is KK-consistent by construction:

```
eps(E) = eps_inf + sum_p A_p/(Ec_p^2 - E^2) + sum_g Gaussian(A_g, En_g, Br_g)
```

* **Sellmeier poles** -- one UV, one IR, seeded from Zelmon single-crystal YAG
  (Appl. Opt. **37**, 4933 (1998)) converted to energy space.
* **Gaussian oscillators** for the Yb(3+) `2F7/2 -> 2F5/2` lines, with the
  analytic KK partner via the Dawson function (the WVASE/CompleteEASE "Gaussian").
  Gaussian rather than Lorentzian because the lines are inhomogeneously
  broadened; a Lorentzian tail would leak absorption across the whole
  transparent window. `python fit/dispersion.py` verifies the analytic `eps1`
  against a numerical KK transform of `eps2` (agreement 1e-10 off resonance,
  2e-5 on resonance, limited by grid discretisation).

**Scattering is deliberately kept out of `eps`.** A ceramic/defect haze
redirects light rather than absorbing it, so it is not KK-consistent and forcing
it into `k` distorts both `n` and the oscillator strengths. It is carried as a
separate intensity loss with a **free exponent**:

```
alpha_sc(lambda) = C * (1000 nm / lambda)^p
```

The exponent is fitted, not fixed at the Rayleigh value. The measured loss falls
only ~20% from 300 to 1650 nm, whereas `p = 4` varies by ~900x across that span,
so a fixed `lambda^-4` term is simply driven to zero and leaves the NIR
transmittance ~2% under-predicted -- exactly the discrepancy the earlier fit
reported. The fitted `p` is a diagnostic of scatterer size.

## Forward model

`fit/slab.py`. A thick, double-side-polished slab in air, fringe-free at the
RC2's resolution, so front and back beams add **incoherently** -- Mueller
matrices are summed, not Jones amplitudes. This is what produces the measured
depolarisation and the angle-fanning of the pseudo-index. An optional
nm-thick surface-roughness layer (Bruggeman 50/50 void) is treated coherently.

`python fit/slab.py` checks, to machine precision:

| test | result |
|---|---|
| lossless slab, `T + R = 1` | exact |
| isotropy: `mm12=mm21`, `mm33=mm44`, `mm34=-mm43` | 0 |
| opaque limit reduces to single-interface Fresnel | 4e-16 |
| `d_rough -> 0` continuity | 5e-15 |

## Two data-handling points that materially change the answer

**1. The reported error bars are not usable as weights.** The RC2 reports
`sigma ~ 1e-4` on the normalised MM elements and `~6e-4` on transmittance. The
real systematic is ~70x larger: `mm12 - mm21` is consistently `-0.0075` although
isotropy requires the two to be equal, and **~10% of the measured `mm12` values
fall below -1**, which a normalised Mueller matrix cannot physically do.

The fit therefore uses the combinations isotropy predicts,

```
s12 = (mm12+mm21)/2   s33 = (mm33+mm44)/2   s34 = (mm34-mm43)/2
```

and takes the **error floor from the antisymmetric combinations**, which the
model predicts to be exactly zero and which consequently measure the instrument
systematic directly. Transmittance gets a 0.5% photometric floor (`--t-sys`).
Left at face value the `T` block outweighs the entire 10-AOI reflection data set
by ~100x and `n` is sacrificed to chase a systematic: with reported errors the
fit returns `n(1550 nm) = 1.73`; with honest floors it returns 1.8081 against the
YAG reference 1.8070.

**2. Oscillators must be confined to the Yb manifold.** With generic bounds one
Gaussian runs off to ~1570 nm with a 1.3 eV width and becomes a smooth background
absorption -- a KK-consistent stand-in for scattering, reintroducing the very
degeneracy the model exists to break. Centres are held to 900-1100 nm and widths
to <= 0.08 eV.

## Results (sample B, d = 1.000 mm assumed)

`n` matches Zelmon single-crystal YAG to **RMS 0.0017** (max deviation 0.0046)
over 350-1650 nm, from a 14-parameter global model.

| parameter | value | 1 sigma |
|---|---|---|
| UV pole amplitude | 258.2 eV^2 | 3.0 |
| UV pole energy | 10.682 eV | 0.058 |
| IR pole amplitude | 0.00433 eV^2 | 0.0026 |
| Gaussian 1 | A=2.603e-5, E=1.3305 eV (931.8 nm), Br=0.0696 eV | 4e-7 / 6e-4 / 1e-3 |
| Gaussian 2 | A=2.831e-5, E=1.3184 eV (940.4 nm), Br=0.0111 eV | 8e-7 / 1e-4 / 4e-4 |
| Gaussian 3 | A=5.53e-6, E=1.2278 eV (1009.8 nm), Br=0.080 eV | 4e-7 / 3e-3 / 9e-3 |
| scatter | C = 3.510 /cm, p = 0.109 | 0.002 / 0.001 |

Per-pass optical depth: scattering `alpha*d = 0.388` (400 nm) falling to `0.335`
(1550 nm); Yb absorption `alpha*d = 0.193` at the 940 nm band, against an
empirical `0.196` read directly off the inverted transmittance.

**The near-grey scatter exponent `p = 0.11` is the main physical result.** The
dominant loss in these samples is not absorption and not Rayleigh scattering: it
is almost wavelength-independent, i.e. Mie/geometric scattering from features
comparable to or larger than the wavelength (grain boundaries, residual pores).
That is what the earlier fit's unexplained 2-3% NIR transmittance gap was.

**Cross-validation.** Sample A, fitted transmission-only with `n` frozen at B's
reflection result, reproduces B's parameters independently:

| | scatter C | p | 931 nm A | 940 nm A | 1009 nm A |
|---|---|---|---|---|---|
| B (R+T) | 3.510 | 0.109 | 2.603e-5 | 2.831e-5 | 5.53e-6 |
| A (T only) | 3.535 | 0.131 | 2.600e-5 | 2.817e-5 | 5.13e-6 |

Agreement is <1% on every oscillator amplitude and 0.1 nm on every line centre.

## Thickness

Unchanged and unavoidable: a fringe-free thick slab measures only the product
`alpha*d`. `d` is an input (`--thickness`, default 1 mm). `n` is independent of
it; `k`, the oscillator amplitudes and the scatter coefficients all scale as
`1/d`. The `alphad_*` columns in the output CSV are the thickness-independent
quantities -- **quote those**, and multiply through once `d` is measured.

## Data screening

Verified from the file headers, not the filenames. CompleteEASE writes
`TransmissionMeas=1` for transmission and omits the key for reflection; the
exports also carry a UTF-8 BOM and CRLF endings that put the header on the
*second* line, so a first-line parse mislabels every transmission file.

* `YbYag_5a_MMt` / `YbYag_5a_UMMt` are named `...t` but are **reflection** scans
  at AOI 30-75 deg.
* **A and B are both Yb-doped.** Both show the 940 nm band at `alpha*d ~ 0.54`
  against a ~0.36 baseline. "Pure" in this repository means *uncoated*, not
  undoped.
* **5% cannot be fitted at all.** Transmission files are corrupt (inf, negative
  and >1 values; the entire 2026-06-23 session returned all-NaN), and *both* 5a
  reflection scans read `mm12 ~ 0`, `mm33 ~ -1` (Psi=45 deg, Delta=180 deg, no
  polarisation change) at AOI 30-50 deg -- a failed alignment. Needs re-measuring.
* **10% cannot be fitted.** Both transmission files contain negative and >1
  transmittance and no reflection scan exists.
* `YbYag_B_UMMt [2026-06-16,175112]` has median `T = 0.835`, sitting at the
  lossless-slab limit -- it looks like a straight-through baseline, not a sample.

Masked: 651-661 nm (grating turret change) and 965-980 nm (Si/InGaAs detector
crossover). The crossover unfortunately sits on the Yb zero-phonon line, and the
inverted `alpha*d` there goes *negative*, confirming it is an artefact rather
than structure. `--keep-zpl` lifts the mask and adds a 969 nm oscillator.

## Known residuals

* `s34 = (mm34-mm43)/2` shows a systematic rise towards the blue, ~0.012 at
  260 nm decaying to ~0.002 in the NIR, which an isotropic model predicts as
  zero. This is the weak anisotropy/strain signal noted in the earlier analysis.
  Fitting it needs a uniaxial substrate; expected `delta_n` is ~1e-3.
* Surface roughness is implemented but left at 0 -- `n` already matches the YAG
  reference to 0.0017 without it, so there is nothing for it to absorb.
