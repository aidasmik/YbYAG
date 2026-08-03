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

Here `n` and `k` come from **one model with 29 parameters** fitted to every
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
degeneracy the model exists to break. Each line is therefore held inside its own
narrow centre window (see the line table below) with widths <= 0.070 eV -- which
also stops neighbouring Yb lines merging into one broad smear.

## Which transmission dips are real

The measured transmittance contains **41 dips with prominence > 0.004**, and the
first version of this model described only the Yb cluster. Before adding
oscillators for the rest, they have to be shown to be *material* absorption --
fitting instrumental structure would inject fake `k` and, through KK, corrupt `n`.

Test: reflection was measured on **2026-06-30**, transmission on **2026-06-18**.
Real absorption changes `k`, which attenuates the backside beam and therefore
shows up in the reflection Mueller matrix and in the depolarisation. An artefact
confined to the transmission channel does not. Correlating the transmittance
ripple against the independently measured reflection observables:

| band | corr(T, depol) | corr(T, s12) | verdict |
|---|---|---|---|
| 300-500 nm | -0.09 | +0.04 | instrumental |
| 500-700 nm | -0.06 | +0.08 | instrumental |
| 700-900 nm | +0.05 | -0.22 | instrumental |
| **900-1060 nm** | **+0.34** | **-0.42** | **real absorption** |
| 1100-1690 nm | +0.06 | +0.20 | instrumental |

Only the Yb manifold is real. The ~30 other dips -- at 349, 403, 467, 486, 544,
581, 602, 611, 756, the dense 800-880 nm series, 1105, 1527 nm -- live purely in
the transmission channel (lamp/grating structure and a ~4 nm-period ripple in the
800-950 nm region). They are reproduced between samples A and B at r = 0.98-0.99
*because both share the same instrument*, so sample-to-sample agreement alone
does not establish they are real; the reflection cross-check does. **They are
deliberately not fitted.**

Checked and excluded as explanations for the 800-950 nm ripple: it is not a
detector-interleave artefact (consecutive-difference sign alternation is 41-44%,
i.e. random, and the even/odd offset is 6e-5), and an FFT in wavenumber returns
inconsistent optical paths per band (151 / 2.8 / 15.9 um), so it is not a single
etalon.

## The Yb manifold: centres are predicted, not fitted

Line centres are **fixed** at the Yb(3+) Stark-level scheme (`YB_STARK_GROUND`,
`YB_STARK_EXCITED` in `run_fit.py`; e.g. Fan et al., IEEE JQE **24**, 924 (1988)):

```
2F7/2 (ground):   0, 565, 612, 785 cm^-1
2F5/2 (excited):  10327, 10634, 10927 cm^-1
```

Every electronic line is a difference of one level from each set, so its position
is *predicted*. Only oscillator strengths and widths are fitted. This matters:
when the centres were free, a line drifted onto whatever residual bump was
nearby and parked itself at 955 nm with a 52 nm width against its bound -- a
"crystal-field line" wider than the entire manifold structure. Fixing the centres
also makes the problem tractable; with free centres neighbouring lines are
mutually degenerate and the fit needed 697 function evaluations instead of 45.

Transitions closer than 2 nm are merged (the band is not resolved below that).
Vibronic sidebands are declared separately and are the only oscillators allowed
to be broad; electronic lines are capped at 0.025 eV.

**The 955 nm feature is the ZPL + 150 cm^-1 phonon replica** -- 150 cm^-1 is the
lowest YAG optical phonon. It is now an explicitly labelled vibronic sideband at
954.5 nm rather than a floating oscillator, and the residual across 950-965 nm
drops below the photometric floor.

**The detector mask was also wrong.** It ran 965-983 nm and swallowed the real
zero-phonon line: at 965-970 nm the inverted optical depth rises smoothly to
`alpha*d = 0.48` at 969 nm, while only at 971-983 nm does the transmittance jump
to 0.60-0.68 -- above the sample's own baseline (~0.58) and therefore unphysical.
With the ZPL masked out the fit had no oscillator for the absorption climbing
into the mask, which is what the spurious 955 nm line was compensating. The mask
is now 971-983 nm and the ZPL is fitted.

## Results (sample B, d = 1.000 mm assumed)

`n` matches Zelmon single-crystal YAG to **RMS 0.0018** (max deviation 0.0047)
over 350-1650 nm. Scattering: `C = 3.505 +/- 0.002 /cm`, `p = 0.1102 +/- 0.0010`.
Sellmeier: UV pole 256.0 eV^2 at 10.640 eV, IR pole 0.0018 eV^2.

Yb(3+) oscillators (kind: e = electronic Stark transition, v = vibronic sideband;
all centres fixed):

| centre (nm) | kind | A | FWHM (eV) | assignment (cm^-1) |
|---|---|---|---|---|
| 875.1 | v | 2.10e-6 +/- 3e-7 | 0.100* | 10927 + 500 |
| 898.7 | v | 3.46e-6 +/- 5e-7 | 0.0368 | 10927 + 200 |
| 915.2 | e | 1.92e-5 +/- 9e-7 | 0.0140 | 0 -> 10927 |
| 930.5 | v | 2.11e-5 +/- 8e-7 | 0.0162 | ZPL + 420 |
| **940.4** | e | **3.72e-5 +/- 1e-6** | 0.0118 | 0 -> 10634 (main pump band) |
| **954.5** | v | **1.46e-5 +/- 4e-7** | 0.0831 | **ZPL + 150** |
| 965.1 | e | 1.6e-6 +/- 1e-6 | 0.0020* | 565 -> 10927 |
| **968.5** | e | **2.09e-5 +/- 1e-6** | 0.0040 | 0 -> 10327 (ZPL) + 612 -> 10927 |
| 993.1 | e | 1.27e-6 +/- 1e-6 | 0.0026* | 565 -> 10634 |
| 997.8 | e | 3.62e-6 +/- 1e-6 | 0.0062 | 612 -> 10634 |
| 1024.4 | e | 1.51e-6 +/- 2e-6 | 0.0191 | 565 -> 10327 |
| 1029.3 | e | 6.89e-6 +/- 2e-6 | 0.0183 | 612 -> 10327 |
| 1048.0 | e | 2.13e-6 +/- 9e-7 | 0.0151 | 785 -> 10327 |

\* width at a bound; that parameter's uncertainty is not meaningful.

**Not detected** and dropped by the automatic prune: 986.0 nm (785->10927) and
1015.3 nm (785->10634). Both originate on the 785 cm^-1 level, whose Boltzmann
population at 300 K is only `exp(-785/208) = 2%`, so weakness is expected. An
oscillator refined to zero amplitude has a width with zero gradient, which makes
the Jacobian singular and returns NaN for *every* uncertainty in the fit -- hence
the prune-and-refit pass (`AMP_DETECT`).

Transmittance residual inside the manifold (880-1080 nm):

| model | RMS | max |
|---|---|---|
| 3 broad Gaussians | 0.0087 | 0.0198 |
| 8 free-centre lines | 0.0026 | 0.0060 |
| **13 fixed Stark/vibronic lines** | **0.0036** | **0.0093** |

The free-centre model fits marginally better but is not physical -- its lines sit
where the residual happens to be, not where Yb(3+) has transitions. The
fixed-centre model is the one to quote.

**The near-grey scatter exponent `p = 0.110 +/- 0.001` remains the main physical
result.** The dominant loss is neither absorption nor Rayleigh scattering: it is
almost wavelength-independent, i.e. Mie/geometric scattering from grain-scale
features. That is what the earlier fit's unexplained 2-3% NIR transmittance gap
was.

**Cross-validation.** Sample A, fitted transmission-only with `n` frozen at B's
reflection result, reproduces B's oscillator strengths to 0.1-0.9%:

| line | B | A |
|---|---|---|
| 915.2 nm | 1.919e-5 | 1.920e-5 |
| 940.4 nm | 3.722e-5 | 3.693e-5 |
| 954.5 nm (ZPL+150) | 1.462e-5 | 1.449e-5 |
| 968.5 nm (ZPL) | 2.092e-5 | 2.090e-5 |
| scatter C, p | 3.505, 0.110 | 3.530, 0.131 |

## A bug worth recording

`_idx("gauss12_amp")` parsed the oscillator index as `int(head[-1])` -- the
trailing *character*. With ten or more oscillators that maps gauss10..gauss13
onto 0..3, silently aliasing four parameters onto four others. Symptoms: every
uncertainty NaN (singular Jacobian), amplitudes driven to absurd values, 697
function evaluations instead of 45 -- while the printed per-line table, which
reads `model["gaussians"]` directly, still looked plausible. Fixed, with a
regression test in `python fit/dispersion.py` that round-trips 14 oscillators.

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

Masked: 651-661 nm (grating turret change) and 971-983 nm (Si/InGaAs detector
crossover). The crossover sits just red of the Yb zero-phonon line; the inverted `alpha*d` inside
971-983 nm goes *negative*, confirming it is an artefact rather than structure.
`--keep-zpl` lifts the mask entirely.

## Known residuals

* `s34 = (mm34-mm43)/2` shows a systematic rise towards the blue, ~0.012 at
  260 nm decaying to ~0.002 in the NIR, which an isotropic model predicts as
  zero. This is the weak anisotropy/strain signal noted in the earlier analysis.
  Fitting it needs a uniaxial substrate; expected `delta_n` is ~1e-3.
* Surface roughness is implemented but left at 0 -- `n` already matches the YAG
  reference to 0.0017 without it, so there is nothing for it to absorb.
