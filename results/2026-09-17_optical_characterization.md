# Yb:YAG optical-characterization checkpoint — 2026-09-17

This note records the conclusions that are currently supported by the measurements discussed during the September 2026 analysis. It deliberately separates measured results, model-dependent inference, and unresolved questions.

## Main conclusion

The uncoated `B` crystal should **not currently be treated as a confirmed 5 at.% Yb:YAG sample**.

The independent optical evidence is more consistent with an effective Yb concentration of roughly

\[
C_{\mathrm{Yb}} \approx 1.5\text{–}2.5\ \mathrm{at.\%},
\]

with the absorption measurements themselves centering closer to **1.6–1.7 at.%**. This is an optical inference, not a direct chemical composition measurement.

A direct composition measurement (EPMA/WDS, ICP-OES/ICP-MS, or calibrated XRF/EDS) is still needed before assigning the concentration definitively.

## 1. Which RC2 transmission measurement is trustworthy?

Three RC2 transmission runs were compared.

| Run | Transparent-region baseline (1060–1100 nm) | Interpretation |
|---|---:|---|
| 2026-06-16 17:51 | ~84.05% | physically consistent with an uncoated YAG plate |
| 2026-06-18 13:14 | ~60.35% | large throughput/alignment loss |
| 2026-06-18 13:22 | ~59.25% | large throughput/alignment loss |

For `n ≈ 1.83`, an uncoated lossless YAG parallel plate is expected to transmit about 84% after Fresnel losses and incoherent internal reflections. The June-16 run therefore has the correct absolute baseline.

The June-18 spectra can be brought almost on top of the June-16 spectrum by multiplying them by approximately 1.405 and 1.430, respectively, outside the anomalous ~970–985 nm region. This strongly suggests a stable optical-throughput/alignment problem rather than a change in the crystal.

![RC2 transmission run comparison](../figures/analysis_2026_09_17/rc2_transmission_runs_comparison.svg)

## 2. Current CompleteEASE B-spline dielectric function

The current B-spline solution is physically consistent with the good June-16 transmission over the main Yb absorption bands.

Important features of the current imaginary dielectric function:

- ~915 nm band: \(\varepsilon_2 \sim 2.5\times10^{-5}\)
- ~941 nm band: \(\varepsilon_2 \sim 4.8\times10^{-5}\)
- ~969 nm ZPL: \(\varepsilon_2 \sim 3.7\times10^{-5}\)
- weak ~1005 nm structure: \(\varepsilon_2 \sim 1.1\times10^{-5}\)
- ~1030 nm absorption is still effectively forced to zero by the B-spline fit

![Current epsilon2](../figures/analysis_2026_09_17/current_b_spline_epsilon2.svg)

The fitted dielectric function reproduces the measured 1-mm transmission very well around 915, 941, and 969 nm. The previous comparison gave an overall 850–1100 nm transmission RMSE of approximately 1.2 percentage points.

![Fit versus good transmission](../figures/analysis_2026_09_17/fit_vs_good_transmission.svg)

### Remaining weakness: 1030 nm

The good transmission measurement indicates weak but non-zero absorption near 1030 nm:

\[
\alpha(1030\ \mathrm{nm}) \approx 0.25\ \mathrm{cm^{-1}}.
\]

The current B-spline sets \(\varepsilon_2\) there to approximately zero, so this weak reabsorption feature is missing from the model.

For \(n\approx1.83\), the missing absorption corresponds to only \(\varepsilon_2\) of order \(7\times10^{-6}\), so it is plausible that the B-spline regularization is suppressing a real but weak feature.

## 3. Absorption-derived Yb concentration

The good June-16 transmission gives approximately:

| Wavelength | Measured \(\alpha\) | Literature \(\sigma_a\) used | Inferred Yb |
|---|---:|---:|---:|
| 941 nm | ~1.88 cm⁻¹ | \(7.89\times10^{-21}\) cm² | ~1.71 at.% |
| 969 nm | ~1.36 cm⁻¹ | \(7.60\times10^{-21}\) cm² | ~1.28 at.% |
| 1030 nm | ~0.25 cm⁻¹ | \(1.18\times10^{-21}\) cm² | ~1.55 at.% |

The 969-nm value is the least trustworthy concentration estimate because the zero-phonon line is narrow and therefore much more sensitive to spectral resolution, wavelength registration, and the nearby measurement/emission artifact.

The agreement between **941 nm (~1.7 at.%)** and **1030 nm (~1.6 at.%)** is especially important because they are spectrally very different transitions.

![Cross section comparison](../figures/analysis_2026_09_17/cross_section_vs_literature.svg)

## 4. PL evidence

The room-temperature PL ratio decreases systematically across the nominal reference samples:

| Sample | \(I_{969}/I_{1030}\) |
|---|---:|
| nominal 5 at.% | ~0.2565 |
| nominal 10 at.% | ~0.2068 |
| nominal 15 at.% | ~0.1851 |
| uncoated / unknown, 300 K | ~0.2857 |

A crude extrapolation of the 5→10 at.% trend gives:

- point ratio \(I_{969}/I_{1030}\): ~2.1 at.%
- integrated 966–972 nm / 1020–1040 nm ratio: ~2.4 at.%

PL therefore supports a concentration below nominal 5 at.%, but it should be treated as secondary evidence because fluorescence line shape depends on temperature, self-absorption, collection geometry, excitation geometry, and coatings.

![PL ratio](../figures/analysis_2026_09_17/pl_ratio_vs_nominal_concentration.svg)

![PL spectral comparison](../figures/analysis_2026_09_17/pl_300K_vs_nominal_concentrations.svg)

## 5. Coated nominal-5% samples

The absolute coated-sample reflectance files are not absolutely calibrated: values exceed 100% in parts of the HR band. They should therefore not be converted directly into an absolute absorption coefficient.

However, after normalizing coated samples A and B outside the absorption bands, the additional optical loss of sample B,

\[
-\ln(R_B/R_A),
\]

has a spectral shape that correlates strongly with the uncoated sample's Yb absorption over 925–975 nm:

\[
r \approx 0.95.
\]

This indicates that the relative depressions near 941 and 969 nm are genuinely Yb-like rather than arbitrary coating structure.

![Coated versus uncoated spectral shape](../figures/analysis_2026_09_17/coated_vs_uncoated_yb_shape.svg)

Under a simplified double-pass model, the 941-nm relative loss would correspond to about 5 at.% if the coated active layer is ~245–250 µm thick. This is plausible for a thin-disk medium, but the result is not an absolute concentration measurement because the coating transfer function is not modeled.

## 6. Combined concentration evidence

![Independent concentration estimates](../figures/analysis_2026_09_17/inferred_yb_concentration_evidence.svg)

The current working interpretation is:

- absorption at 941 nm: ~1.7 at.%
- weak 1030-nm ground-state reabsorption: ~1.6 at.%
- 969-nm absorption: ~1.3 at.% but resolution-sensitive
- PL ratios: ~2.1–2.4 at.% by crude extrapolation

Therefore:

\[
\boxed{C_{\mathrm{Yb}}\approx1.5\text{–}2.5\ \mathrm{at.\%}}
\]

is a defensible working range, with ~1.6–2.0 at.% especially plausible from the absorption data.

## 7. Literature values used

The repository stores the numerical reference values in `data/analysis_2026_09_17/literature_reference_values.csv`.

### Brown et al. (2005) — single-crystal Yb:YAG

D. C. Brown, R. L. Cone, Y. Sun, and R. W. Equall, *Yb:YAG Absorption at Ambient and Cryogenic Temperatures*, IEEE JSTQE 11, 604–612 (2005). DOI: [10.1109/JSTQE.2005.850236](https://doi.org/10.1109/JSTQE.2005.850236)

Key 300-K values:

- doping equivalence: 1 at.% Yb = \(1.385\times10^{20}\) ions cm⁻³
- 941-nm pump-band peak: 941.025 nm
- 941-nm FWHM: ~18.3 nm
- ZPL peak: 968.825 nm
- ZPL FWHM: ~2.8 nm

Brown et al. also showed that normalized absorption cross sections were essentially concentration independent across their measured crystal concentrations.

### Liu et al. (2007) — Yb:YAG crystal rod

Q. Liu, X. Fu, M. Gong, and L. Huang, *Effects of the temperature dependence of absorption coefficients in edge-pumped Yb:YAG slab lasers*, JOSA B 24, 2081–2089 (2007). DOI: [10.1364/JOSAB.24.002081](https://doi.org/10.1364/JOSAB.24.002081)

Their experiment used a 4.5 at.% Yb:YAG crystal rod. At 23 °C they report:

- peak wavelength: 941.2 nm
- peak absorption cross section: \(7.89\times10^{-21}\) cm²
- FWHM: ~19.8 nm

### Körner et al. (2012) — temperature-dependent Yb:YAG spectra

J. Körner et al., *Measurement of temperature-dependent absorption and emission spectra of Yb:YAG, Yb:LuAG, and Yb:CaF2 between 20 °C and 200 °C and predictions on their influence on laser performance*, JOSA B 29, 2493–2502 (2012). DOI: [10.1364/JOSAB.29.002493](https://doi.org/10.1364/JOSAB.29.002493)

Analytic fits to these spectra, reproduced in M. Zeyen's ETH thesis, give at 300 K approximately:

- \(\sigma_a(941\,\mathrm{nm}) \approx 8.15\times10^{-21}\) cm²
- \(\sigma_a(969\,\mathrm{nm}) \approx 7.60\times10^{-21}\) cm²
- \(\sigma_a(1030\,\mathrm{nm}) \approx 1.18\times10^{-21}\) cm²

The weak 1030-nm ground-state absorption is particularly useful here because the measured \(\alpha(1030)\) independently implies ~1.6 at.% Yb.

### De Vido et al. (2020) — high-resolution ZPL reference

M. De Vido, A. Wojtusiak, and K. Ertel, *High-resolution absorption measurement at the zero phonon line of Yb:YAG between 80 K and 300 K*, Optical Materials Express 10, 717–723 (2020). DOI: [10.1364/OME.386436](https://doi.org/10.1364/OME.386436)

At 300 K:

- ZPL center: 969.04 nm
- peak absorption cross section: \(0.8\times10^{-20}\) cm²
- FWHM: 2.38 nm

Note: this experiment used a 1.1 at.% Yb:YAG ceramic; it is used here mainly as a high-resolution ZPL reference.

### Pirri et al. (2018) — Yb:YAG comparison spectrum

A. Pirri et al., *A Comprehensive Characterization of a 10 at.% Yb:YSAG Laser Ceramic Sample*, Materials 11, 837 (2018). DOI: [10.3390/ma11050837](https://doi.org/10.3390/ma11050837)

The paper includes a Yb:YAG ceramic comparison spectrum with:

- main YAG absorption peak near 939.4 nm
- ZPL near 968.93 nm

The repository's `figures/YbYAG_literature_band_comparison.png` is a redrawn comparison using these published spectral positions.

## 8. What should be done next?

1. Treat the uncoated sample concentration as **unknown**, not fixed to 5 at.%.
2. Retain the current B-spline as the empirical dielectric-function reference.
3. Add enough freedom around ~1030 nm to recover the weak measured reabsorption without distorting the 941/969 bands.
4. Re-measure absolute transmission with careful beam centering/reference checks.
5. Obtain an independent direct composition measurement.
6. Once concentration is known, recompute the absolute absorption cross section from the measured dielectric function.

## Data provenance

- PL and coated-reflection source files were already present in this repository.
- The analysis files added with this checkpoint contain processed/derived data from the September 2026 discussion.
- Literature figures are **not copied** into the repository; comparison plots are redrawn from cited numerical values and spectral positions.
