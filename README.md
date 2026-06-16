# Yb:YAG

[Reflection notebook](YbYAG_reflection.ipynb)

![Yb:YAG reflection measurements](figures/YbYAG_reflection.png)

[Optical plotting script](scripts/plot_optical_spectra.py)

![Transmission and absorbance](figures/YbYAG_transmission_absorbance.png)

![Pure sample A and B transmission and absorbance comparison](figures/YbYAG_pure_A_B_abs_trans_comparison.png)

![Pure sample B offset to match sample A](figures/YbYAG_pure_A_B_offset_alignment.png)

## Pure A/B coating fit

[Coating fit script](scripts/fit_pure_ab_coating.py)

![Pure sample A and B coating match fit](figures/YbYAG_pure_A_B_coating_fit.png)

Fit summary: [data/YbYAG_pure_A_B_coating_fit_summary.csv](data/YbYAG_pure_A_B_coating_fit_summary.csv).

Using the 1100-1200 nm baseline region, sample A matches an uncoated polished YAG slab. Sample B does not match a normal uncoated or AR-coated slab; the simplest optical fit is one bare YAG side plus one side with about 22.3% effective reflection/loss.

## Literature comparison

[Comparison script](scripts/compare_with_literature.py)

![Measured Yb:YAG absorbance with published spectral positions](figures/YbYAG_literature_band_comparison.png)

The published 939.4 and 968.93 nm absorption peaks are unresolved in the broad high-absorbance plateau.

Spectral positions from [Pirri et al., Materials 11 (2018) 837](https://doi.org/10.3390/ma11050837).

## Photoluminescence

[PL plotting script](scripts/plot_pl_spectra.py)

![Yb:YAG photoluminescence comparison](figures/YbYAG_photoluminescence.png)

PL spectra for all Yb concentrations are shown on the same plot.

![Combined measured spectra](figures/YbYAG_combined_spectra.png)

Full measured PL, reflection, and absorbance spectra on an aligned wavelength axis.
