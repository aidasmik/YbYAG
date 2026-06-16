#!/usr/bin/env python3

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from plot_optical_spectra import (
    COLORS,
    FIGURES,
    ROOT,
    apply_updated_pure_spectra,
    load_export,
)


N_YAG_1060NM = 1.81523
FIT_WINDOW_NM = (1100.0, 1200.0)
SUMMARY = ROOT / "data" / "YbYAG_pure_A_B_coating_fit_summary.csv"


def surface_reflectance(index):
    return ((index - 1.0) / (index + 1.0)) ** 2


def incoherent_slab_transmission(r_front, r_back):
    return (1.0 - r_front) * (1.0 - r_back) / (1.0 - r_front * r_back)


def symmetric_reflectance_from_transmission(transmission):
    return (1.0 - transmission) / (1.0 + transmission)


def index_from_reflectance(reflectance):
    root = np.sqrt(reflectance)
    return (1.0 + root) / (1.0 - root)


def one_side_reflectance_from_transmission(transmission, known_side_reflectance):
    numerator = 1.0 - known_side_reflectance - transmission
    denominator = 1.0 - known_side_reflectance - transmission * known_side_reflectance
    return numerator / denominator


def median_in_window(spectrum, window=FIT_WINDOW_NM):
    wavelength, signal = spectrum
    lo, hi = window
    mask = (wavelength >= lo) & (wavelength <= hi) & np.isfinite(signal)
    return float(np.median(signal[mask]))


def build_summary(transmission):
    yag_surface_r = surface_reflectance(N_YAG_1060NM)
    bare_yag_t = incoherent_slab_transmission(yag_surface_r, yag_surface_r)
    rows = []

    for sample in ("Pure sample A", "Pure sample B"):
        median_t = median_in_window(transmission[sample]) / 100.0
        symmetric_r = symmetric_reflectance_from_transmission(median_t)
        rows.append(
            {
                "sample": sample,
                "fit_window_nm": f"{FIT_WINDOW_NM[0]:.0f}-{FIT_WINDOW_NM[1]:.0f}",
                "median_transmission_fraction": median_t,
                "median_transmission_percent": 100.0 * median_t,
                "absorbance_from_transmission": -np.log10(median_t),
                "symmetric_surface_reflectance_fit": symmetric_r,
                "equivalent_symmetric_index": index_from_reflectance(symmetric_r),
                "one_side_reflectance_fit_assuming_other_side_bare_yag": one_side_reflectance_from_transmission(
                    median_t,
                    yag_surface_r,
                ),
                "bare_yag_surface_reflectance_model": yag_surface_r,
                "bare_yag_slab_transmission_model": bare_yag_t,
            }
        )

    return rows


def save_summary(rows):
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    with SUMMARY.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def plot_fit(transmission, rows):
    FIGURES.mkdir(parents=True, exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")

    yag_surface_r = surface_reflectance(N_YAG_1060NM)
    bare_yag_t = incoherent_slab_transmission(yag_surface_r, yag_surface_r)
    b_row = next(row for row in rows if row["sample"] == "Pure sample B")
    b_one_side_r = b_row["one_side_reflectance_fit_assuming_other_side_bare_yag"]
    b_fit_t = incoherent_slab_transmission(yag_surface_r, b_one_side_r)

    figure, axis = plt.subplots(figsize=(10.5, 6.0))
    for sample in ("Pure sample A", "Pure sample B"):
        wavelength, signal = transmission[sample]
        axis.plot(
            wavelength,
            signal,
            color=COLORS[sample],
            linewidth=2.0,
            label=sample,
        )

    axis.axhline(
        100.0 * bare_yag_t,
        color="#222222",
        linewidth=1.8,
        linestyle="--",
        label=f"Uncoated YAG slab model ({100.0 * bare_yag_t:.2f}%T)",
    )
    axis.axhline(
        100.0 * b_fit_t,
        color="#D55E00",
        linewidth=1.8,
        linestyle="--",
        label=f"B fit: one side R={100.0 * b_one_side_r:.1f}%",
    )

    axis.axvspan(
        FIT_WINDOW_NM[0],
        FIT_WINDOW_NM[1],
        color="#999999",
        alpha=0.18,
        label="fit window",
    )

    axis.set(
        title="Pure Yb:YAG A/B coating match fit",
        xlabel="Wavelength (nm)",
        ylabel="Transmission (%)",
        xlim=(800, 1200),
        ylim=(45, 90),
    )
    axis.legend(frameon=True, loc="lower right")
    figure.tight_layout()
    figure.savefig(
        FIGURES / "YbYAG_pure_A_B_coating_fit.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(figure)


def main():
    transmission = load_export(ROOT / "YbYag_T.csv")
    absorbance = load_export(ROOT / "YbYag_ABS.csv")
    apply_updated_pure_spectra(transmission, absorbance)

    rows = build_summary(transmission)
    save_summary(rows)
    plot_fit(transmission, rows)

    for row in rows:
        print(
            f"{row['sample']}: T={row['median_transmission_percent']:.3f}% "
            f"R_sym={100.0 * row['symmetric_surface_reflectance_fit']:.3f}% "
            f"n_equiv={row['equivalent_symmetric_index']:.4f} "
            f"R_one_side={100.0 * row['one_side_reflectance_fit_assuming_other_side_bare_yag']:.3f}%"
        )


if __name__ == "__main__":
    main()
