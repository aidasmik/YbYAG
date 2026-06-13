#!/usr/bin/env python3

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PL_DIR = ROOT / "PL"
CLEAN_DIR = ROOT / "data" / "pl"
FIGURE_PATH = ROOT / "figures" / "YbYAG_photoluminescence.png"

SAMPLES = {
    "5% Yb": ("YbYag_5_PL.csv", 5, 200),
    "10% Yb": ("YbYag_10_PL.csv", 10, 1000),
    "15% Yb": ("YbYag_15_PL.csv", 15, 1000),
}

COLORS = {
    "5% Yb": "#0072B2",
    "10% Yb": "#D55E00",
    "15% Yb": "#009E73",
}

LITERATURE_PEAK_NM = 1030.0
LITERATURE_FWHM_NM = 8.5


def parse_number(value):
    value = value.strip()
    if not value or value == "--":
        return np.nan
    return float(value)


def load_pl(path, correction_scale):
    rows = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.reader(handle):
            if len(row) < 3:
                continue
            wavelength, response, raw = map(parse_number, row[:3])
            if not np.all(np.isfinite((wavelength, response, raw))):
                continue

            # A few Origin CSV cells lost their decimal separator on export.
            if response > 1000:
                response /= 1000

            corrected = raw / response * correction_scale
            rows.append((wavelength, response, raw, corrected))

    values = np.asarray(rows)
    order = np.argsort(values[:, 0])
    return values[order]


def peak_in_range(wavelength, intensity, lower, upper):
    selected = (wavelength >= lower) & (wavelength <= upper)
    local_wavelength = wavelength[selected]
    local_intensity = intensity[selected]
    index = np.argmax(local_intensity)
    return local_wavelength[index], local_intensity[index]


def main_peak_fwhm(wavelength, intensity):
    peak_wavelength, peak_intensity = peak_in_range(
        wavelength, intensity, 1010, 1045
    )
    half_maximum = peak_intensity / 2
    selected = (wavelength >= 1010) & (wavelength <= 1045)
    local_wavelength = wavelength[selected]
    local_intensity = intensity[selected]
    above_half = local_wavelength[local_intensity >= half_maximum]
    return peak_wavelength, above_half[-1] - above_half[0]


def write_clean_csv(label, values):
    output = CLEAN_DIR / f"{label.split('%')[0]}pct_Yb_PL_clean.csv"
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "wavelength_nm",
                "lamp_response",
                "raw_pl_counts",
                "corrected_pl_arb_u",
            ]
        )
        writer.writerows(values)


def analyze(values):
    wavelength = values[:, 0]
    intensity = values[:, 3]
    zero_line_nm, zero_line_intensity = peak_in_range(
        wavelength, intensity, 960, 980
    )
    main_peak_nm, main_peak_intensity = peak_in_range(
        wavelength, intensity, 1010, 1045
    )
    _, fwhm_nm = main_peak_fwhm(wavelength, intensity)
    return {
        "zero_line_nm": zero_line_nm,
        "main_peak_nm": main_peak_nm,
        "fwhm_nm": fwhm_nm,
        "zero_to_main_ratio": zero_line_intensity / main_peak_intensity,
    }


def plot(spectra, metrics):
    figure, axes = plt.subplots(2, 2, figsize=(12, 9))

    for label, values in spectra.items():
        wavelength = values[:, 0]
        intensity = values[:, 3]
        axes[0, 0].plot(
            wavelength,
            intensity,
            color=COLORS[label],
            linewidth=1.8,
            label=label,
        )
        axes[0, 1].plot(
            wavelength,
            intensity / np.max(intensity),
            color=COLORS[label],
            linewidth=1.8,
            label=label,
        )

    axes[0, 0].set(
        title="Response-corrected photoluminescence",
        xlabel="Wavelength (nm)",
        ylabel="PL intensity (arb. u.)",
        xlim=(950, 1125),
    )
    axes[0, 0].legend()

    axes[0, 1].axvline(
        968.93,
        color="#333333",
        linestyle=":",
        linewidth=1.3,
        label="968.93 nm zero line",
    )
    axes[0, 1].axvspan(
        LITERATURE_PEAK_NM - LITERATURE_FWHM_NM / 2,
        LITERATURE_PEAK_NM + LITERATURE_FWHM_NM / 2,
        color="#E69F00",
        alpha=0.18,
        label="Published 1030 nm band",
    )
    axes[0, 1].set(
        title="Normalized spectral shape",
        xlabel="Wavelength (nm)",
        ylabel="Normalized PL",
        xlim=(950, 1080),
        ylim=(0, 1.05),
    )
    axes[0, 1].legend(fontsize=8)

    concentrations = np.array(
        [SAMPLES[label][1] for label in SAMPLES], dtype=float
    )
    main_peaks = np.array([metrics[label]["main_peak_nm"] for label in SAMPLES])
    fwhm = np.array([metrics[label]["fwhm_nm"] for label in SAMPLES])
    ratios = np.array(
        [metrics[label]["zero_to_main_ratio"] for label in SAMPLES]
    )

    axes[1, 0].errorbar(
        concentrations,
        main_peaks,
        yerr=fwhm / 2,
        fmt="o",
        color="#0072B2",
        capsize=5,
        label="Measured peak +/- FWHM/2",
    )
    axes[1, 0].axhspan(
        LITERATURE_PEAK_NM - LITERATURE_FWHM_NM / 2,
        LITERATURE_PEAK_NM + LITERATURE_FWHM_NM / 2,
        color="#E69F00",
        alpha=0.18,
        label="Published 1030 nm, FWHM 8.5 nm",
    )
    axes[1, 0].set(
        title="Main emission band",
        xlabel="Nominal Yb concentration (%)",
        ylabel="Peak wavelength (nm)",
        xticks=concentrations,
        ylim=(1022, 1037),
    )
    axes[1, 0].legend(fontsize=8)

    axes[1, 1].plot(
        concentrations,
        ratios,
        marker="o",
        linewidth=1.8,
        color="#7A5195",
    )
    for concentration, ratio in zip(concentrations, ratios):
        axes[1, 1].text(
            concentration,
            ratio + 0.004,
            f"{ratio:.3f}",
            ha="center",
            fontsize=9,
        )
    axes[1, 1].set(
        title="Zero-line emission relative to 1030 nm",
        xlabel="Nominal Yb concentration (%)",
        ylabel=r"$I_{969}/I_{1030}$",
        xticks=concentrations,
        ylim=(0.15, 0.28),
    )

    figure.suptitle("Yb:YAG photoluminescence under 942 nm excitation")
    figure.text(
        0.5,
        0.008,
        "Published band: Pirri et al., Materials 11 (2018) 837. "
        "Absolute intensities reflect the recorded measurement settings.",
        ha="center",
        fontsize=8,
    )
    figure.tight_layout(rect=(0, 0.035, 1, 0.97))
    figure.savefig(FIGURE_PATH, dpi=200, bbox_inches="tight")
    plt.close(figure)


def main():
    CLEAN_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")

    spectra = {}
    metrics = {}
    for label, (filename, _, correction_scale) in SAMPLES.items():
        values = load_pl(PL_DIR / filename, correction_scale)
        spectra[label] = values
        metrics[label] = analyze(values)
        write_clean_csv(label, values)

    plot(spectra, metrics)


if __name__ == "__main__":
    main()
