#!/usr/bin/env python3

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "figures"
SUMMARY = ROOT / "data" / "optical_summary.csv"
HC_EV_NM = 1239.841984

SAMPLES = [
    "Instrument baseline",
    "5% Yb",
    "10% Yb",
    "15% Yb",
    "Pure sample B",
    "Pure sample A",
]

COLORS = {
    "Instrument baseline": "#888888",
    "5% Yb": "#0072B2",
    "10% Yb": "#D55E00",
    "15% Yb": "#009E73",
    "Pure sample B": "#CC79A7",
    "Pure sample A": "#7A5195",
}

LINESTYLES = {
    "Instrument baseline": ":",
    "5% Yb": "-",
    "10% Yb": "-",
    "15% Yb": "-",
    "Pure sample B": "--",
    "Pure sample A": "-.",
}

FEATURE_WINDOWS = [
    ("Feature 1", 810, 840),
    ("Feature 2", 845, 875),
    ("Feature 3", 880, 910),
]


def load_export(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.reader(handle))

    numeric_rows = []
    for row in rows[2:]:
        if len(row) < 12 or not all(value.strip() for value in row[:12]):
            continue
        try:
            numeric_rows.append([float(value) for value in row[:12]])
        except ValueError:
            continue

    values = np.asarray(numeric_rows)
    spectra = {}
    for index, sample in enumerate(SAMPLES):
        wavelength = values[:, 2 * index]
        signal = values[:, 2 * index + 1]
        order = np.argsort(wavelength)
        spectra[sample] = (wavelength[order], signal[order])
    return spectra


def align_spectra(first, second):
    first_wavelength, first_signal = first
    second_wavelength, second_signal = second
    wavelength, first_index, second_index = np.intersect1d(
        first_wavelength,
        second_wavelength,
        return_indices=True,
    )
    return wavelength, first_signal[first_index], second_signal[second_index]


def find_local_feature(wavelength, absorbance, low, high):
    smoothed = np.convolve(absorbance, np.ones(5) / 5, mode="same")
    region = np.flatnonzero((wavelength >= low) & (wavelength <= high))
    interior = region[1:-1]
    candidates = interior[
        (smoothed[interior] > smoothed[interior - 1])
        & (smoothed[interior] >= smoothed[interior + 1])
    ]
    if candidates.size:
        peak_index = candidates[np.argmax(smoothed[candidates])]
    else:
        peak_index = region[np.argmax(smoothed[region])]
    return wavelength[peak_index], absorbance[peak_index]


def add_energy_axis(axis, location="top"):
    def transform(wavelength):
        with np.errstate(divide="ignore", invalid="ignore"):
            return HC_EV_NM / np.asarray(wavelength)

    energy_axis = axis.secondary_xaxis(location, functions=(transform, transform))
    energy_axis.set_xlabel("Photon energy (eV)")


def plot_overview(transmission, absorbance):
    figure, axes = plt.subplots(2, 1, figsize=(12, 9), sharex=True)

    for sample in SAMPLES:
        wavelength, signal = transmission[sample]
        axes[0].plot(
            wavelength,
            signal,
            color=COLORS[sample],
            linestyle=LINESTYLES[sample],
            linewidth=1.8,
            label=sample,
        )

        wavelength, signal = absorbance[sample]
        axes[1].plot(
            wavelength,
            signal,
            color=COLORS[sample],
            linestyle=LINESTYLES[sample],
            linewidth=1.8,
            label=sample,
        )

    axes[0].set(
        title="Transmission",
        ylabel="Transmission (%)",
        xlim=(800, 1200),
    )
    axes[1].set(
        title="Absorbance",
        xlabel="Wavelength (nm)",
        ylabel="Absorbance",
        xlim=(800, 1200),
    )
    axes[0].legend(ncol=3, frameon=True)
    add_energy_axis(axes[0])

    figure.suptitle("Yb:YAG optical spectra")
    figure.tight_layout()
    figure.savefig(
        FIGURES / "YbYAG_transmission_absorbance.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(figure)


def analyze_doped_samples(transmission, absorbance):
    summary_rows = []
    comparison_data = {}

    for concentration, sample in ((5, "5% Yb"), (10, "10% Yb"), (15, "15% Yb")):
        wavelength, transmission_signal, direct_absorbance = align_spectra(
            transmission[sample],
            absorbance[sample],
        )
        calculated_absorbance = -np.log10(np.clip(transmission_signal, 1e-9, None) / 100)
        offset = np.median(direct_absorbance - calculated_absorbance)
        adjusted_absorbance = calculated_absorbance + offset
        residual = direct_absorbance - adjusted_absorbance

        features = [
            find_local_feature(wavelength, direct_absorbance, low, high)
            for _, low, high in FEATURE_WINDOWS
        ]

        row = {
            "yb_percent": concentration,
            "transmission_absorbance_correlation": np.corrcoef(
                direct_absorbance,
                calculated_absorbance,
            )[0, 1],
            "absorbance_offset": offset,
            "residual_rmse": np.sqrt(np.mean(residual**2)),
        }
        for (feature_name, _, _), (feature_wavelength, _) in zip(
            FEATURE_WINDOWS,
            features,
        ):
            key = feature_name.lower().replace(" ", "_")
            row[f"{key}_nm"] = feature_wavelength
            row[f"{key}_ev"] = HC_EV_NM / feature_wavelength

        summary_rows.append(row)
        comparison_data[sample] = {
            "concentration": concentration,
            "wavelength": wavelength,
            "direct": direct_absorbance,
            "adjusted": adjusted_absorbance,
            "residual": residual,
            "features": features,
        }

    return summary_rows, comparison_data


def write_summary(rows):
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    with SUMMARY.open("w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def plot_comparisons(comparison_data):
    figure = plt.figure(figsize=(13, 9), constrained_layout=True)
    grid = figure.add_gridspec(2, 2)
    feature_axis = figure.add_subplot(grid[0, 0])
    agreement_axis = figure.add_subplot(grid[0, 1])
    residual_axis = figure.add_subplot(grid[1, :])

    concentrations = [
        values["concentration"]
        for values in comparison_data.values()
    ]
    for feature_index, (feature_name, _, _) in enumerate(FEATURE_WINDOWS):
        positions = [
            values["features"][feature_index][0]
            for values in comparison_data.values()
        ]
        feature_axis.plot(
            concentrations,
            positions,
            marker="o",
            linewidth=1.8,
            label=feature_name,
        )

    feature_axis.set(
        title="Low-wavelength feature positions",
        xlabel="Yb concentration (%)",
        ylabel="Wavelength (nm)",
        xticks=concentrations,
    )
    feature_axis.legend(frameon=True)

    for sample, values in comparison_data.items():
        agreement_axis.scatter(
            values["adjusted"][::4],
            values["direct"][::4],
            s=15,
            alpha=0.7,
            color=COLORS[sample],
            label=sample,
        )
        residual_axis.plot(
            values["wavelength"],
            values["residual"],
            color=COLORS[sample],
            linewidth=1.5,
            label=sample,
        )

    limits = (0, 3)
    agreement_axis.plot(limits, limits, color="#444444", linestyle=":", linewidth=1)
    agreement_axis.set(
        title="Direct vs transmission-derived absorbance",
        xlabel="Derived absorbance, offset adjusted",
        ylabel="Direct absorbance",
        xlim=limits,
        ylim=limits,
    )
    agreement_axis.legend(frameon=True)

    residual_axis.axhline(0, color="#444444", linestyle=":", linewidth=1)
    residual_axis.set(
        title="Absorbance consistency residual",
        xlabel="Wavelength (nm)",
        ylabel="Direct - derived absorbance",
        xlim=(800, 1200),
    )
    residual_axis.legend(ncol=3, frameon=True)
    add_energy_axis(residual_axis)

    figure.suptitle("Yb:YAG optical comparisons")
    figure.savefig(
        FIGURES / "YbYAG_optical_comparisons.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(figure)


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")

    transmission = load_export(ROOT / "YbYag_T.csv")
    absorbance = load_export(ROOT / "YbYag_ABS.csv")

    plot_overview(transmission, absorbance)
    summary_rows, comparison_data = analyze_doped_samples(
        transmission,
        absorbance,
    )
    write_summary(summary_rows)
    plot_comparisons(comparison_data)


if __name__ == "__main__":
    main()
