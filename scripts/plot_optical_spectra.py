#!/usr/bin/env python3

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "figures"
HC_EV_NM = 1239.841984

EXPORT_SAMPLES = [
    "Instrument baseline",
    "5% Yb",
    "10% Yb",
    "15% Yb",
    "Pure sample B",
    "Pure sample A",
]

PLOT_SAMPLES = EXPORT_SAMPLES[1:]

COLORS = {
    "5% Yb": "#0072B2",
    "10% Yb": "#D55E00",
    "15% Yb": "#009E73",
    "Pure sample B": "#CC79A7",
    "Pure sample A": "#7A5195",
}

LINESTYLES = {
    "5% Yb": "-",
    "10% Yb": "-",
    "15% Yb": "-",
    "Pure sample B": "--",
    "Pure sample A": "-.",
}


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
    for index, sample in enumerate(EXPORT_SAMPLES):
        wavelength = values[:, 2 * index]
        signal = values[:, 2 * index + 1]
        order = np.argsort(wavelength)
        spectra[sample] = (wavelength[order], signal[order])
    return spectra


def add_energy_axis(axis, location="top"):
    def transform(wavelength):
        with np.errstate(divide="ignore", invalid="ignore"):
            return HC_EV_NM / np.asarray(wavelength)

    energy_axis = axis.secondary_xaxis(location, functions=(transform, transform))
    energy_axis.set_xlabel("Photon energy (eV)")


def plot_overview(transmission, absorbance):
    figure, axes = plt.subplots(2, 1, figsize=(12, 9), sharex=True)

    for sample in PLOT_SAMPLES:
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


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")

    transmission = load_export(ROOT / "YbYag_T.csv")
    absorbance = load_export(ROOT / "YbYag_ABS.csv")

    plot_overview(transmission, absorbance)


if __name__ == "__main__":
    main()
