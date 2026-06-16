#!/usr/bin/env python3

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "figures"
HC_EV_NM = 1239.841984
UPDATED_PURE_SPECTRA = ROOT / "data" / "YbYAG_A3_B1_last_spectra.csv"

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


def load_updated_pure_spectra(path):
    values = {
        "Pure sample A": {"wavelength": [], "absorbance": [], "transmission": []},
        "Pure sample B": {"wavelength": [], "absorbance": [], "transmission": []},
    }
    column_map = {
        "Pure sample A": ("A3 Abs", "A3 Transmittance (%)"),
        "Pure sample B": ("B1 Abs", "B1 Transmittance (%)"),
    }

    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            wavelength = float(row["Wavelength (nm)"])
            for sample, (absorbance_column, transmission_column) in column_map.items():
                values[sample]["wavelength"].append(wavelength)
                values[sample]["absorbance"].append(float(row[absorbance_column]))
                values[sample]["transmission"].append(float(row[transmission_column]))

    spectra = {"transmission": {}, "absorbance": {}}
    for sample, sample_values in values.items():
        wavelength = np.asarray(sample_values["wavelength"], dtype=float)
        order = np.argsort(wavelength)
        for key in ("transmission", "absorbance"):
            signal = np.asarray(sample_values[key], dtype=float)
            spectra[key][sample] = (wavelength[order], signal[order])
    return spectra


def apply_updated_pure_spectra(transmission, absorbance, path=UPDATED_PURE_SPECTRA):
    pure_spectra = load_updated_pure_spectra(path)
    for sample in ("Pure sample A", "Pure sample B"):
        transmission[sample] = pure_spectra["transmission"][sample]
        absorbance[sample] = pure_spectra["absorbance"][sample]


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


def plot_transmission_only(transmission):
    wavelength_figure, wavelength_axis = plt.subplots(figsize=(9, 5.4))
    energy_figure, energy_axis = plt.subplots(figsize=(9, 5.4))

    for sample in PLOT_SAMPLES:
        wavelength, signal = transmission[sample]
        wavelength_axis.plot(
            wavelength,
            signal,
            color=COLORS[sample],
            linestyle=LINESTYLES[sample],
            marker="o",
            markersize=2,
            linewidth=1.3,
            label=sample.replace("% Yb", "%").replace("Pure sample ", ""),
        )

        energy = HC_EV_NM / wavelength
        order = np.argsort(energy)
        energy_axis.plot(
            energy[order],
            signal[order],
            color=COLORS[sample],
            linestyle=LINESTYLES[sample],
            marker="o",
            markersize=2,
            linewidth=1.3,
            label=sample.replace("% Yb", "%").replace("Pure sample ", ""),
        )

    for axis, xlabel in (
        (wavelength_axis, "Wavelength (nm)"),
        (energy_axis, "Energy (eV)"),
    ):
        axis.set(
            title="Yb:YAG transmission spectrum",
            xlabel=xlabel,
            ylabel="Transmission (%T)",
        )
        axis.legend(title="Sample")

    wavelength_axis.set_xlim(780, 1220)
    energy_axis.set_xlim(1.0, 1.56)

    wavelength_figure.tight_layout()
    energy_figure.tight_layout()
    wavelength_figure.savefig(
        ROOT / "YbYAG_transmission_plot_points_small.png",
        dpi=160,
        bbox_inches="tight",
    )
    energy_figure.savefig(
        ROOT / "YbYAG_transmission_plot_eV.png",
        dpi=160,
        bbox_inches="tight",
    )
    plt.close(wavelength_figure)
    plt.close(energy_figure)


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")

    transmission = load_export(ROOT / "YbYag_T.csv")
    absorbance = load_export(ROOT / "YbYag_ABS.csv")
    apply_updated_pure_spectra(transmission, absorbance)

    plot_overview(transmission, absorbance)
    plot_transmission_only(transmission)


if __name__ == "__main__":
    main()
