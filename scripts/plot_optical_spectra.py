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


def apply_academic_style():
    plt.style.use("default")
    plt.rcParams.update(
        {
            "figure.dpi": 120,
            "savefig.dpi": 300,
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "legend.fontsize": 9,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "axes.linewidth": 0.9,
            "axes.grid": True,
            "grid.color": "#D9D9D9",
            "grid.linewidth": 0.55,
            "grid.alpha": 0.8,
            "legend.frameon": False,
        }
    )


def format_axis(axis):
    axis.tick_params(direction="in", which="both", top=True, right=True)
    axis.minorticks_on()


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
    energy_axis.tick_params(direction="in", which="both")


def plot_spectra(
    spectra,
    samples,
    title,
    ylabel,
    output_path,
    xlim=(800, 1200),
    energy_axis=False,
):
    figure, axis = plt.subplots(figsize=(7.0, 4.5))

    for sample in samples:
        wavelength, signal = spectra[sample]
        axis.plot(
            wavelength,
            signal,
            color=COLORS[sample],
            linestyle=LINESTYLES[sample],
            linewidth=1.55,
            label=sample,
        )

    axis.set(
        title=title,
        xlabel="Wavelength (nm)",
        ylabel=ylabel,
        xlim=xlim,
    )
    if energy_axis:
        add_energy_axis(axis)
    format_axis(axis)
    axis.legend(ncol=2)

    figure.tight_layout()
    figure.savefig(output_path, bbox_inches="tight")
    plt.close(figure)


def plot_overview(transmission, absorbance):
    plot_spectra(
        transmission,
        PLOT_SAMPLES,
        "Yb:YAG transmission spectra",
        "Transmission (%)",
        FIGURES / "YbYAG_transmission.png",
        energy_axis=True,
    )
    plot_spectra(
        absorbance,
        PLOT_SAMPLES,
        "Yb:YAG absorbance spectra",
        "Absorbance",
        FIGURES / "YbYAG_absorbance.png",
    )


def plot_transmission_only(transmission):
    wavelength_figure, wavelength_axis = plt.subplots(figsize=(7.0, 4.5))
    energy_figure, energy_axis = plt.subplots(figsize=(7.0, 4.5))

    for sample in PLOT_SAMPLES:
        wavelength, signal = transmission[sample]
        wavelength_axis.plot(
            wavelength,
            signal,
            color=COLORS[sample],
            linestyle=LINESTYLES[sample],
            linewidth=1.45,
            label=sample.replace("% Yb", "%").replace("Pure sample ", ""),
        )

        energy = HC_EV_NM / wavelength
        order = np.argsort(energy)
        energy_axis.plot(
            energy[order],
            signal[order],
            color=COLORS[sample],
            linestyle=LINESTYLES[sample],
            linewidth=1.45,
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
        format_axis(axis)

    wavelength_axis.set_xlim(780, 1220)
    energy_axis.set_xlim(1.0, 1.56)

    wavelength_figure.tight_layout()
    energy_figure.tight_layout()
    wavelength_figure.savefig(
        FIGURES / "YbYAG_transmission_wavelength.png",
        bbox_inches="tight",
    )
    energy_figure.savefig(
        FIGURES / "YbYAG_transmission_energy.png",
        bbox_inches="tight",
    )
    plt.close(wavelength_figure)
    plt.close(energy_figure)


def plot_pure_ab_comparison(transmission, absorbance):
    samples = ("Pure sample A", "Pure sample B")
    plot_spectra(
        transmission,
        samples,
        "Pure Yb:YAG A/B transmission",
        "Transmission (%)",
        FIGURES / "YbYAG_pure_A_B_transmission.png",
        energy_axis=True,
    )
    plot_spectra(
        absorbance,
        samples,
        "Pure Yb:YAG A/B absorbance",
        "Absorbance",
        FIGURES / "YbYAG_pure_A_B_absorbance.png",
    )


def offset_to_match(reference, target):
    reference_wavelength, reference_signal = reference
    target_wavelength, target_signal = target
    interpolated_target = np.interp(
        reference_wavelength,
        target_wavelength,
        target_signal,
        left=np.nan,
        right=np.nan,
    )
    finite = np.isfinite(reference_signal) & np.isfinite(interpolated_target)
    offset = float(np.nanmedian(reference_signal[finite] - interpolated_target[finite]))
    return reference_wavelength, interpolated_target + offset, offset


def plot_pure_ab_offset_alignment(transmission, absorbance):
    panels = (
        (
            transmission["Pure sample A"],
            transmission["Pure sample B"],
            "Pure Yb:YAG B transmission offset to match A",
            "Transmission (%)",
            "percentage points",
            FIGURES / "YbYAG_pure_A_B_transmission_offset_alignment.png",
        ),
        (
            absorbance["Pure sample A"],
            absorbance["Pure sample B"],
            "Pure Yb:YAG B absorbance offset to match A",
            "Absorbance",
            "absorbance",
            FIGURES / "YbYAG_pure_A_B_absorbance_offset_alignment.png",
        ),
    )

    for reference, target, title, ylabel, unit, output_path in panels:
        figure, axis = plt.subplots(figsize=(7.0, 4.5))
        reference_wavelength, reference_signal = reference
        target_wavelength, target_signal = target
        shifted_wavelength, shifted_signal, offset = offset_to_match(reference, target)

        axis.plot(
            reference_wavelength,
            reference_signal,
            color=COLORS["Pure sample A"],
            linestyle=LINESTYLES["Pure sample A"],
            linewidth=2.0,
            label="Pure sample A",
        )
        axis.plot(
            target_wavelength,
            target_signal,
            color=COLORS["Pure sample B"],
            linestyle=":",
            alpha=0.55,
            linewidth=2.0,
            label="Pure sample B original",
        )
        axis.plot(
            shifted_wavelength,
            shifted_signal,
            color=COLORS["Pure sample B"],
            linestyle=LINESTYLES["Pure sample B"],
            linewidth=1.65,
            label="Pure sample B + offset",
        )
        axis.set(
            title=f"{title} ({offset:+.4g} {unit})",
            xlabel="Wavelength (nm)",
            ylabel=ylabel,
            xlim=(800, 1200),
        )
        format_axis(axis)
        axis.legend(loc="best")

        figure.tight_layout()
        figure.savefig(output_path, bbox_inches="tight")
        plt.close(figure)


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)
    apply_academic_style()

    transmission = load_export(ROOT / "YbYag_T.csv")
    absorbance = load_export(ROOT / "YbYag_ABS.csv")
    apply_updated_pure_spectra(transmission, absorbance)

    plot_overview(transmission, absorbance)
    plot_transmission_only(transmission)
    plot_pure_ab_comparison(transmission, absorbance)
    plot_pure_ab_offset_alignment(transmission, absorbance)


if __name__ == "__main__":
    main()
