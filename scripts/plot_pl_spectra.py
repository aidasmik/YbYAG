#!/usr/bin/env python3

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PL_DIR = ROOT / "PL"
CLEAN_DIR = ROOT / "data" / "pl"
FIGURE_PATH = ROOT / "figures" / "YbYAG_photoluminescence.png"
ABSORBANCE_PATH = ROOT / "YbYag_ABS.csv"
REFLECTION_DIR = ROOT / "data" / "reflection"

SAMPLES = {
    "5% Yb": ("YbYag_5_PL.csv", 5, 200),
    "10% Yb": ("YbYag_10_PL.csv", 10, 1000),
    "15% Yb": ("YbYag_15_PL.csv", 15, 1000),
}

REFLECTION_FILES = {
    "5% Yb": {"A": "YbYag_5_A_R.txt", "B": "YbYag_5_B_R.txt"},
    "10% Yb": {"A": "YbYag_10_A_R.txt", "B": "YbYag_10_B_R.txt"},
    "15% Yb": {"A": "YbYag_15_A_R.txt", "B": "YbYag_15_B_R.txt"},
}

COLORS = {
    "5% Yb": "#0072B2",
    "10% Yb": "#D55E00",
    "15% Yb": "#009E73",
}

SIGNAL_COLORS = {
    "PL emission": "#0072B2",
    "Absorbance": "#D55E00",
    "Reflection": "#009E73",
}


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


def load_absorbance(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.reader(handle))

    spectra = {}
    for sample_index, label in enumerate(SAMPLES, start=1):
        values = []
        for row in rows[2:]:
            wavelength_index = 2 * sample_index
            signal_index = wavelength_index + 1
            if len(row) <= signal_index:
                continue
            try:
                values.append(
                    (float(row[wavelength_index]), float(row[signal_index]))
                )
            except ValueError:
                continue
        values = np.asarray(values)
        spectra[label] = values[np.argsort(values[:, 0])]
    return spectra


def load_reflection(path):
    values = []
    with path.open(encoding="utf-8-sig") as handle:
        for line in handle:
            fields = line.strip().replace('"', "").split(",")
            if len(fields) != 4:
                continue
            try:
                wavelength = float(f"{fields[0]}.{fields[1]}")
                reflection = float(f"{fields[2]}.{fields[3]}")
            except ValueError:
                continue
            values.append((wavelength, reflection))
    return np.asarray(values)


def load_higher_reflection_sides():
    spectra = {}
    selected_sides = {}
    for label, sides in REFLECTION_FILES.items():
        side_spectra = {
            side: load_reflection(REFLECTION_DIR / filename)
            for side, filename in sides.items()
        }
        selected_side = max(
            side_spectra,
            key=lambda side: np.mean(side_spectra[side][:, 1]),
        )
        spectra[label] = side_spectra[selected_side]
        selected_sides[label] = selected_side
    return spectra, selected_sides


def normalize_in_range(values, lower=950, upper=1125):
    wavelength = values[:, 0]
    signal = values[:, 1]
    selected = (wavelength >= lower) & (wavelength <= upper)
    wavelength = wavelength[selected]
    signal = signal[selected]
    signal = signal - np.min(signal)
    maximum = np.max(signal)
    if maximum > 0:
        signal = signal / maximum
    return wavelength, signal


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


def plot(pl_spectra, absorbance, reflection, reflection_sides):
    figure, axes = plt.subplots(2, 2, figsize=(12, 9))

    for label, values in pl_spectra.items():
        wavelength = values[:, 0]
        intensity = values[:, 3]
        axes[0, 0].plot(
            wavelength,
            intensity,
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

    overlay_axes = (axes[0, 1], axes[1, 0], axes[1, 1])
    for axis, label in zip(overlay_axes, SAMPLES):
        pl_values = np.column_stack(
            (pl_spectra[label][:, 0], pl_spectra[label][:, 3])
        )
        datasets = (
            ("PL emission", pl_values),
            ("Absorbance", absorbance[label]),
            ("Reflection", reflection[label]),
        )
        for signal_label, values in datasets:
            wavelength, signal = normalize_in_range(values)
            legend_label = signal_label
            if signal_label == "Reflection":
                legend_label += f" (side {reflection_sides[label]})"
            axis.plot(
                wavelength,
                signal,
                color=SIGNAL_COLORS[signal_label],
                linewidth=1.7,
                label=legend_label,
            )
        axis.set(
            title=label,
            xlabel="Wavelength (nm)",
            ylabel="Normalized signal",
            xlim=(950, 1125),
            ylim=(-0.03, 1.05),
        )
        axis.legend(fontsize=8)

    figure.suptitle(
        "Yb:YAG emission, absorbance, and reflection "
        "(PL excitation: 942 nm)"
    )
    figure.text(
        0.5,
        0.008,
        "Overlays are independently normalized for spectral-position comparison. "
        "Absolute PL intensities reflect the recorded measurement settings.",
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

    pl_spectra = {}
    for label, (filename, _, correction_scale) in SAMPLES.items():
        values = load_pl(PL_DIR / filename, correction_scale)
        pl_spectra[label] = values
        write_clean_csv(label, values)

    absorbance = load_absorbance(ABSORBANCE_PATH)
    reflection, reflection_sides = load_higher_reflection_sides()
    plot(pl_spectra, absorbance, reflection, reflection_sides)


if __name__ == "__main__":
    main()
