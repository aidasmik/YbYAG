#!/usr/bin/env python3

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

try:
    from scripts.plot_optical_spectra import apply_academic_style, format_axis
except ModuleNotFoundError:
    from plot_optical_spectra import apply_academic_style, format_axis


ROOT = Path(__file__).resolve().parents[1]
PL_DIR = ROOT / "PL"
CLEAN_DIR = ROOT / "data" / "pl"
FIGURE_PATH = ROOT / "figures" / "YbYAG_photoluminescence.png"
TEMP_FIGURE_PATH = ROOT / "figures" / "YbYAG_temperature_photoluminescence.png"
REFLECTION_FIGURE_PATH = ROOT / "figures" / "YbYAG_reflection.png"
ABSORBANCE_FIGURE_PATH = ROOT / "figures" / "YbYAG_doped_absorbance.png"
ABSORBANCE_PATH = ROOT / "YbYag_ABS.csv"
REFLECTION_DIR = ROOT / "data" / "reflection"

SAMPLES = {
    "5% Yb": ("YbYag_5_PL.csv", 5, 200),
    "10% Yb": ("YbYag_10_PL.csv", 10, 1000),
    "15% Yb": ("YbYag_15_PL.csv", 15, 1000),
}

TEMPERATURE_SAMPLES = {
    "4 K": "4K.csv",
    "50 K": "50K.csv",
    "100 K": "100K.csv",
    "200 K": "200K.csv",
    "300 K": "300K.csv",
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

TEMP_COLORS = {
    "4 K": "#332288",
    "50 K": "#117733",
    "100 K": "#44AA99",
    "200 K": "#DDCC77",
    "300 K": "#CC6677",
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


def load_temperature_pl(path):
    rows = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.reader(handle):
            if len(row) < 4:
                continue
            try:
                wavelength, response, raw, intensity = map(parse_number, row[:4])
            except ValueError:
                continue
            if not np.all(np.isfinite((wavelength, response, raw, intensity))):
                continue

            # A few Origin CSV cells lost their decimal separator on export.
            while intensity > 100:
                intensity /= 1000

            rows.append((wavelength, response, raw, intensity))

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


def write_temperature_clean_csv(label, values):
    output = CLEAN_DIR / f"{label.replace(' ', '')}_PL_clean.csv"
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "wavelength_nm",
                "lamp_response",
                "raw_pl_counts",
                "intensity_arb_u",
            ]
        )
        writer.writerows(values)


def plot_pl_panels(pl_spectra):
    figure, axis = plt.subplots(figsize=(7.0, 4.5))

    for label, values in pl_spectra.items():
        wavelength = values[:, 0]
        intensity = values[:, 3]
        axis.plot(
            wavelength,
            intensity,
            color=COLORS[label],
            linewidth=1.55,
            label=label,
        )
    axis.set(
        title="Yb:YAG photoluminescence under 942 nm excitation",
        xlabel="Wavelength (nm)",
        ylabel="PL intensity (arb. u.)",
        xlim=(950, 1140),
    )
    format_axis(axis)
    axis.legend()

    figure.tight_layout()
    figure.savefig(FIGURE_PATH, bbox_inches="tight")
    plt.close(figure)


def plot_temperature_pl(temp_spectra):
    figure, axis = plt.subplots(figsize=(7.0, 4.5))

    for label, values in temp_spectra.items():
        wavelength = values[:, 0]
        intensity = values[:, 3]
        color = TEMP_COLORS[label]

        axis.plot(
            wavelength,
            intensity,
            color=color,
            linewidth=1.55,
            label=label,
        )

    axis.set(
        title="Yb:YAG temperature-dependent photoluminescence",
        xlabel="Wavelength (nm)",
        ylabel="Intensity (arb. u.)",
    )
    format_axis(axis)
    axis.legend(frameon=True)

    figure.tight_layout()
    figure.savefig(TEMP_FIGURE_PATH, bbox_inches="tight")
    plt.close(figure)


def plot_reflection_spectra(reflection, reflection_sides):
    figure, axis = plt.subplots(figsize=(7.0, 4.5))
    for label in SAMPLES:
        values = reflection[label]
        axis.plot(
            values[:, 0],
            values[:, 1],
            color=COLORS[label],
            linewidth=1.5,
            label=f"{label} (side {reflection_sides[label]})",
        )

    axis.set(
        title="Yb:YAG reflection spectra",
        xlabel="Wavelength (nm)",
        ylabel="Reflection (%)",
        xlim=(800, 1200),
    )
    format_axis(axis)
    axis.legend()
    figure.tight_layout()
    figure.savefig(REFLECTION_FIGURE_PATH, bbox_inches="tight")
    plt.close(figure)


def plot_absorbance_spectra(absorbance):
    figure, axis = plt.subplots(figsize=(7.0, 4.5))
    for label in SAMPLES:
        values = absorbance[label]
        axis.plot(
            values[:, 0],
            values[:, 1],
            color=COLORS[label],
            linewidth=1.5,
            label=label,
        )

    axis.set(
        title="Yb:YAG absorbance spectra",
        xlabel="Wavelength (nm)",
        ylabel="Absorbance",
        xlim=(800, 1200),
    )
    format_axis(axis)
    axis.legend()
    figure.tight_layout()
    figure.savefig(ABSORBANCE_FIGURE_PATH, bbox_inches="tight")
    plt.close(figure)


def main():
    CLEAN_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    apply_academic_style()

    pl_spectra = {}
    for label, (filename, _, correction_scale) in SAMPLES.items():
        values = load_pl(PL_DIR / filename, correction_scale)
        pl_spectra[label] = values
        write_clean_csv(label, values)

    temp_spectra = {}
    for label, filename in TEMPERATURE_SAMPLES.items():
        path = PL_DIR / filename
        if not path.exists():
            continue
        values = load_temperature_pl(path)
        temp_spectra[label] = values
        write_temperature_clean_csv(label, values)

    absorbance = load_absorbance(ABSORBANCE_PATH)
    reflection, reflection_sides = load_higher_reflection_sides()
    plot_pl_panels(pl_spectra)
    if temp_spectra:
        plot_temperature_pl(temp_spectra)
    plot_reflection_spectra(reflection, reflection_sides)
    plot_absorbance_spectra(absorbance)


if __name__ == "__main__":
    main()
