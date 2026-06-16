#!/usr/bin/env python3

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from plot_optical_spectra import FIGURES, ROOT, load_export
from plot_pl_spectra import load_reflection


SUMMARY = ROOT / "data" / "YbYAG_coated_sample_coating_match_summary.csv"
CANDIDATE_SUMMARY = ROOT / "data" / "YbYAG_coated_sample_material_candidates.csv"
FIT_FIGURE = FIGURES / "YbYAG_coated_sample_coating_match.png"

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

MATERIAL_PAIRS = {
    "Ta2O5/SiO2": (2.10, 1.45),
    "Nb2O5/SiO2": (2.25, 1.45),
    "HfO2/SiO2": (1.90, 1.45),
    "TiO2/SiO2": (2.35, 1.45),
    "Al2O3/SiO2": (1.76, 1.45),
}


def median_in_window(wavelength, signal, lo, hi):
    mask = (wavelength >= lo) & (wavelength <= hi) & np.isfinite(signal)
    return float(np.median(signal[mask]))


def moving_average(values, width=11):
    kernel = np.ones(width, dtype=float) / width
    return np.convolve(values, kernel, mode="same")


def lower_half_max_edge(wavelength, normalized_reflectance):
    mask = (wavelength >= 800.0) & (wavelength <= 1000.0)
    x = wavelength[mask]
    y = moving_average(normalized_reflectance[mask])
    crossing = np.flatnonzero(y >= 0.5)
    if crossing.size == 0:
        return np.nan
    index = int(crossing[0])
    if index == 0:
        return float(x[index])
    x0, x1 = x[index - 1], x[index]
    y0, y1 = y[index - 1], y[index]
    if y1 == y0:
        return float(x1)
    return float(x0 + (0.5 - y0) * (x1 - x0) / (y1 - y0))


def dbr_half_band_fraction(high_index, low_index):
    contrast = (high_index - low_index) / (high_index + low_index)
    return 2.0 * np.arcsin(contrast) / np.pi


def summarize_reflection():
    rows = []
    for sample, sides in REFLECTION_FILES.items():
        for side, filename in sides.items():
            values = load_reflection(ROOT / "data" / "reflection" / filename)
            wavelength = values[:, 0]
            reflectance = values[:, 1]
            plateau = median_in_window(wavelength, reflectance, 1000.0, 1150.0)
            normalized = reflectance / plateau
            rows.append(
                {
                    "sample": sample,
                    "side": side,
                    "plateau_reflectance_percent_1000_1150": plateau,
                    "lower_50_percent_edge_nm": lower_half_max_edge(
                        wavelength,
                        normalized,
                    ),
                    "normalized_reflectance_800_850": median_in_window(
                        wavelength,
                        normalized,
                        800.0,
                        850.0,
                    ),
                    "normalized_reflectance_900_930": median_in_window(
                        wavelength,
                        normalized,
                        900.0,
                        930.0,
                    ),
                    "normalized_reflectance_1020_1060": median_in_window(
                        wavelength,
                        normalized,
                        1020.0,
                        1060.0,
                    ),
                    "normalized_reflectance_1100_1200": median_in_window(
                        wavelength,
                        normalized,
                        1100.0,
                        1200.0,
                    ),
                }
            )
    return rows


def summarize_transmission():
    spectra = load_export(ROOT / "YbYag_T.csv")
    rows = {}
    for sample in REFLECTION_FILES:
        wavelength, transmission = spectra[sample]
        rows[sample] = {
            "transmission_percent_1030": median_in_window(
                wavelength,
                transmission,
                1028.0,
                1032.0,
            ),
            "transmission_percent_1064": median_in_window(
                wavelength,
                transmission,
                1062.0,
                1066.0,
            ),
            "transmission_percent_1100_1200": median_in_window(
                wavelength,
                transmission,
                1100.0,
                1200.0,
            ),
        }
    return rows


def candidate_matches(lower_edge_nm):
    matches = []
    for material_pair, (high_index, low_index) in MATERIAL_PAIRS.items():
        half_fraction = dbr_half_band_fraction(high_index, low_index)
        center = lower_edge_nm / (1.0 - half_fraction)
        upper = center * (1.0 + half_fraction)
        score = abs(center - 1030.0) / 50.0
        if upper < 1200.0:
            score += (1200.0 - upper) / 80.0
        if lower_edge_nm < 860.0 or lower_edge_nm > 920.0:
            score += 1.0
        matches.append(
            {
                "candidate_material_pair": material_pair,
                "dbr_center_from_measured_lower_edge_nm": center,
                "dbr_upper_edge_estimate_nm": upper,
                "dbr_fractional_half_width": half_fraction,
                "relative_score_lower_is_better": score,
            }
        )
    return sorted(matches, key=lambda row: row["relative_score_lower_is_better"])


def write_summary(reflection_rows, transmission_rows, match_rows):
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(reflection_rows[0]) + [
        "transmission_percent_1030",
        "transmission_percent_1064",
        "transmission_percent_1100_1200",
        "best_candidate_material_pair",
        "best_candidate_center_nm",
        "best_candidate_upper_edge_nm",
    ]
    best = match_rows[0]
    with SUMMARY.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in reflection_rows:
            output = dict(row)
            output.update(transmission_rows[row["sample"]])
            output.update(
                {
                    "best_candidate_material_pair": best["candidate_material_pair"],
                    "best_candidate_center_nm": best[
                        "dbr_center_from_measured_lower_edge_nm"
                    ],
                    "best_candidate_upper_edge_nm": best["dbr_upper_edge_estimate_nm"],
                }
            )
            writer.writerow(output)

    with CANDIDATE_SUMMARY.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(match_rows[0]),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(match_rows)


def plot_match(reflection_rows, match_rows):
    plt.style.use("seaborn-v0_8-whitegrid")
    FIGURES.mkdir(parents=True, exist_ok=True)

    figure, axis = plt.subplots(figsize=(11.0, 6.3))
    for sample, sides in REFLECTION_FILES.items():
        for side, filename in sides.items():
            values = load_reflection(ROOT / "data" / "reflection" / filename)
            wavelength = values[:, 0]
            reflectance = values[:, 1]
            plateau = median_in_window(wavelength, reflectance, 1000.0, 1150.0)
            normalized = reflectance / plateau
            axis.plot(
                wavelength,
                normalized,
                color=COLORS[sample],
                linestyle="-" if side == "A" else "--",
                linewidth=1.6,
                alpha=0.85,
                label=f"{sample} side {side}",
            )

    lower_edges = np.asarray(
        [row["lower_50_percent_edge_nm"] for row in reflection_rows],
        dtype=float,
    )
    measured_lower = float(np.nanmedian(lower_edges))
    best = match_rows[0]
    axis.axvline(
        measured_lower,
        color="#222222",
        linewidth=1.6,
        linestyle=":",
        label=f"median 50% lower edge {measured_lower:.0f} nm",
    )
    axis.axvspan(
        measured_lower,
        best["dbr_upper_edge_estimate_nm"],
        color="#999999",
        alpha=0.18,
        label=(
            f"{best['candidate_material_pair']} DBR band estimate, "
            f"center {best['dbr_center_from_measured_lower_edge_nm']:.0f} nm"
        ),
    )
    for wavelength, label in ((1030.0, "1030 nm"), (1064.0, "1064 nm")):
        axis.axvline(wavelength, color="#555555", linewidth=1.0, alpha=0.55)
        axis.text(
            wavelength + 4.0,
            0.08,
            label,
            rotation=90,
            va="bottom",
            ha="left",
            fontsize=9,
            color="#333333",
        )

    axis.set(
        title="Coated Yb:YAG samples: dielectric HR coating match",
        xlabel="Wavelength (nm)",
        ylabel="Reflectance normalized to 1000-1150 nm plateau",
        xlim=(800, 1200),
        ylim=(0, 1.25),
    )
    axis.legend(ncol=2, frameon=True, fontsize=8)
    figure.tight_layout()
    figure.savefig(FIT_FIGURE, dpi=200, bbox_inches="tight")
    plt.close(figure)


def main():
    reflection_rows = summarize_reflection()
    lower_edge_nm = float(
        np.nanmedian([row["lower_50_percent_edge_nm"] for row in reflection_rows])
    )
    match_rows = candidate_matches(lower_edge_nm)
    transmission_rows = summarize_transmission()

    write_summary(reflection_rows, transmission_rows, match_rows)
    plot_match(reflection_rows, match_rows)

    print(f"Median lower 50% HR edge: {lower_edge_nm:.1f} nm")
    for row in match_rows:
        print(
            f"{row['candidate_material_pair']}: "
            f"center={row['dbr_center_from_measured_lower_edge_nm']:.1f} nm, "
            f"upper={row['dbr_upper_edge_estimate_nm']:.1f} nm, "
            f"score={row['relative_score_lower_is_better']:.3f}"
        )


if __name__ == "__main__":
    main()
