#!/usr/bin/env python3

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

try:
    from scripts.plot_optical_spectra import COLORS, FIGURES, ROOT, load_export
except ModuleNotFoundError:
    from plot_optical_spectra import COLORS, FIGURES, ROOT, load_export


REFLECTION_FILES = {
    "Pure sample A": ROOT / "data" / "reflection" / "YbYag_A_A.txt",
    "Pure sample B": ROOT / "data" / "reflection" / "YbYag_B_A.txt",
}

PURE_COLORS = {
    "Pure sample A": COLORS["Pure sample A"],
    "Pure sample B": COLORS["Pure sample B"],
}


def load_reflection(path: Path):
    wavelength = []
    reflection = []

    for line in path.read_text(encoding="ascii").splitlines()[2:]:
        fields = line.strip().split(",")
        if len(fields) != 4:
            continue
        wavelength.append(float(f"{fields[0]}.{fields[1]}"))
        reflection.append(float(f"{fields[2]}.{fields[3]}"))

    return np.asarray(wavelength), np.asarray(reflection)


def yag_refractive_index(wavelength_nm):
    wavelength_um = np.asarray(wavelength_nm) / 1000
    wavelength_squared = wavelength_um**2
    n_squared = (
        1
        + 2.28200 * wavelength_squared / (wavelength_squared - 0.01185)
        + 3.27644 * wavelength_squared / (wavelength_squared - 282.734)
    )
    return np.sqrt(n_squared)


def fresnel_spectra(wavelength_nm):
    refractive_index = yag_refractive_index(wavelength_nm)
    surface_reflection = ((refractive_index - 1) / (refractive_index + 1)) ** 2

    # Incoherent multiple reflections in a lossless, plane-parallel slab.
    slab_reflection = 2 * surface_reflection / (1 + surface_reflection)
    slab_transmission = (1 - surface_reflection) / (1 + surface_reflection)
    return 100 * surface_reflection, 100 * slab_reflection, 100 * slab_transmission


def plot_ybyag_band_comparison(absorbance):
    figure, axes = plt.subplots(
        2,
        1,
        figsize=(11, 7),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )

    for sample in ("5% Yb", "10% Yb", "15% Yb"):
        wavelength, signal = absorbance[sample]
        axes[0].plot(
            wavelength,
            signal,
            color=COLORS[sample],
            linewidth=1.8,
            label=sample,
        )

    absorption_peaks = (939.4, 968.93)
    for wavelength in absorption_peaks:
        axes[0].axvline(wavelength, color="#333333", linestyle=":", linewidth=1.3)

    axes[0].text(939.4, 2.92, "939.4", ha="right", va="top", fontsize=9)
    axes[0].text(968.93, 2.92, "968.93 nm", ha="left", va="top", fontsize=9)
    axes[0].set(
        title="Measured Yb:YAG absorbance and published spectral positions",
        ylabel="Measured absorbance",
        ylim=(0, 3.05),
    )
    axes[0].legend(ncol=3)

    axes[1].set_yticks([1, 0], ["Absorption", "Emission"])
    axes[1].scatter(
        absorption_peaks,
        [1, 1],
        marker="|",
        s=350,
        linewidths=2.5,
        color="#333333",
    )
    axes[1].scatter(
        [1030],
        [1],
        marker="o",
        s=45,
        facecolors="none",
        edgecolors="#333333",
    )
    axes[1].barh(
        0,
        width=8.5,
        left=1030 - 8.5 / 2,
        height=0.34,
        color="#E69F00",
        alpha=0.55,
    )
    axes[1].scatter(
        [1030, 1050],
        [0, 0],
        marker="D",
        s=[55, 35],
        color="#E69F00",
    )
    axes[1].set(
        xlabel="Wavelength (nm)",
        xlim=(880, 1100),
        ylim=(-0.55, 1.55),
    )
    axes[1].text(1034, 1.15, "weak near 1030 nm", fontsize=8)
    axes[1].text(1030, -0.34, "1030 nm, FWHM 8.5 nm", ha="center", fontsize=8)
    axes[1].text(1050, 0.18, "secondary", ha="center", fontsize=8)

    figure.text(
        0.5,
        0.005,
        "Published positions: Pirri et al., Materials 11 (2018) 837. "
        "Emission markers are spectral positions, not measured absorbance.",
        ha="center",
        fontsize=8,
    )
    figure.tight_layout(rect=(0, 0.035, 1, 1))
    figure.savefig(
        FIGURES / "YbYAG_literature_band_comparison.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(figure)


def plot_yag_fresnel_comparison(transmission):
    wavelength = np.linspace(800, 1200, 401)
    surface_r, slab_r, slab_t = fresnel_spectra(wavelength)
    reflection = {
        sample: load_reflection(path) for sample, path in REFLECTION_FILES.items()
    }

    figure, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True)

    for sample in ("Pure sample A", "Pure sample B"):
        measured_wavelength, measured_signal = reflection[sample]
        axes[0].plot(
            measured_wavelength,
            measured_signal,
            color=PURE_COLORS[sample],
            linewidth=1.7,
            label=sample,
        )

    axes[0].plot(
        wavelength,
        surface_r,
        color="#222222",
        linestyle=":",
        linewidth=1.8,
        label="Theory: one surface",
    )
    axes[0].plot(
        wavelength,
        slab_r,
        color="#222222",
        linestyle="--",
        linewidth=1.8,
        label="Theory: two-surface slab",
    )
    axes[0].set(
        title="Pure YAG reflection compared with Fresnel theory",
        ylabel="Reflection (%)",
        ylim=(0, 25),
    )
    axes[0].legend(ncol=2)

    for sample in ("Pure sample A", "Pure sample B"):
        measured_wavelength, measured_signal = transmission[sample]
        axes[1].plot(
            measured_wavelength,
            measured_signal,
            color=PURE_COLORS[sample],
            linewidth=1.7,
            label=sample,
        )

    axes[1].plot(
        wavelength,
        slab_t,
        color="#222222",
        linestyle="--",
        linewidth=1.8,
        label="Theory: lossless slab",
    )
    axes[1].set(
        title="Pure YAG transmission compared with Fresnel theory",
        xlabel="Wavelength (nm)",
        ylabel="Transmission (%)",
        xlim=(800, 1200),
        ylim=(50, 90),
    )
    axes[1].legend(ncol=3)

    figure.text(
        0.5,
        0.005,
        "Theory uses the YAG refractive-index model of Zelmon et al. "
        "and assumes an uncoated, nonabsorbing plane-parallel slab.",
        ha="center",
        fontsize=8,
    )
    figure.tight_layout(rect=(0, 0.035, 1, 1))
    figure.savefig(
        FIGURES / "YAG_fresnel_comparison.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(figure)


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")

    absorbance = load_export(ROOT / "YbYag_ABS.csv")
    transmission = load_export(ROOT / "YbYag_T.csv")

    plot_ybyag_band_comparison(absorbance)
    plot_yag_fresnel_comparison(transmission)


if __name__ == "__main__":
    main()
