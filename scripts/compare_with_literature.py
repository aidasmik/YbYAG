#!/usr/bin/env python3

import matplotlib.pyplot as plt

try:
    from scripts.plot_optical_spectra import COLORS, FIGURES, ROOT, load_export
except ModuleNotFoundError:
    from plot_optical_spectra import COLORS, FIGURES, ROOT, load_export


def plot_ybyag_band_comparison(absorbance):
    figure, axis = plt.subplots(figsize=(11, 6))

    for sample in ("5% Yb", "10% Yb", "15% Yb"):
        wavelength, signal = absorbance[sample]
        axis.plot(
            wavelength,
            signal,
            color=COLORS[sample],
            linewidth=1.8,
            label=sample,
        )

    for wavelength in (939.4, 968.93):
        axis.axvline(
            wavelength,
            color="#333333",
            linestyle=":",
            linewidth=1.5,
        )

    axis.axvline(
        1030,
        color="#E69F00",
        linestyle="--",
        linewidth=1.5,
    )
    axis.axvspan(
        1030 - 8.5 / 2,
        1030 + 8.5 / 2,
        color="#E69F00",
        alpha=0.18,
    )
    axis.axvline(
        1050,
        color="#E69F00",
        linestyle=":",
        linewidth=1.2,
    )

    axis.text(939.4, 2.92, "Abs. 939.4 nm", ha="right", va="top", fontsize=9)
    axis.text(968.93, 2.92, "Abs. 968.93 nm", ha="left", va="top", fontsize=9)
    axis.text(
        1030,
        0.07,
        "Emission 1030 nm\nFWHM 8.5 nm",
        ha="center",
        va="bottom",
        fontsize=9,
        color="#9C6B00",
    )
    axis.text(
        1050,
        0.07,
        "Emission\n~1050 nm",
        ha="left",
        va="bottom",
        fontsize=9,
        color="#9C6B00",
    )

    axis.set(
        title="Measured Yb:YAG absorbance and published spectral positions",
        xlabel="Wavelength (nm)",
        ylabel="Measured absorbance",
        xlim=(880, 1100),
        ylim=(0, 3.05),
    )
    axis.legend(ncol=3)

    figure.text(
        0.5,
        0.005,
        "Published positions: Pirri et al., Materials 11 (2018) 837.",
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


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")

    absorbance = load_export(ROOT / "YbYag_ABS.csv")
    plot_ybyag_band_comparison(absorbance)


if __name__ == "__main__":
    main()
