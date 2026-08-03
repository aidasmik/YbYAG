"""
Reader for J.A. Woollam CompleteEASE Mueller-matrix .dat exports (Yb:YAG, RC2).

The exports in `YbYag txt/` carry a UTF-8 BOM and CRLF line endings, and the
`VASEmethod[...]` header therefore sits on the *second* line -- a naive
`head -1` / first-line header parse silently misses it and mislabels every
transmission file as reflection.  `read_dat` handles both.

File structure (whitespace separated), one repeating block per (wavelength, angle):

    E     <wl_nm> <angle> <Psi> <Delta> <errPsi> <errDelta>
    dPolE <wl_nm> <angle> <depol> <err>
    mm11  <wl_nm> <angle> <value> <sigma>                    # raw intensity R or T
    mm12  <wl_nm> <angle> <value> <value> <sigma> <sigma>    # normalised to mm11
    ...
    mm44  <wl_nm> <angle> <value> <value> <sigma> <sigma>

Note on the measurement type: CompleteEASE writes `TransmissionMeas=1` into the
header for transmission scans and omits the key entirely for reflection.  Several
files in this data set are named `...MMt...` but are actually *reflection*
measurements (e.g. `YbYag_5a_MMt`, `YbYag_5a_UMMt` -- no flag, AOI 30-75 deg), so
the header, not the filename, decides.
"""

import re
import numpy as np
import pandas as pd

MM_KEYS = [f"mm{i}{j}" for i in range(1, 5) for j in range(1, 5)]


def read_dat(path):
    """Parse a Woollam .dat export.

    Returns
    -------
    df : DataFrame, one row per (wl, angle), columns wl, angle, Psi, Delta,
         errPsi, errDelta, depol, mm11 (+ _err) and the normalised mm12..mm44
         (+ _err).
    meta : dict of the VASEmethod header keys, plus
           meta['is_transmission'] -- bool taken from TransmissionMeas.
    """
    rows, meta = {}, {}
    with open(path, "r", encoding="utf-8-sig", errors="replace") as fh:
        for line in fh:
            parts = line.split()
            if not parts:
                continue
            tag = parts[0]
            if tag.startswith("VASEmethod"):
                m = re.search(r"\[(.*)\]", line)
                if m:
                    for kv in m.group(1).split(","):
                        if "=" in kv:
                            k, v = kv.split("=", 1)
                            meta[k.strip()] = v.strip()
                continue
            if tag == "nm":
                meta["units"] = "nm"
                continue
            if len(parts) < 4:
                continue
            try:
                wl, angle = float(parts[1]), float(parts[2])
            except ValueError:
                continue
            d = rows.setdefault((wl, angle), {"wl": wl, "angle": angle})
            vals = parts[3:]
            if tag == "E":
                d["Psi"], d["Delta"] = float(vals[0]), float(vals[1])
                if len(vals) >= 4:
                    d["errPsi"], d["errDelta"] = float(vals[2]), float(vals[3])
            elif tag == "dPolE":
                d["depol"] = float(vals[0])
            elif tag == "mm11":
                d["mm11"] = float(vals[0])
                if len(vals) >= 2:
                    d["mm11_err"] = float(vals[1])
            elif re.fullmatch(r"mm\d\d", tag):
                # normalised element: value value sigma sigma (value is repeated)
                d[tag] = float(vals[0])
                if len(vals) >= 3:
                    d[tag + "_err"] = float(vals[2])
                elif len(vals) >= 2:
                    d[tag + "_err"] = float(vals[1])

    meta["is_transmission"] = str(meta.get("TransmissionMeas", "0")).strip() == "1"
    df = pd.DataFrame(list(rows.values())).sort_values(["angle", "wl"])
    return df.reset_index(drop=True), meta


def grid(df, col, wls, angles):
    """Pivot one column onto a dense [wavelength, angle] array (NaN where absent)."""
    if col not in df.columns:
        return np.full((len(wls), len(angles)), np.nan)
    p = df.pivot_table(index="wl", columns="angle", values=col)
    return p.reindex(index=wls, columns=angles).values
