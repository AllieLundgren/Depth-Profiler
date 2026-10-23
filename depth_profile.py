"""Parsing and processing for depth profile CSV exports.

Two export formats are supported:

* Repeated raster - first header cell is "Non Quantified Element"; one column
  per raster pass named like ``..._rep_12.h5``, followed by summary columns
  (Average Value / Standard Deviation) that are ignored here.
* Single spot time series - first row starts with "Time" and holds the time
  stamps (s); the second row is a header ("Element, Atomic Mass, ...") and each
  following row is one isotope.

Both are parsed into a :class:`Profile`: a DataFrame indexed by x (raster number
or time) with one column per isotope.
"""

from __future__ import annotations

import csv
import io
import re
from dataclasses import dataclass

import numpy as np
import pandas as pd

RASTER = "raster"
SPOT = "spot"

_REP_RE = re.compile(r"rep_(\d+)", re.IGNORECASE)

# Rows at the same mass whose signals differ by less than this (median relative
# difference) are treated as the same measured peak reported under several
# element assignments, and are counted once in the "Total" normalization.
DUPLICATE_TOLERANCE = 0.05


@dataclass
class Isotope:
    label: str    # unique display label, e.g. "57Fe" or "?-175"
    element: str  # element symbol, or "?" if unassigned
    mass: int


@dataclass
class Profile:
    kind: str                 # RASTER or SPOT
    data: pd.DataFrame        # index = x (raster number or time), columns = isotope labels
    isotopes: list[Isotope]   # same order as data.columns

    @property
    def x_label(self) -> str:
        return "Raster number" if self.kind == RASTER else "Time (s)"


# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #

def _to_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def _to_mass(value: str) -> int | None:
    """Parse a mass number; None unless it is a positive whole number."""
    try:
        mass = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(mass) or mass < 1 or abs(mass - round(mass)) > 1e-6:
        return None
    return int(round(mass))


def _is_blank(row: list[str]) -> bool:
    return not any(cell.strip() for cell in row)


def _data_block(rows: list[list[str]]) -> list[list[str]]:
    """Rows of the first data block: from the first non-blank row up to the next
    blank row. Anything after that (e.g. tables added in a spreadsheet) is ignored."""
    block: list[list[str]] = []
    for row in rows:
        if _is_blank(row):
            if block:
                break
            continue
        block.append(row)
    return block


def _make_label(element: str, mass: int) -> str:
    if element in ("", "?"):
        return f"?-{mass}"
    return f"{mass}{element}"


def _build_profile(kind: str, x: list[float], rows: list[tuple[str, int, list[float]]]) -> Profile:
    isotopes: list[Isotope] = []
    columns: dict[str, list[float]] = {}
    for element, mass, values in rows:
        label = _make_label(element, mass)
        base, n = label, 2
        while label in columns:  # keep labels unique if the export repeats a row
            label = f"{base} ({n})"
            n += 1
        isotopes.append(Isotope(label=label, element=element or "?", mass=mass))
        columns[label] = values
    data = pd.DataFrame(columns, index=pd.Index(x, name="x"))
    return Profile(kind=kind, data=data, isotopes=isotopes)


def _parse_raster(rows: list[list[str]]) -> Profile:
    header = rows[0]
    rep_cols = []  # (column index, rep number)
    for i, name in enumerate(header):
        m = _REP_RE.search(name)
        if m:
            rep_cols.append((i, int(m.group(1))))
    if not rep_cols:
        raise ValueError("No raster columns (named like '..._rep_0.h5') were found in the header.")

    parsed = []
    for row in rows[1:]:
        if len(row) < 2:
            continue
        mass = _to_mass(row[1])
        if mass is None:
            continue
        values = [_to_float(row[i]) if i < len(row) else np.nan for i, _ in rep_cols]
        parsed.append((row[0].strip(), mass, values))
    return _build_profile(RASTER, [rep for _, rep in rep_cols], parsed)


def _parse_spot(rows: list[list[str]]) -> Profile:
    times = [_to_float(v) for v in rows[0][2:]]
    # Drop trailing empty cells.
    while times and np.isnan(times[-1]):
        times.pop()
    if not times:
        raise ValueError("The 'Time' row contains no time values.")

    n = len(times)
    parsed = []
    for row in rows[1:]:
        if len(row) < 3:
            continue
        mass = _to_mass(row[1])
        if mass is None:  # skips the "Element, Atomic Mass, ..." header row
            continue
        values = [_to_float(v) for v in row[2:2 + n]]
        values += [np.nan] * (n - len(values))
        parsed.append((row[0].strip(), mass, values))
    return _build_profile(SPOT, times, parsed)


def parse_csv(content: bytes | str) -> Profile:
    """Detect the file format and parse it into a Profile."""
    if isinstance(content, bytes):
        content = content.decode("utf-8-sig", errors="replace")
    rows = _data_block(list(csv.reader(io.StringIO(content))))
    if not rows:
        raise ValueError("The file is empty.")

    first = rows[0][0].strip().lower()
    if first == "time":
        profile = _parse_spot(rows)
    elif any(_REP_RE.search(cell) for cell in rows[0]):
        profile = _parse_raster(rows)
    else:
        raise ValueError(
            "Unrecognised file format. Expected either a repeated raster export "
            "(columns named '..._rep_N.h5') or a spot time series export (first row starting with 'Time')."
        )
    if profile.data.empty or not profile.isotopes:
        raise ValueError("No isotope rows were found in the file.")
    return profile


# --------------------------------------------------------------------------- #
# Processing
# --------------------------------------------------------------------------- #

def strongest_isotopes(profile: Profile, n: int = 8) -> list[str]:
    """Labels of the n isotopes with the highest mean signal."""
    return profile.data.mean().sort_values(ascending=False).head(n).index.tolist()


def _same_signal(a: pd.Series, b: pd.Series) -> bool:
    scale = np.maximum(a.abs(), b.abs())
    rel = ((a - b).abs() / scale.where(scale > 0)).dropna()
    return not rel.empty and rel.median() < DUPLICATE_TOLERANCE


def total_signal(profile: Profile, data: pd.DataFrame | None = None) -> pd.Series:
    """Sum of all isotope signals at each x, counting each measured mass once.

    Rows that share a mass and carry (nearly) the same signal - e.g. 50Cr and
    50V reported from the same peak - are counted once. Rows at the same mass
    with clearly different signals are summed.
    """
    data = profile.data if data is None else data
    by_mass: dict[int, list[str]] = {}
    for iso in profile.isotopes:
        by_mass.setdefault(iso.mass, []).append(iso.label)

    kept: list[str] = []
    for labels in by_mass.values():
        unique: list[str] = []
        for label in labels:
            if not any(_same_signal(data[label], data[u]) for u in unique):
                unique.append(label)
        kept.extend(unique)
    return data[kept].sum(axis=1, min_count=1)


def depth_axis(x: pd.Index, total_depth_um: float) -> pd.Index:
    """Map x linearly onto 0..total_depth_um, assuming a constant ablation rate."""
    x0, x1 = float(x.min()), float(x.max())
    if x1 == x0:
        return pd.Index(np.zeros(len(x)), name="Depth (µm)")
    return pd.Index((x.to_numpy(dtype=float) - x0) / (x1 - x0) * total_depth_um, name="Depth (µm)")


def process(
    profile: Profile,
    x_range: tuple[float, float] | None = None,
    smoothing: int = 1,
    normalization: str | None = None,
    total_depth_um: float | None = None,
) -> pd.DataFrame:
    """Return the processed signals for every isotope.

    Steps, in order: depth conversion (calibrated on the full, untrimmed
    range), trimming to ``x_range`` (in original x units), centered moving
    average smoothing over ``smoothing`` points, then normalization.

    ``normalization`` is None for raw signal, "Total", or an isotope label.
    """
    data = profile.data.copy()

    if total_depth_um:
        depth = depth_axis(data.index, total_depth_um)
    else:
        depth = None

    if x_range is not None:
        mask = (data.index >= x_range[0]) & (data.index <= x_range[1])
        data = data[mask]
        if depth is not None:
            depth = depth[mask]

    if smoothing and smoothing > 1:
        data = data.rolling(window=int(smoothing), center=True, min_periods=1).mean()

    if normalization == "Total":
        denom = total_signal(profile, data)
        data = data.div(denom.where(denom != 0), axis=0)
    elif normalization:
        denom = data[normalization]
        data = data.div(denom.where(denom != 0), axis=0)

    if depth is not None:
        data.index = depth
    else:
        data.index.name = profile.x_label
    return data
