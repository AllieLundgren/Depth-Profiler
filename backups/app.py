"""Depth profile plotter - Streamlit app.

Run locally with:  streamlit run app.py
"""

import hashlib
import io

import matplotlib

matplotlib.use("Agg")  # headless backend; must be set before pyplot is imported
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import streamlit as st

from depth_profile import (
    RASTER,
    SPOT,
    align_results,
    average_results,
    parse_csv,
    process,
    strongest_isotopes,
    union_isotopes,
)

NONE_OPTION = "None (raw signal)"
TOTAL_OPTION = "Total"
DEFAULT_SELECTED = 8
MAX_SPOTS = 12

VIEW_AVERAGE = "Average"
VIEW_SINGLE = "Single spot"
VIEW_COMPARE = "Compare spots"

st.set_page_config(page_title="Depth Profile Plotter", layout="wide")


@st.cache_data(show_spinner=False)
def load_profile(content: bytes):
    return parse_csv(content)


def checkbox_key(set_id: str, label: str) -> str:
    return f"iso::{set_id}::{label}"


def set_all(set_id: str, labels: list[str], value: bool) -> None:
    for label in labels:
        st.session_state[checkbox_key(set_id, label)] = value


def palette_for(n: int) -> list[str]:
    palette = px.colors.qualitative.Plotly if n <= 10 else px.colors.qualitative.Alphabet
    return [palette[i % len(palette)] for i in range(n)]


def rgba(hex_color: str, alpha: float) -> str:
    h = hex_color.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    return f"rgba({r},{g},{b},{alpha})"


def band_limits(mean: pd.Series, sd: pd.Series, log_y: bool) -> tuple[pd.Series, pd.Series]:
    lower, upper = mean - sd, mean + sd
    if log_y:  # a log axis cannot show values <= 0
        lower = lower.where(lower > 0)
    return lower, upper


# A trace is a dict: label, x, y, color, width, dash ("solid"/"dot"), and
# optionally lower/upper for a shaded band. The same list drives the
# interactive Plotly chart and the static PNG.

def plotly_figure(traces, x_title, y_title, log_y, legend_title) -> go.Figure:
    fig = go.Figure()
    for t in traces:
        if t.get("lower") is not None:
            fig.add_trace(go.Scatter(
                x=t["x"], y=t["upper"], mode="lines", line=dict(width=0),
                showlegend=False, hoverinfo="skip", legendgroup=t["label"],
            ))
            fig.add_trace(go.Scatter(
                x=t["x"], y=t["lower"], mode="lines", line=dict(width=0),
                fill="tonexty", fillcolor=rgba(t["color"], 0.18),
                showlegend=False, hoverinfo="skip", legendgroup=t["label"],
            ))
        fig.add_trace(go.Scatter(
            x=t["x"], y=t["y"], mode="lines", name=t["label"], legendgroup=t["label"],
            line=dict(width=t["width"], color=t["color"], dash=t["dash"]),
        ))
    fig.update_layout(
        height=620,
        xaxis_title=x_title,
        yaxis_title=y_title,
        yaxis_type="log" if log_y else "linear",
        hovermode="x unified",
        legend=dict(title=legend_title),
        margin=dict(l=10, r=10, t=30, b=10),
    )
    return fig


def plot_png(traces, x_title, y_title, log_y, legend_title) -> bytes:
    """Render the current plot as a static PNG (matplotlib, no browser needed)."""
    fig, ax = plt.subplots(figsize=(10, 6), dpi=200)
    for t in traces:
        if t.get("lower") is not None:
            ax.fill_between(t["x"], t["lower"], t["upper"], color=t["color"], alpha=0.18, linewidth=0)
        ax.plot(t["x"], t["y"], label=t["label"], color=t["color"], linewidth=t["width"] * 0.75,
                linestyle=":" if t["dash"] == "dot" else "-")
    if log_y:
        ax.set_yscale("log")
    ax.set_xlabel(x_title)
    ax.set_ylabel(y_title)
    ax.grid(True, which="major", alpha=0.3)
    ax.legend(title=legend_title, loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False,
              ncol=max(1, -(-len(traces) // 30)), fontsize=8)
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight")
    plt.close(fig)
    return buf.getvalue()


# --------------------------------------------------------------------------- #
# File input
# --------------------------------------------------------------------------- #

st.title("Depth Profile Plotter")

uploads = st.file_uploader(
    "Upload a depth profile CSV: one repeated raster file, or up to "
    f"{MAX_SPOTS} single spot time series files",
    type=["csv"],
    accept_multiple_files=True,
)
if not uploads:
    st.info(
        "Upload a CSV to begin. Two formats are recognised automatically:\n\n"
        "- **Repeated raster**: one column per raster pass (`..._rep_N.h5`); the x-axis is raster number.\n"
        "- **Single spot**: a time series whose first row starts with `Time`; the x-axis is time in seconds. "
        f"Upload several spot files together (up to {MAX_SPOTS}) to see each spot and their average."
    )
    st.stop()

profiles, names, hashes, errors = [], [], [], []
for up in uploads:
    content = up.getvalue()
    try:
        profiles.append(load_profile(content))
        names.append(up.name)
        hashes.append(hashlib.md5(content).hexdigest())
    except ValueError as err:
        errors.append(f"**{up.name}**: {err}")
for msg in errors:
    st.error(f"Could not read {msg}")
if not profiles:
    st.stop()

kinds = {p.kind for p in profiles}
if len(profiles) > 1 and RASTER in kinds:
    st.error(
        "Repeated raster files can only be loaded one at a time. "
        "Multiple files are supported for single spot time series only."
    )
    st.stop()
if len(profiles) > MAX_SPOTS:
    st.error(f"Too many files: {len(profiles)} uploaded, the maximum is {MAX_SPOTS} spots.")
    st.stop()

multi = len(profiles) > 1
kind = profiles[0].kind
x_label = profiles[0].x_label
set_id = hashlib.md5("".join(sorted(hashes)).encode()).hexdigest()[:12]

isotopes = union_isotopes(profiles)
labels = [iso.label for iso in isotopes]
by_mass = [iso.label for iso in sorted(isotopes, key=lambda i: (i.mass, i.label))]

if multi:
    longest = max(len(p.data) for p in profiles)
    st.caption(f"{len(profiles)} single spot files · {len(labels)} isotopes · up to {longest} points each")
else:
    kind_name = "Repeated raster" if kind == RASTER else "Single spot time series"
    st.caption(f"**{names[0]}** · {kind_name} · {len(labels)} isotopes · {len(profiles[0].data)} points")

# Default selection for a newly uploaded file set: the strongest signals.
if st.session_state.get("loaded_set") != set_id:
    defaults = set(strongest_isotopes(profiles, DEFAULT_SELECTED))
    for label in labels:
        st.session_state.setdefault(checkbox_key(set_id, label), label in defaults)
    st.session_state["loaded_set"] = set_id

# Streamlit forgets the state of widgets that are not drawn in a run. The
# isotope checkboxes are hidden in the Compare view, so re-store their values
# each run to keep the selection when switching views.
for label in labels:
    key = checkbox_key(set_id, label)
    if key in st.session_state:
        st.session_state[key] = st.session_state[key]

# --------------------------------------------------------------------------- #
# Sidebar settings
# --------------------------------------------------------------------------- #

# Settings are keyed to the loaded file set so new files start from defaults
# and the widgets never show values that differ from what is being plotted.
with st.sidebar:
    st.header("Settings")

    view = VIEW_SINGLE
    included = [0]
    if multi:
        view = st.radio("View", [VIEW_AVERAGE, VIEW_SINGLE, VIEW_COMPARE], key=f"view::{set_id}")
        if view == VIEW_SINGLE:
            spot = st.selectbox("Spot", range(len(names)), format_func=lambda i: names[i],
                                key=f"spot::{set_id}")
            included = [spot]
        else:
            with st.expander("Spots included in the average", expanded=False):
                included = [
                    i for i, name in enumerate(names)
                    if st.checkbox(name, value=True, key=f"incl::{set_id}::{i}")
                ]
        if view == VIEW_AVERAGE:
            show_sd = st.toggle("Show ±1 SD band", value=True, key=f"sd::{set_id}")
        if view == VIEW_COMPARE:
            compare_iso = st.selectbox(
                "Isotope to compare", by_mass, key=f"cmp::{set_id}",
                index=by_mass.index(strongest_isotopes(profiles, 1)[0]),
            )
        st.divider()

    normalization = st.selectbox(
        "Normalization",
        [NONE_OPTION, TOTAL_OPTION, *by_mass],
        key=f"norm::{set_id}",
        help=(
            "Divide every signal, point by point, by the chosen isotope's signal. "
            "'Total' divides by the sum of all isotopes (each measured mass counted once). "
            "With several spots, each spot is normalized before averaging."
        ),
    )

    smoothing = st.slider(
        "Smoothing (moving average window, points)", min_value=1, max_value=25, value=1,
        key=f"smooth::{set_id}",
        help="Centered moving average. 1 = no smoothing.",
    )

    # Range options come from the longest file so every spot can be trimmed.
    x_values = max((p.data.index for p in profiles), key=len).tolist()
    fmt = (lambda v: f"{int(v)}") if kind == RASTER else (lambda v: f"{v:.2f} s")
    x_range = st.select_slider(
        f"{x_label} range",
        options=x_values,
        value=(x_values[0], x_values[-1]),
        format_func=fmt,
        key=f"xrange::{set_id}",
        help="Trim the start or end of the profile (e.g. a spike in the first point).",
    )

    use_depth = st.checkbox("Convert x-axis to depth", key=f"use_depth::{set_id}")
    total_depth = None
    if use_depth:
        total_depth = st.number_input(
            "Total depth reached (µm)", min_value=0.0, value=1.0, step=0.1, format="%.3f",
            key=f"depth::{set_id}",
            help=(
                "Depth at the last point of the full profile. A constant ablation rate is assumed"
                + (", and every spot is assumed to reach this depth." if multi else ".")
            ),
        )
        if total_depth <= 0:
            total_depth = None

    log_y = st.toggle("Log y-axis", value=True, key=f"log::{set_id}")

# --------------------------------------------------------------------------- #
# Main area: plot + isotope selection
# --------------------------------------------------------------------------- #

plot_col, select_col = st.columns([3, 1])

if view != VIEW_COMPARE:
    with select_col:
        st.subheader("Isotopes")
        b1, b2 = st.columns(2)
        b1.button("Select all", on_click=set_all, args=(set_id, labels, True))
        b2.button("Clear all", on_click=set_all, args=(set_id, labels, False))

        # Group isotopes by element, elements ordered by their lightest isotope,
        # with unassigned masses last.
        groups: dict[str, list] = {}
        for iso in sorted(isotopes, key=lambda i: (i.element == "?", i.mass, i.label)):
            groups.setdefault(iso.element, []).append(iso)

        with st.container(height=560):
            for element, group in groups.items():
                st.markdown(f"**{'Unassigned' if element == '?' else element}**")
                cols = st.columns(2)
                for i, iso in enumerate(group):
                    cols[i % 2].checkbox(iso.label, key=checkbox_key(set_id, iso.label))

selected = [label for label in labels if st.session_state.get(checkbox_key(set_id, label))]

norm_arg = None if normalization == NONE_OPTION else normalization
if norm_arg is None:
    y_title = "Signal"
elif norm_arg == TOTAL_OPTION:
    y_title = "Fraction of total signal"
else:
    y_title = f"Signal / {norm_arg}"


def processed(i: int) -> pd.DataFrame:
    return process(profiles[i], x_range=x_range, smoothing=smoothing,
                   normalization=norm_arg, total_depth_um=total_depth)


with plot_col:
    if view != VIEW_SINGLE and not included:
        st.warning("Include at least one spot (sidebar: Spots included in the average).")
        st.stop()
    if view != VIEW_COMPARE and not selected:
        st.warning("Select at least one isotope to plot.")
        st.stop()

    base_name = names[0].rsplit(".", 1)[0] if not multi else "spots"
    legend_title = "Isotope"

    if view == VIEW_SINGLE:
        i = included[0]
        result = processed(i)
        base_name = names[i].rsplit(".", 1)[0]
        traces = [
            dict(label=label, x=result.index, y=result[label], color=color, width=2, dash="solid")
            for label, color in zip(selected, palette_for(len(selected)))
            if label in result.columns
        ]
        table = result[[t["label"] for t in traces]]

    else:
        results = {i: processed(i) for i in range(len(profiles))}
        cols_needed = selected if view == VIEW_AVERAGE else [compare_iso]
        aligned_all = align_results([results[i] for i in range(len(profiles))], cols_needed)
        mean, sd, n = average_results([aligned_all[i] for i in included])
        x = mean.index

        if view == VIEW_AVERAGE:
            base_name = "spots_average"
            traces = []
            for label, color in zip(selected, palette_for(len(selected))):
                t = dict(label=label, x=x, y=mean[label], color=color, width=2, dash="solid")
                if show_sd and len(included) > 1:
                    t["lower"], t["upper"] = band_limits(mean[label], sd[label], log_y)
                traces.append(t)
            table = pd.concat(
                {f"{label} mean": mean[label] for label in selected}
                | {f"{label} SD": sd[label] for label in selected},
                axis=1,
            )[[f"{label} {s}" for label in selected for s in ("mean", "SD")]]
            missing = [f"{label} ({n[label]}/{len(included)})" for label in selected
                       if n[label] < len(included)]
            st.caption(
                f"Average of {len(included)} of {len(profiles)} spots · "
                f"common range {x.min():.3g} – {x.max():.3g} ({x.name})"
                + (f" · present in only some spots: {', '.join(missing)}" if missing else "")
            )

        else:  # VIEW_COMPARE
            base_name = f"spots_compare_{compare_iso}"
            legend_title = compare_iso
            colors = palette_for(len(profiles))
            traces = []
            for i, name in enumerate(names):
                r = results[i]
                if compare_iso not in r.columns:
                    continue
                short = name.rsplit(".", 1)[0]
                excluded = i not in included
                traces.append(dict(
                    label=f"{short} (excluded)" if excluded else short,
                    x=r.index, y=r[compare_iso], color="#9e9e9e" if excluded else colors[i],
                    width=1.2, dash="dot" if excluded else "solid",
                ))
            traces.append(dict(label=f"Average ({len(included)})", x=x, y=mean[compare_iso],
                               color="#000000", width=3.5, dash="solid"))
            table = pd.concat(
                {names[i].rsplit(".", 1)[0]: aligned_all[i][compare_iso] for i in range(len(profiles))}
                | {"mean": mean[compare_iso], "SD": sd[compare_iso]},
                axis=1,
            )
            st.caption(
                f"{compare_iso}: each spot and the average of {len(included)} included spots. "
                "Excluded spots are dotted grey."
            )

    x_title = traces[0]["x"].name if traces else x_label
    fig = plotly_figure(traces, x_title, y_title, log_y, legend_title)
    st.plotly_chart(fig, config={
        "toImageButtonOptions": {"format": "png", "filename": base_name, "scale": 3},
        "displaylogo": False,
    })
    if log_y and any((np.asarray(t["y"], dtype=float) <= 0).any() for t in traces):
        st.caption("Zero or negative values are not shown on the log axis.")

    d1, d2, d3 = st.columns(3)
    d1.download_button(
        "Download plot (PNG)",
        plot_png(traces, x_title, y_title, log_y, legend_title),
        file_name=f"{base_name}_plot.png",
        mime="image/png",
    )
    d2.download_button(
        "Download processed data (CSV)",
        table.to_csv().encode("utf-8"),
        file_name=f"{base_name}_processed.csv",
        mime="text/csv",
    )
    d3.download_button(
        "Download interactive plot (HTML)",
        fig.to_html(include_plotlyjs="cdn").encode("utf-8"),
        file_name=f"{base_name}_plot.html",
        mime="text/html",
    )
