"""Depth profile plotter - Streamlit app.

Run locally with:  streamlit run app.py
"""

import hashlib
import io

import matplotlib

matplotlib.use("Agg")  # headless backend; must be set before pyplot is imported
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px
import streamlit as st

from depth_profile import RASTER, parse_csv, process, strongest_isotopes

NONE_OPTION = "None (raw signal)"
TOTAL_OPTION = "Total"
DEFAULT_SELECTED = 8

st.set_page_config(page_title="Depth Profile Plotter", layout="wide")


@st.cache_data(show_spinner=False)
def load_profile(content: bytes):
    return parse_csv(content)


def checkbox_key(file_id: str, label: str) -> str:
    return f"iso::{file_id}::{label}"


def set_all(file_id: str, labels: list[str], value: bool) -> None:
    for label in labels:
        st.session_state[checkbox_key(file_id, label)] = value


def plot_png(result, selected, colors, x_title, y_title, log_y) -> bytes:
    """Render the current plot as a static PNG (matplotlib, no browser needed)."""
    fig, ax = plt.subplots(figsize=(10, 6), dpi=200)
    for label, color in zip(selected, colors):
        ax.plot(result.index, result[label], label=label, color=color, linewidth=1.5)
    if log_y:
        ax.set_yscale("log")
    ax.set_xlabel(x_title)
    ax.set_ylabel(y_title)
    ax.grid(True, which="major", alpha=0.3)
    ax.legend(title="Isotope", loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False,
              ncol=max(1, -(-len(selected) // 30)), fontsize=8)
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight")
    plt.close(fig)
    return buf.getvalue()


# --------------------------------------------------------------------------- #
# File input
# --------------------------------------------------------------------------- #

st.title("Depth Profile Plotter")

uploaded = st.file_uploader(
    "Upload a depth profile CSV (repeated raster or single spot time series)", type=["csv"]
)
if uploaded is None:
    st.info(
        "Upload a CSV to begin. Two formats are recognised automatically:\n\n"
        "- **Repeated raster**: one column per raster pass (`..._rep_N.h5`); the x-axis is raster number.\n"
        "- **Single spot**: a time series whose first row starts with `Time`; the x-axis is time in seconds."
    )
    st.stop()

content = uploaded.getvalue()
file_id = hashlib.md5(content).hexdigest()[:12]
try:
    profile = load_profile(content)
except ValueError as err:
    st.error(f"Could not read **{uploaded.name}**: {err}")
    st.stop()

labels = [iso.label for iso in profile.isotopes]
kind_name = "Repeated raster" if profile.kind == RASTER else "Single spot time series"
st.caption(f"**{uploaded.name}** · {kind_name} · {len(labels)} isotopes · {len(profile.data)} points")

# Default selection for a newly uploaded file: the strongest signals.
if st.session_state.get("loaded_file") != file_id:
    defaults = set(strongest_isotopes(profile, DEFAULT_SELECTED))
    for label in labels:
        st.session_state.setdefault(checkbox_key(file_id, label), label in defaults)
    st.session_state["loaded_file"] = file_id

# --------------------------------------------------------------------------- #
# Sidebar settings
# --------------------------------------------------------------------------- #

# Settings are keyed to the loaded file so a new file starts from defaults and
# the widgets never show values that differ from what is being plotted.
with st.sidebar:
    st.header("Settings")

    normalization = st.selectbox(
        "Normalization",
        [NONE_OPTION, TOTAL_OPTION, *[i.label for i in sorted(profile.isotopes, key=lambda i: (i.mass, i.label))]],
        key=f"norm::{file_id}",
        help=(
            "Divide every signal, point by point, by the chosen isotope's signal. "
            "'Total' divides by the sum of all isotopes (each measured mass counted once)."
        ),
    )

    smoothing = st.slider(
        "Smoothing (moving average window, points)", min_value=1, max_value=100, value=1,
        key=f"smooth::{file_id}",
        help="Centered moving average. 1 = no smoothing.",
    )

    x_values = profile.data.index.tolist()
    fmt = (lambda v: f"{int(v)}") if profile.kind == RASTER else (lambda v: f"{v:.2f} s")
    x_range = st.select_slider(
        f"{profile.x_label} range",
        options=x_values,
        value=(x_values[0], x_values[-1]),
        format_func=fmt,
        key=f"xrange::{file_id}",
        help="Trim the start or end of the profile (e.g. a spike in the first point).",
    )

    use_depth = st.checkbox("Convert x-axis to depth", key=f"use_depth::{file_id}")
    total_depth = None
    if use_depth:
        total_depth = st.number_input(
            "Total depth reached (µm)", min_value=0.0, value=1.0, step=0.1, format="%.3f",
            key=f"depth::{file_id}",
            help="Depth at the last point of the full profile. A constant ablation rate is assumed.",
        )
        if total_depth <= 0:
            total_depth = None

    log_y = st.toggle("Log y-axis", value=True, key=f"log::{file_id}")

# --------------------------------------------------------------------------- #
# Main area: plot + isotope selection
# --------------------------------------------------------------------------- #

plot_col, select_col = st.columns([3, 1])

with select_col:
    st.subheader("Isotopes")
    b1, b2 = st.columns(2)
    b1.button("Select all", on_click=set_all, args=(file_id, labels, True))
    b2.button("Clear all", on_click=set_all, args=(file_id, labels, False))

    # Group isotopes by element, elements ordered by their lightest isotope,
    # with unassigned masses last.
    groups: dict[str, list] = {}
    for iso in sorted(profile.isotopes, key=lambda i: (i.element == "?", i.mass, i.label)):
        groups.setdefault(iso.element, []).append(iso)

    with st.container(height=560):
        for element, isotopes in groups.items():
            st.markdown(f"**{'Unassigned' if element == '?' else element}**")
            cols = st.columns(2)
            for i, iso in enumerate(isotopes):
                cols[i % 2].checkbox(iso.label, key=checkbox_key(file_id, iso.label))

selected = [label for label in labels if st.session_state.get(checkbox_key(file_id, label))]

norm_arg = None if normalization == NONE_OPTION else normalization
result = process(
    profile,
    x_range=x_range,
    smoothing=smoothing,
    normalization=norm_arg,
    total_depth_um=total_depth,
)

if norm_arg is None:
    y_title = "Signal"
elif norm_arg == TOTAL_OPTION:
    y_title = "Fraction of total signal"
else:
    y_title = f"Signal / {norm_arg}"
x_title = result.index.name

with plot_col:
    if not selected:
        st.warning("Select at least one isotope to plot.")
    else:
        palette = px.colors.qualitative.Plotly if len(selected) <= 10 else px.colors.qualitative.Alphabet
        colors = [palette[i % len(palette)] for i in range(len(selected))]
        fig = go.Figure()
        for label, color in zip(selected, colors):
            fig.add_trace(go.Scatter(
                x=result.index, y=result[label], mode="lines", name=label,
                line=dict(width=2, color=color),
            ))
        fig.update_layout(
            height=620,
            xaxis_title=x_title,
            yaxis_title=y_title,
            yaxis_type="log" if log_y else "linear",
            hovermode="x unified",
            legend=dict(title="Isotope"),
            margin=dict(l=10, r=10, t=30, b=10),
        )
        st.plotly_chart(fig, config={
            "toImageButtonOptions": {"format": "png", "filename": "depth_profile", "scale": 3},
            "displaylogo": False,
        })
        if log_y and (result[selected] <= 0).any().any():
            st.caption("Zero or negative values are not shown on the log axis.")

        base_name = uploaded.name.rsplit(".", 1)[0]
        d1, d2, d3 = st.columns(3)
        d1.download_button(
            "Download plot (PNG)",
            plot_png(result, selected, colors, x_title, y_title, log_y),
            file_name=f"{base_name}_plot.png",
            mime="image/png",
        )
        d2.download_button(
            "Download processed data (CSV)",
            result[selected].to_csv().encode("utf-8"),
            file_name=f"{base_name}_processed.csv",
            mime="text/csv",
        )
        d3.download_button(
            "Download interactive plot (HTML)",
            fig.to_html(include_plotlyjs="cdn").encode("utf-8"),
            file_name=f"{base_name}_plot.html",
            mime="text/html",
        )
