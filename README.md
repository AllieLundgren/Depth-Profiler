# Depth Profile Plotter

A Streamlit app for plotting depth profiles from CSV exports.

## Supported files

The format is detected automatically:

| Format | How it is recognised | x-axis |
|---|---|---|
| **Repeated raster** | Header has one column per pass, named like `..._rep_12.h5` (summary columns such as *Average Value* / *Standard Deviation* are ignored) | Raster number |
| **Single spot time series** | First row starts with `Time` and holds the time stamps; next row is `Element, Atomic Mass, ...` | Time (s) |

Each data row is one isotope: element symbol, mass, then the signal at each point. Rows with an unassigned element (`?`) are labelled like `?-175`.

Only the first block of data is read: everything after the first blank row (e.g. tables added later in a spreadsheet) is ignored, as are rows whose mass is not a whole number.

## Features

- **Isotope selection**: checkboxes grouped by element, with *Select all* / *Clear all*. The 8 strongest signals are selected when a file is loaded.
- **Normalization**:
  - *None*: raw signal.
  - *Total*: each point divided by the sum of all isotopes at that point. Rows at the same mass that carry the same signal (e.g. 50Cr and 50V reported from one peak, median difference < 5%) are counted once; rows at the same mass with different signals are summed.
  - *Any isotope*: each point divided by that isotope's signal at the same point (e.g. normalize to 57Fe).
- **Smoothing**: centered moving average, window 1–25 points (1 = off).
- **Range trim**: drop points at the start or end (e.g. a spike in the first point).
- **Depth conversion**: enter the total depth reached (µm) to convert the x-axis to depth. A constant ablation rate is assumed: the first point of the full profile is 0 µm and the last is the depth entered. Trimming does not change this calibration.
- **Log / linear** y-axis.
- **Downloads**: the plot as a high-resolution PNG, processed data (CSV) for the selected isotopes, and the interactive plot (HTML).

Processing order: depth calibration → trim → smoothing → normalization.

## Run locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Deploy on Streamlit Community Cloud

1. Push this repository to GitHub.
2. Go to <https://share.streamlit.io>, choose **Create app**, and select the repository.
3. Set the main file to `app.py` and deploy.

## Files

- `app.py`: Streamlit user interface.
- `depth_profile.py`: CSV parsing and processing (format detection, normalization, smoothing, depth conversion).
- `requirements.txt`: Python dependencies.

CSV files are excluded by `.gitignore` so measurement data is not committed by accident.
