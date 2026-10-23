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

## Multiple spots

Upload up to 12 single spot files together (repeated raster files are loaded one at a time). A **View** selector appears in the sidebar:

- **Average**: the mean of all included spots for each selected isotope, with an optional shaded ±1 SD band (sample standard deviation).
- **Single spot**: one spot on its own, chosen from a dropdown.
- **Compare spots**: one isotope, with every spot as a thin line and the average as a thick black line. Excluded spots are shown dotted grey.

How the average is calculated:

- Each spot is processed on its own first (trim, smoothing, normalization), then the spots are averaged point by point. Normalizing before averaging means the average is the typical ratio across spots, and no single spot with a stronger overall signal dominates it.
- Spots are matched by time stamp (all spots are assumed to start at t = 0). If they have different lengths, the average covers only the range that every spot shares.
- With depth conversion on, every spot is assumed to reach the same total depth.
- An isotope missing from some files is averaged over the files that contain it; the caption lists such isotopes.
- Untick a spot under **Spots included in the average** to leave it out without re-uploading.

The CSV download contains the mean and SD columns (Average), the processed spot (Single spot), or every spot plus the mean and SD for the chosen isotope (Compare spots).

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
- `depth_profile.py`: CSV parsing and processing (format detection, normalization, smoothing, depth conversion, multi-spot averaging).
- `requirements.txt`: Python dependencies.

CSV files are excluded by `.gitignore` so measurement data is not committed by accident.
