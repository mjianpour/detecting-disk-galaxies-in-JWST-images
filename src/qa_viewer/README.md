# QA viewer for the detection / morphology pipeline

A Streamlit app for visually QA-ing pipeline outputs: catching bad deblends,
unstable Sérsic fits, and misclassified sources.

## Run it

```bash
# from the repo root (venv already has streamlit/astropy/numpy/pandas/matplotlib)
source .venv/bin/activate
cd src/qa_viewer
python make_mock_data.py       # once — generates demo data in mock_data/
streamlit run app.py
```

The app opens pointed at the mock data and is fully usable immediately.

## What's in here

| File | Role |
|---|---|
| `app.py` | Streamlit UI only — no FITS/catalog parsing |
| `pipeline_io.py` | All pipeline-output I/O (catalog, cutouts, mosaic) |
| `make_mock_data.py` | Generates `mock_data/` (catalog + cutouts + mosaic) |

## Pointing it at real pipeline outputs

Everything is configured live in the sidebar ("Pipeline outputs"); no code
changes needed if your outputs follow these conventions:

- **Catalog** — CSV or FITS table. `source_id` is required; the app uses
  `x, y, sersic_n, gini, m20, concentration, asymmetry, is_disk, flag,
  flag_sersic` when present and degrades gracefully when they're missing.
- **Cutouts** — one FITS file per source at `{cutouts_dir}/{source_id}.fits`,
  science image in the `SCI` HDU (falls back to the first 2D image HDU).
- **Segmentation masks** — two conventions, chosen in the sidebar:
  - *SEG HDU in same file* (default): integer label map in a `SEG` HDU of the
    cutout file (falls back to the first integer-typed HDU after the science).
  - *Parallel directory*: `{seg_cutouts_dir}/{source_id}.fits`.
  - Labels: pixels equal to `source_id` are treated as the source's own
    segment; other non-zero labels are drawn as neighboring sources (useful
    for deblend QA). A binary 0/1 mask also works.
- **Mosaic** (optional) — full-pointing FITS (`SCI` HDU); automatically
  downsampled for display if large.

If your column names or file layout differ, adapt `pipeline_io.py` — the UI
never touches files directly.

## Views

- **Catalog & detail** — filterable table (Sérsic n / Gini / M20 ranges, disk
  classification, hide-flagged toggle). Selecting a row loads that source's
  cutout with live stretch controls (ZScale or percentile interval; linear /
  asinh / log stretch; asinh softening), a segmentation-boundary overlay, an
  optional label-map panel, all statmorph metrics, both disk-cut criteria
  evaluated explicitly, and a Gini–M20 diagram with the merger/disk locus
  line (`G = −0.14·M20 + 0.33`) and the selected source highlighted.
- **Mosaic overview** — downsampled pointing with sources marked (disk = blue
  circle, non-disk = orange triangle; shape carries the distinction as well as
  color). Markers aren't clickable (that would need plotly/bokeh — kept out to
  stay minimal); use the "Jump to source" selector instead.

FITS reads are cached with `st.cache_data`, keyed on file mtime — so
re-running the pipeline and overwriting an output file is picked up on the
next interaction without restarting the app.
