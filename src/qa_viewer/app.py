"""Streamlit QA viewer for the JWST disk-galaxy detection/morphology pipeline.

Run from this directory:  streamlit run app.py
All pipeline-output reading goes through pipeline_io.py; this file is UI only.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st
from astropy.visualization import (
    AsinhStretch, ImageNormalize, LinearStretch, LogStretch,
    PercentileInterval, ZScaleInterval,
)

from pipeline_io import DataConfig, load_catalog, load_cutout, load_mosaic

# --- chart styling (validated palette; see dataviz notes in README) ---------
COLOR_DISK = "#2a78d6"      # categorical slot 1 (blue)
COLOR_NONDISK = "#eb6834"   # categorical slot 2 (orange)
COLOR_SEG_OWN = "#2a78d6"   # boundary of the selected source's own label
COLOR_SEG_OTHER = "#eb6834" # boundaries of neighboring labels (deblend QA)
SURFACE, INK, INK_2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e1e0d9"

_HERE = os.path.dirname(os.path.abspath(__file__))
MOCK_DIR = os.path.join(_HERE, "mock_data")
DEFAULT_MOSAIC = os.path.normpath(os.path.join(
    _HERE, "..", "..", "images",
    "hlsp_ceers_jwst_nircam_nircam7_f444w_dr0.6_i2d.fits"))

METRIC_COLS = ["sersic_n", "gini", "m20", "concentration", "asymmetry"]
FLAG_COLS = ["flag", "flag_sersic"]

st.set_page_config(page_title="JWST morphology QA", layout="wide")


# --- cached I/O wrappers (mtime busts the cache while the pipeline iterates) -
@st.cache_data(show_spinner=False)
def cached_catalog(path: str, mtime: float) -> pd.DataFrame:
    return load_catalog(path)


@st.cache_data(show_spinner=False)
def cached_cutout(config: DataConfig, source_id, mtime: float):
    return load_cutout(config, source_id)


@st.cache_data(show_spinner="Loading mosaic…")
def cached_mosaic(path: str, mtime: float, max_size: int = 1500):
    return load_mosaic(path, max_size)


def mtime_or_zero(path: str) -> float:
    try:
        return os.path.getmtime(path)
    except OSError:
        return 0.0


# --- sidebar: data locations + display options ------------------------------
def sidebar_config() -> DataConfig:
    st.sidebar.header("Pipeline outputs")
    catalog_path = st.sidebar.text_input(
        "Catalog (CSV or FITS table)", os.path.join(MOCK_DIR, "catalog.csv"))
    cutouts_dir = st.sidebar.text_input(
        "Cutouts directory", os.path.join(MOCK_DIR, "cutouts"))
    convention = st.sidebar.radio(
        "Segmentation mask convention",
        ["same_file", "parallel_dir"],
        format_func=lambda c: {"same_file": "SEG HDU in same file",
                               "parallel_dir": "Parallel directory"}[c],
        help="Where the per-source segmentation mask lives relative to the "
             "science cutout.")
    seg_dir = None
    if convention == "parallel_dir":
        seg_dir = st.sidebar.text_input(
            "Seg cutouts directory", os.path.join(MOCK_DIR, "seg_cutouts"))
    mosaic_path = st.sidebar.text_input(
        "Full mosaic (optional)",
        DEFAULT_MOSAIC if os.path.exists(DEFAULT_MOSAIC)
        else os.path.join(MOCK_DIR, "mosaic.fits"))
    return DataConfig(catalog_path=catalog_path, cutouts_dir=cutouts_dir,
                      mask_convention=convention, seg_cutouts_dir=seg_dir,
                      mosaic_path=mosaic_path or None)


def sidebar_display_options() -> dict:
    st.sidebar.header("Image display")
    interval = st.sidebar.selectbox("Interval", ["ZScale", "Percentile"])
    percentile = st.sidebar.slider(
        "Percentile clip", 90.0, 100.0, 99.5, 0.1,
        disabled=(interval != "Percentile"))
    stretch = st.sidebar.selectbox("Stretch", ["Asinh", "Linear", "Log"])
    asinh_a = st.sidebar.slider(
        "Asinh softening (a)", 0.01, 1.0, 0.1, 0.01,
        disabled=(stretch != "Asinh"))
    cmap = st.sidebar.selectbox("Colormap", ["gray", "viridis", "magma"])
    return dict(interval=interval, percentile=percentile, stretch=stretch,
                asinh_a=asinh_a, cmap=cmap)


def make_norm(data: np.ndarray, opts: dict) -> ImageNormalize:
    finite = data[np.isfinite(data)]
    interval = (ZScaleInterval() if opts["interval"] == "ZScale"
                else PercentileInterval(opts["percentile"]))
    stretch = {"Linear": LinearStretch(),
               "Asinh": AsinhStretch(opts["asinh_a"]),
               "Log": LogStretch()}[opts["stretch"]]
    vmin, vmax = interval.get_limits(finite)
    return ImageNormalize(vmin=vmin, vmax=vmax, stretch=stretch, clip=True)


# --- catalog filtering ------------------------------------------------------
def range_filter(df, col, label, step=0.01):
    if col not in df.columns:
        return df
    lo, hi = float(df[col].min()), float(df[col].max())
    if not np.isfinite(lo) or lo == hi:
        return df
    sel = st.slider(label, lo, hi, (lo, hi), step=step)
    return df[df[col].between(*sel) | df[col].isna()]


def filter_catalog(catalog: pd.DataFrame) -> pd.DataFrame:
    with st.expander("Filters", expanded=True):
        c1, c2, c3, c4 = st.columns(4)
        df = catalog
        with c1:
            df = range_filter(df, "sersic_n", "Sérsic n")
        with c2:
            df = range_filter(df, "gini", "Gini")
        with c3:
            df = range_filter(df, "m20", "M20")
        with c4:
            if "is_disk" in df.columns:
                disk_sel = st.radio("Classification",
                                    ["All", "Disks", "Non-disks"])
                if disk_sel == "Disks":
                    df = df[df["is_disk"]]
                elif disk_sel == "Non-disks":
                    df = df[~df["is_disk"]]
            present_flags = [c for c in FLAG_COLS if c in df.columns]
            if present_flags and st.checkbox("Hide flagged fits"):
                mask = np.zeros(len(df), dtype=bool)
                for c in present_flags:
                    mask |= (df[c].fillna(0) > 0).to_numpy()
                df = df[~mask]
    return df


# --- detail view ------------------------------------------------------------
def plot_cutout(sci, seg, source_id, opts, overlay, show_mask_panel):
    ncols = 2 if (show_mask_panel and seg is not None) else 1
    fig, axes = plt.subplots(1, ncols, figsize=(5 * ncols, 5))
    axes = np.atleast_1d(axes)
    fig.patch.set_facecolor(SURFACE)

    norm = make_norm(sci, opts)
    ax = axes[0]
    ax.imshow(sci, origin="lower", cmap=opts["cmap"], norm=norm)
    ax.set_title(f"Source {source_id} — science", color=INK, fontsize=11)

    if overlay and seg is not None:
        own = (seg == source_id)
        others = (seg > 0) & ~own
        if own.any():
            ax.contour(own, levels=[0.5], colors=COLOR_SEG_OWN, linewidths=1.4)
        if others.any():
            ax.contour(others, levels=[0.5], colors=COLOR_SEG_OTHER,
                       linewidths=1.2, linestyles="--")

    if ncols == 2:
        axm = axes[1]
        seg_show = np.ma.masked_equal(seg, 0)
        axm.imshow(sci, origin="lower", cmap="gray", norm=norm)
        axm.imshow(seg_show, origin="lower", cmap="tab20", alpha=0.5,
                   interpolation="nearest")
        axm.set_title("segmentation labels", color=INK, fontsize=11)

    for a in axes:
        a.set_facecolor(SURFACE)
        a.tick_params(colors=INK_2, labelsize=8)
        for spine in a.spines.values():
            spine.set_color(GRID)
    fig.tight_layout()
    return fig


def metrics_panel(row: pd.Series):
    st.subheader(f"Source {row['source_id']}")

    if "is_disk" in row.index:
        st.markdown("**Disk**" if row["is_disk"] else "**Not a disk**")

    cols = st.columns(2)
    for i, col in enumerate([c for c in METRIC_COLS if c in row.index]):
        cols[i % 2].metric(col, f"{row[col]:.3f}")

    # evaluate the two classification criteria explicitly for QA
    if {"sersic_n", "gini", "m20"}.issubset(row.index):
        n_ok = row["sersic_n"] < 2.5
        locus = -0.14 * row["m20"] + 0.33
        g_ok = row["gini"] < locus
        st.caption(
            f"{'✅' if n_ok else '❌'} Sérsic n = {row['sersic_n']:.2f} "
            f"{'<' if n_ok else '≥'} 2.5  \n"
            f"{'✅' if g_ok else '❌'} Gini = {row['gini']:.3f} "
            f"{'<' if g_ok else '≥'} −0.14·M20 + 0.33 = {locus:.3f}")

    present_flags = [c for c in FLAG_COLS if c in row.index]
    if present_flags:
        st.markdown("**Fit quality**")
        for c in present_flags:
            val = int(row[c]) if np.isfinite(row[c]) else -1
            icon = "✅" if val == 0 else "⚠️"
            st.text(f"{icon} {c} = {val}")


def gini_m20_figure(catalog: pd.DataFrame, selected_id):
    fig, ax = plt.subplots(figsize=(5.5, 4))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    has_disk = "is_disk" in catalog.columns
    disks = catalog[catalog["is_disk"]] if has_disk else catalog.iloc[0:0]
    others = catalog[~catalog["is_disk"]] if has_disk else catalog
    # shape doubles as the non-color channel: disks ●, non-disks ▲
    ax.scatter(others["m20"], others["gini"], s=34, marker="^",
               facecolors=COLOR_NONDISK, edgecolors=SURFACE, linewidths=0.8,
               label="non-disk", zorder=3)
    ax.scatter(disks["m20"], disks["gini"], s=34, marker="o",
               facecolors=COLOR_DISK, edgecolors=SURFACE, linewidths=0.8,
               label="disk", zorder=3)

    sel = catalog[catalog["source_id"] == selected_id]
    if not sel.empty:
        ax.scatter(sel["m20"], sel["gini"], s=170, facecolors="none",
                   edgecolors=INK, linewidths=1.6, zorder=4,
                   label=f"source {selected_id}")

    xs = np.linspace(catalog["m20"].min() - 0.2, catalog["m20"].max() + 0.2, 2)
    ax.plot(xs, -0.14 * xs + 0.33, color=INK_2, lw=1, ls="--",
            label="G = −0.14·M20 + 0.33", zorder=2)

    ax.invert_xaxis()  # convention: M20 decreases to the right
    ax.set_xlabel("M20", color=INK_2, fontsize=9)
    ax.set_ylabel("Gini", color=INK_2, fontsize=9)
    ax.tick_params(colors=INK_2, labelsize=8)
    ax.grid(color=GRID, lw=0.6)
    for spine in ax.spines.values():
        spine.set_color(GRID)
    leg = ax.legend(fontsize=8, framealpha=0.9, facecolor=SURFACE,
                    edgecolor=GRID)
    for t in leg.get_texts():
        t.set_color(INK_2)
    fig.tight_layout()
    return fig


def detail_view(config: DataConfig, catalog: pd.DataFrame,
                filtered: pd.DataFrame, opts: dict):
    ids = filtered["source_id"].tolist()
    if not ids:
        st.info("No sources pass the current filters.")
        return
    if st.session_state.get("detail_select") not in ids:
        st.session_state["detail_select"] = ids[0]

    c_sel, c_over, c_panel = st.columns([2, 1, 1], vertical_alignment="bottom")
    sid = c_sel.selectbox("Source ID", ids, key="detail_select")
    overlay = c_over.checkbox("Segmentation overlay", value=True)
    show_mask = c_panel.checkbox("Mask panel", value=False)

    sci_path = os.path.join(config.cutouts_dir, f"{sid}.fits")
    if not os.path.exists(sci_path):
        st.error(f"Cutout not found: {sci_path}")
        return
    sci, seg, _header = cached_cutout(config, sid, mtime_or_zero(sci_path))

    col_img, col_metrics = st.columns([3, 2])
    with col_img:
        fig = plot_cutout(sci, seg, sid, opts, overlay, show_mask)
        st.pyplot(fig, width="stretch")
        plt.close(fig)
        if overlay and seg is not None:
            st.caption("Solid blue: this source's segmentation boundary. "
                       "Dashed orange: neighboring labels (check deblending).")
        elif seg is None:
            st.warning("No segmentation mask found for this source.")
    with col_metrics:
        metrics_panel(catalog[catalog["source_id"] == sid].iloc[0])

    if {"gini", "m20"}.issubset(catalog.columns):
        st.markdown("##### Gini–M20 plane")
        fig = gini_m20_figure(catalog, sid)
        st.pyplot(fig, width="content")
        plt.close(fig)


# --- mosaic overview --------------------------------------------------------
def jump_to_source():
    sid = st.session_state["mosaic_jump"]
    st.session_state["detail_select"] = sid
    st.session_state["_last_table_sel"] = sid


def mosaic_view(config: DataConfig, catalog: pd.DataFrame):
    if not config.mosaic_path or not os.path.exists(config.mosaic_path):
        st.info("No mosaic file configured (set the path in the sidebar).")
        return
    image, stride = cached_mosaic(config.mosaic_path,
                                  mtime_or_zero(config.mosaic_path))
    if stride > 1:
        st.caption(f"Downsampled ×{stride} for display.")

    fig, ax = plt.subplots(figsize=(9, 9))
    fig.patch.set_facecolor(SURFACE)
    finite = image[np.isfinite(image)]
    vmin, vmax = ZScaleInterval().get_limits(finite)
    ax.imshow(image, origin="lower", cmap="gray", vmin=vmin, vmax=vmax)

    if {"x", "y"}.issubset(catalog.columns):
        has_disk = "is_disk" in catalog.columns
        disks = catalog[catalog["is_disk"]] if has_disk else catalog.iloc[0:0]
        others = catalog[~catalog["is_disk"]] if has_disk else catalog
        ax.scatter(others["x"] / stride, others["y"] / stride, s=60,
                   marker="^", facecolors="none", edgecolors=COLOR_NONDISK,
                   linewidths=1.2, label="non-disk")
        ax.scatter(disks["x"] / stride, disks["y"] / stride, s=60, marker="o",
                   facecolors="none", edgecolors=COLOR_DISK, linewidths=1.2,
                   label="disk")
        sel = catalog[catalog["source_id"] == st.session_state.get("detail_select")]
        if not sel.empty:
            ax.scatter(sel["x"] / stride, sel["y"] / stride, s=220,
                       facecolors="none", edgecolors=INK, linewidths=1.6,
                       label="selected")
        leg = ax.legend(fontsize=9, framealpha=0.9, facecolor=SURFACE,
                        edgecolor=GRID)
        for t in leg.get_texts():
            t.set_color(INK_2)
    ax.tick_params(colors=INK_2, labelsize=8)
    for spine in ax.spines.values():
        spine.set_color(GRID)
    st.pyplot(fig, width="stretch")
    plt.close(fig)

    st.selectbox("Jump to source (then open the Catalog & detail tab)",
                 catalog["source_id"].tolist(), key="mosaic_jump",
                 on_change=jump_to_source)


# --- main -------------------------------------------------------------------
def main():
    st.title("JWST disk-galaxy pipeline QA")
    config = sidebar_config()
    opts = sidebar_display_options()

    if not os.path.exists(config.catalog_path):
        st.error(
            f"Catalog not found: `{config.catalog_path}`.\n\n"
            "Generate the demo data first: `python make_mock_data.py`")
        st.stop()
    catalog = cached_catalog(config.catalog_path,
                             mtime_or_zero(config.catalog_path))

    tab_catalog, tab_mosaic = st.tabs(["Catalog & detail", "Mosaic overview"])

    with tab_catalog:
        filtered = filter_catalog(catalog)
        st.caption(f"{len(filtered)} / {len(catalog)} sources pass filters. "
                   "Select a row to inspect it below.")
        event = st.dataframe(
            filtered, hide_index=True, on_select="rerun",
            selection_mode="single-row", height=280, key="catalog_table")
        if event.selection.rows:
            sel_sid = filtered.iloc[event.selection.rows[0]]["source_id"]
            # act only on selection *changes* so the selectbox still works
            if st.session_state.get("_last_table_sel") != sel_sid:
                st.session_state["_last_table_sel"] = sel_sid
                st.session_state["detail_select"] = sel_sid
        st.divider()
        detail_view(config, catalog, filtered, opts)

    with tab_mosaic:
        mosaic_view(config, catalog)


if __name__ == "__main__":
    main()
