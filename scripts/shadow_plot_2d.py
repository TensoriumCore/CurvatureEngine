#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np


SCALAR_FIELDS = {
    "linear_r": "Linear red",
    "linear_g": "Linear green",
    "linear_b": "Linear blue",
    "intensity": "Observed intensity",
    "optical_depth": "Optical depth",
    "effective_temperature": "Effective temperature",
    "luminance": "Linear luminance",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot a 2D map from CurvatureEngine shadow CSV exports."
    )
    parser.add_argument(
        "--input",
        default="Output/accretion_disk_data.csv",
        help="Path to the CSV exported by Shadow mode.",
    )
    parser.add_argument(
        "--field",
        default="intensity",
        choices=["rgb", *SCALAR_FIELDS.keys()],
        help="Field to plot. Use 'rgb' for a color image built from linear_r/g/b.",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output PNG path. Defaults to Output/shadow_<field>.png.",
    )
    parser.add_argument(
        "--cmap",
        default="inferno",
        help="Matplotlib colormap for scalar fields.",
    )
    parser.add_argument(
        "--log",
        action="store_true",
        help="Use logarithmic scaling for scalar fields.",
    )
    parser.add_argument(
        "--clip-percentile",
        type=float,
        default=99.5,
        help="Upper percentile used for robust normalization.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=180,
        help="Output DPI.",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Custom figure title.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display the figure interactively in addition to saving it.",
    )
    parser.add_argument(
        "--pixel-space",
        action="store_true",
        help="Plot axes in pixel coordinates instead of screen coordinates.",
    )
    return parser.parse_args()


def load_shadow_csv(csv_path: Path) -> np.ndarray:
    if not csv_path.exists():
        raise FileNotFoundError(f"Missing input CSV: {csv_path}")
    data = np.genfromtxt(csv_path, delimiter=",", names=True, dtype=np.float64)
    if data.size == 0:
        raise ValueError(f"Empty CSV: {csv_path}")
    return np.atleast_1d(data)


def infer_grid_shape(data: np.ndarray) -> tuple[int, int]:
    width = int(np.max(data["pixel_x"])) + 1
    height = int(np.max(data["pixel_y"])) + 1
    return width, height


def build_scalar_grid(data: np.ndarray, field: str) -> np.ndarray:
    width, height = infer_grid_shape(data)
    grid = np.full((height, width), np.nan, dtype=np.float64)
    x = data["pixel_x"].astype(np.int64)
    y = data["pixel_y"].astype(np.int64)

    if field == "luminance":
        values = (
            0.2126 * data["linear_r"]
            + 0.7152 * data["linear_g"]
            + 0.0722 * data["linear_b"]
        )
    else:
        values = data[field]

    grid[y, x] = values
    return grid


def build_rgb_grid(data: np.ndarray) -> np.ndarray:
    width, height = infer_grid_shape(data)
    grid = np.zeros((height, width, 3), dtype=np.float64)
    x = data["pixel_x"].astype(np.int64)
    y = data["pixel_y"].astype(np.int64)
    grid[y, x, 0] = data["linear_r"]
    grid[y, x, 1] = data["linear_g"]
    grid[y, x, 2] = data["linear_b"]
    return grid


def compute_screen_extent(data: np.ndarray) -> list[float]:
    xs = np.unique(data["screen_x"])
    ys = np.unique(data["screen_y"])
    dx = float(np.median(np.diff(xs))) if xs.size > 1 else 1.0
    dy = float(np.median(np.diff(ys))) if ys.size > 1 else 1.0
    return [
        float(xs.min() - 0.5 * dx),
        float(xs.max() + 0.5 * dx),
        float(ys.min() - 0.5 * dy),
        float(ys.max() + 0.5 * dy),
    ]


def compute_pixel_extent(data: np.ndarray) -> list[float]:
    width, height = infer_grid_shape(data)
    return [-0.5, width - 0.5, -0.5, height - 0.5]


def robust_upper(values: np.ndarray, percentile: float) -> float:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return 1.0
    upper = float(np.percentile(finite, percentile))
    if not np.isfinite(upper) or upper <= 0.0:
        upper = float(np.max(finite)) if finite.size else 1.0
    return upper if upper > 0.0 else 1.0


def normalize_rgb(rgb: np.ndarray, percentile: float) -> np.ndarray:
    upper = robust_upper(rgb, percentile)
    rgb = np.clip(rgb / upper, 0.0, 1.0)
    return np.power(rgb, 1.0 / 2.2)


def scalar_norm(values: np.ndarray, log_scale: bool, percentile: float):
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return None

    upper = robust_upper(finite, percentile)
    if log_scale:
        positive = finite[finite > 0.0]
        if positive.size == 0:
            raise ValueError("Log scale requested but field has no positive values.")
        lower = max(float(np.min(positive)), upper * 1.0e-8)
        return mcolors.LogNorm(vmin=lower, vmax=upper)

    lower = float(np.nanmin(finite))
    return mcolors.Normalize(vmin=lower, vmax=upper)


def make_title(field: str, custom_title: str | None) -> str:
    if custom_title:
        return custom_title
    if field == "rgb":
        return "Accretion Disk RGB"
    return SCALAR_FIELDS[field]


def plot_rgb(data: np.ndarray, args: argparse.Namespace, output_path: Path) -> None:
    rgb = build_rgb_grid(data)
    rgb = np.flipud(normalize_rgb(rgb, args.clip_percentile))
    extent = (
        compute_pixel_extent(data)
        if args.pixel_space
        else compute_screen_extent(data)
    )

    fig, ax = plt.subplots(figsize=(8.5, 5.0), constrained_layout=True)
    ax.imshow(rgb, origin="lower", extent=extent, interpolation="nearest")
    ax.set_title(make_title("rgb", args.title))
    ax.set_xlabel("pixel_x" if args.pixel_space else "screen_x")
    ax.set_ylabel("pixel_y" if args.pixel_space else "screen_y")
    ax.set_aspect("equal")
    fig.savefig(output_path, dpi=args.dpi)
    if args.show:
        plt.show()
    plt.close(fig)


def plot_scalar(data: np.ndarray, args: argparse.Namespace, output_path: Path) -> None:
    field_grid = build_scalar_grid(data, args.field)
    display_grid = np.flipud(field_grid)
    extent = (
        compute_pixel_extent(data)
        if args.pixel_space
        else compute_screen_extent(data)
    )
    norm = scalar_norm(display_grid, args.log, args.clip_percentile)

    fig, ax = plt.subplots(figsize=(8.5, 5.0), constrained_layout=True)
    im = ax.imshow(
        display_grid,
        origin="lower",
        extent=extent,
        cmap=args.cmap,
        norm=norm,
        interpolation="nearest",
    )
    ax.set_title(make_title(args.field, args.title))
    ax.set_xlabel("pixel_x" if args.pixel_space else "screen_x")
    ax.set_ylabel("pixel_y" if args.pixel_space else "screen_y")
    ax.set_aspect("equal")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(SCALAR_FIELDS[args.field])
    fig.savefig(output_path, dpi=args.dpi)
    if args.show:
        plt.show()
    plt.close(fig)


def main() -> int:
    args = parse_args()
    input_path = Path(args.input)
    output_path = (
        Path(args.output)
        if args.output
        else Path("Output") / f"shadow_{args.field}.png"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)

    data = load_shadow_csv(input_path)

    if args.field == "rgb":
        plot_rgb(data, args, output_path)
    else:
        plot_scalar(data, args, output_path)

    print(f"Saved 2D plot to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
