import numpy as np
import dask
from dask import delayed
from skimage.feature import blob_log


def detect_blobs_in_tile(image, y, x, tile_size, overlap, image_min, image_max, min_sigma, max_sigma, num_sigma,
                         threshold):
    height, width = image.shape
    y_end = min(y + tile_size, height)
    x_end = min(x + tile_size, width)

    actual_h = y_end - y
    actual_w = x_end - x

    tile = image[y:y_end, x:x_end].astype(float)
    tile = (tile - image_min) / (image_max - image_min + 1e-8)

    if actual_h < tile_size or actual_w < tile_size:
        padded = np.zeros((tile_size, tile_size), dtype=tile.dtype)
        padded[:actual_h, :actual_w] = tile
        tile = padded

    blobs = blob_log(tile, min_sigma=min_sigma, max_sigma=max_sigma, num_sigma=num_sigma, threshold=threshold)

    if len(blobs) > 0:
        blobs = blobs[(blobs[:, 0] < actual_h) & (blobs[:, 1] < actual_w)]

        blobs[:, 0] += y
        blobs[:, 1] += x

        half_overlap = overlap // 2
        y_min_keep = y + half_overlap if y > 0 else y
        x_min_keep = x + half_overlap if x > 0 else x
        y_max_keep = y_end - half_overlap if y_end < height else y_end
        x_max_keep = x_end - half_overlap if x_end < width else x_end

        mask = (
                (blobs[:, 0] >= y_min_keep) & (blobs[:, 0] < y_max_keep) &
                (blobs[:, 1] >= x_min_keep) & (blobs[:, 1] < x_max_keep)
        )
        blobs = blobs[mask]

    return blobs


def detect_blobs_tiled(image, tile_size=1024, overlap=50, min_sigma=2, max_sigma=5, num_sigma=5, threshold=0.1,
                       n_workers=8):
    print(f"Detecting blobs in tiles of {tile_size}x{tile_size} with {overlap}px overlap using {n_workers} workers...")

    height, width = image.shape
    image_min = image.min()
    image_max = image.max()

    y_starts = list(range(0, height, tile_size - overlap))
    x_starts = list(range(0, width, tile_size - overlap))
    total_tiles = len(y_starts) * len(x_starts)
    print(f"Total tiles: {total_tiles}")

    tasks = [
        delayed(detect_blobs_in_tile)(
            image, y, x, tile_size, overlap, image_min, image_max,
            min_sigma, max_sigma, num_sigma, threshold
        )
        for y in y_starts
        for x in x_starts
    ]

    results = dask.compute(*tasks, scheduler='threads', num_workers=n_workers)

    all_blobs = [r for r in results if len(r) > 0]
    if all_blobs:
        all_blobs = np.concatenate(all_blobs, axis=0)
    else:
        all_blobs = np.zeros((0, 3))

    print(f"Found {len(all_blobs)} blobs total")
    return all_blobs
