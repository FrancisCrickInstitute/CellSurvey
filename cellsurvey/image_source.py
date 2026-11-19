"""Image sources: where stage 1 gets the input image from.

Every source is turned into an :class:`ImageSource`, a lazy ``(c, y, x)`` dask
array plus channel names, physical pixel size and provenance. Stage 1 converts
it to a SpatialData object with :meth:`ImageSource.to_spatialdata` and writes
``<prefix>.zarr``; everything after that reads only the local Zarr, so later
stages don't know or care where the image came from.

Sources:

- local multichannel TIFFs, read exactly as before (``sopa.io.ome_tif`` and
  the OME pixel size via bioio);
- ``omero://host[:port]/Image:<id>``: see :mod:`cellsurvey.omero_source`
  (omero-py is imported only for these).
"""

import os
from contextlib import nullcontext
from dataclasses import dataclass, field
from typing import Callable

import dask
import dask.array as da
import sopa
from sopa.constants import SopaAttrs
from spatialdata import SpatialData
from spatialdata.models import Image2DModel
from xarray import DataArray

from cellsurvey.segmentation import get_pixel_size

# Same layout as sopa.io.ome_tif, so downstream stages see an identical image.
IMAGE_CHUNKS = (1, 1024, 1024)
IMAGE_SCALE_FACTORS = [2, 2, 2, 2]

# Key under SpatialData.attrs (saved in the Zarr) for CellSurvey's own metadata.
ATTRS_KEY = 'cellsurvey'


@dataclass
class ImageSource:
    """A lazily read multichannel image and its metadata.

    Attributes:
        data: Lazy ``(c, y, x)`` dask array with an integer dtype.
        channel_names: One name per channel.
        name: Image element name in the SpatialData object.
        pixel_size_um: Physical pixel size in µm/px from the source's metadata,
            or ``None`` if the source has none.
        source: Provenance, recorded in the Zarr (e.g. ``{'type': 'file', 'path': ...}``).
        readers: Thread count for reading while writing the Zarr; ``None``
            keeps dask's default scheduler (used for local files).
        closer: Releases the source (e.g. OMERO pixel stores and connection);
            called by :meth:`close`. ``None`` for local files.
    """
    data: da.Array
    channel_names: list
    name: str
    pixel_size_um: float | None
    source: dict = field(default_factory=dict)
    readers: int | None = None
    closer: Callable | None = field(default=None, repr=False)

    def close(self):
        """Release the source; the lazy ``data`` can't be read afterwards. Safe to call twice."""
        closer, self.closer = self.closer, None
        if closer is not None:
            closer()

    def to_spatialdata(self, pixel_size_um):
        """Build the SpatialData object for ``<prefix>.zarr``.

        Mirrors ``sopa.io.ome_tif`` (same chunks, pyramid and Sopa attribute)
        and also records ``attrs['cellsurvey'] = {pixel_size_um, source}``.
        ``pixel_size_um`` is the effective value (``--pixel-size`` wins over
        the source's own metadata).
        """
        image = DataArray(self.data.rechunk(IMAGE_CHUNKS), dims=['c', 'y', 'x'],
                          name=self.name, coords={'c': self.channel_names})
        image = Image2DModel.parse(image, c_coords=self.channel_names,
                                   chunks=IMAGE_CHUNKS, scale_factors=IMAGE_SCALE_FACTORS)
        attrs = {
            SopaAttrs.CELL_SEGMENTATION: self.name,
            ATTRS_KEY: {'pixel_size_um': float(pixel_size_um), 'source': dict(self.source)},
        }
        return SpatialData(images={self.name: image}, attrs=attrs)

    def write(self, sdata, zarr_path):
        """Write ``sdata`` to ``zarr_path``, reading the source with ``readers`` threads."""
        ctx = (dask.config.set(scheduler='threads', num_workers=self.readers)
               if self.readers else nullcontext())
        with ctx:
            sdata.write(zarr_path, overwrite=True)


def stored_pixel_size(sdata):
    """Return the pixel size (µm/px) stored by :meth:`ImageSource.to_spatialdata`, or ``None``.

    ``None`` for Zarrs written before it was stored, or if the stored value is
    not a positive number.
    """
    value = sdata.attrs.get(ATTRS_KEY, {}).get('pixel_size_um')
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if value > 0 else None


def is_omero(path):
    return str(path).lower().startswith('omero://')


def open_image_source(path, readers=None):
    """Open the input image at ``path`` as an :class:`ImageSource`.

    ``omero://host[:port]/Image:<id>`` opens an OMERO image (``readers`` =
    parallel tile readers; see :mod:`cellsurvey.omero_source`). Anything else
    is a local multichannel TIFF, read with ``sopa.io.ome_tif`` (dask_image/
    tifffile, 3D ``(c, y, x)`` or 4D ``(1, y, x, c)``) and the OME physical
    pixel size via bioio. Call :meth:`ImageSource.close` when done.
    """
    if is_omero(path):
        from cellsurvey.omero_source import open_omero_image
        return open_omero_image(path, readers=readers)
    image = sopa.io.ome_tif(path, as_image=True)  # rechunked, integer dtype, channel coords
    return ImageSource(
        data=image.data,
        channel_names=[str(c) for c in image.coords['c'].values],
        name=image.name,
        pixel_size_um=get_pixel_size(path),
        source={'type': 'file', 'path': os.path.abspath(path)},
    )
