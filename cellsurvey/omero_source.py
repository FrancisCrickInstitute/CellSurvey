"""Read an image from an OMERO server as an :class:`~cellsurvey.image_source.ImageSource`.

Used by :func:`cellsurvey.image_source.open_image_source` for
``omero://host[:port]/Image:<id>`` inputs (``omero://host/<id>`` also works).
omero-py is imported only here, so ``cellsurvey`` imports without it.

- **Login**, tried in order (never a password on the command line):
  ``OMERO_SESSION_KEY``; the session saved by ``omero login`` for the same
  host; ``OMERO_USER`` + ``OMERO_PASSWORD``.
- **Groups**: the image is looked up across all groups (client-side
  ``omero.group=-1``), and pixels are read with a call context set to the
  image's own group. ``setGroupForSession`` is deliberately not used: with a
  reused ``omero login`` session it would also switch the group of the user's
  CLI session.
- **Pixels**: full resolution, Z=0, T=0, read lazily in 1024² tiles. Each dask
  thread uses its own ``RawPixelsStore`` (``getTile`` + ``np.frombuffer``;
  omero-py's ``getTiles`` opens a store per call and decodes with ``struct``,
  too slow for whole slides).
- **Read once**: stage 1 writes the image to ``<prefix>.zarr`` and closes the
  source; nothing later contacts OMERO.
"""

import os
import re
import threading

import numpy as np
import dask.array as da

DEFAULT_PORT = 4064
DEFAULT_READERS = 8      # parallel tile readers (~300 MB/s on the Crick server)
TILE = 1024              # multiple of the server's 512 px native tiles
TILE_RETRIES = 3
# Images below this fraction of the largest image's area in their fileset are
# thumbnails / label / overview photos, not the scan.
MIN_FILESET_AREA_FRACTION = 0.01
MAX_FILESET_IMAGES_CHECKED = 100

_URL_RE = re.compile(r'^omero://(?P<host>[^/:\s]+)(?::(?P<port>\d+))?/(?:image[:/])?(?P<id>\d+)/?$',
                     re.IGNORECASE)


def parse_omero_url(url):
    """Return ``(host, port, image_id)`` from ``omero://host[:port]/Image:<id>`` (or ``/<id>``)."""
    m = _URL_RE.match(str(url).strip())
    if m is None:
        raise ValueError(f"Not a valid OMERO image URL: {url!r}. Expected "
                         f"omero://<host>[:<port>]/Image:<id>, e.g. omero://omero.example.org/Image:1234")
    return m['host'], int(m['port'] or DEFAULT_PORT), int(m['id'])


def _try_connect(conn, **kwargs):
    try:
        return bool(conn.connect(**kwargs))
    except Exception:
        return False


def connect(host, port=DEFAULT_PORT):
    """Log in to OMERO. Returns ``(conn, owns_session, how)``.

    ``owns_session`` is True only when we created the session (user/password);
    then it is ended on close. Reused sessions are left open.
    """
    from omero.gateway import BlitzGateway

    tried = []
    key = os.environ.get('OMERO_SESSION_KEY')
    if key:
        conn = BlitzGateway(host=host, port=port)
        if _try_connect(conn, sUuid=key):
            return conn, False, 'OMERO_SESSION_KEY'
        tried.append('OMERO_SESSION_KEY (rejected or expired)')

    try:
        from omero.util.sessions import SessionsStore
        s_host, s_user, s_uuid, s_port = SessionsStore().get_current()
    except Exception:
        s_host = s_user = s_uuid = s_port = None
    if s_uuid and s_host and s_host.strip().lower() == host.lower():
        conn = BlitzGateway(host=host, port=int(s_port or port))
        if _try_connect(conn, sUuid=s_uuid):
            return conn, False, f"saved `omero login` session (user {s_user})"
        tried.append(f"saved `omero login` session for {s_host} (expired?)")
    elif s_uuid:
        tried.append(f"saved `omero login` session is for {s_host}, not {host}")

    user, password = os.environ.get('OMERO_USER'), os.environ.get('OMERO_PASSWORD')
    if user and password:
        conn = BlitzGateway(user, password, host=host, port=port, secure=True)
        if _try_connect(conn):
            return conn, True, f"OMERO_USER/OMERO_PASSWORD (user {user})"
        tried.append('OMERO_USER/OMERO_PASSWORD (login failed)')

    detail = f" Tried: {'; '.join(tried)}." if tried else ""
    raise RuntimeError(
        f"Could not log in to OMERO at {host}:{port}. Run `omero login -s {host} -u <user>` first "
        f"(sessions expire after a period of inactivity), or set OMERO_SESSION_KEY, or OMERO_USER "
        f"and OMERO_PASSWORD.{detail}")


def _check_fileset(image, image_id, host):
    """Stop if ``image`` is a tiny member of its fileset (thumbnail/label/overview); warn if not the largest."""
    fileset = image.getFileset()
    if fileset is None:
        return
    members = list(fileset.copyImages())
    if len(members) <= 1 or len(members) > MAX_FILESET_IMAGES_CHECKED:
        return
    area = lambda im: (im.getSizeX() or 0) * (im.getSizeY() or 0)
    largest = max(members, key=area)
    if largest.getId() == image_id or area(largest) <= area(image):
        return
    describe = lambda im: f"Image:{im.getId()} '{im.getName()}' ({im.getSizeX()} x {im.getSizeY()} px)"
    if area(image) < MIN_FILESET_AREA_FRACTION * area(largest):
        raise ValueError(
            f"{describe(image)} is a small image in its fileset (a thumbnail, label or overview photo?). "
            f"The largest is {describe(largest)}. Did you mean omero://{host}/Image:{largest.getId()}?")
    print(f"WARNING: {describe(image)} is not the largest image in its fileset; the largest is "
          f"{describe(largest)}. Continuing with Image:{image_id}.")


class _TileReader:
    """Reads ``(y, x)`` tiles of one channel at full resolution, one RawPixelsStore per thread."""

    def __init__(self, conn, pixels_id, ctx, size_x, size_y, dtype):
        self._conn, self._pixels_id, self._ctx = conn, pixels_id, ctx
        self._size = (size_x, size_y)
        self.dtype = np.dtype(dtype)
        self._wire = self.dtype.newbyteorder('>')   # OMERO sends big-endian bytes
        self._local = threading.local()
        self._stores = []
        self._lock = threading.Lock()
        self._store()   # open one now: fails early and checks the resolution level

    def _store(self):
        rps = getattr(self._local, 'rps', None)
        if rps is not None:
            return rps
        rps = self._conn.createRawPixelsStore()
        rps.setPixelsId(self._pixels_id, True, self._ctx)
        levels = rps.getResolutionLevels()
        if levels > 1:
            # Levels are numbered smallest = 0 .. largest = n-1.
            rps.setResolutionLevel(levels - 1)
            try:
                full = rps.getResolutionDescriptions()[0]   # largest first
                if (full.sizeX, full.sizeY) != self._size:
                    raise RuntimeError(f"OMERO full-resolution level is {full.sizeX} x {full.sizeY}, "
                                       f"expected {self._size[0]} x {self._size[1]}")
            except AttributeError:
                pass   # older servers without getResolutionDescriptions
        self._local.rps = rps
        with self._lock:
            self._stores.append(rps)
        return rps

    def _drop_store(self):
        rps = getattr(self._local, 'rps', None)
        self._local.rps = None
        if rps is not None:
            try:
                rps.close()
            except Exception:
                pass

    def read(self, c, y0, y1, x0, x1):
        h, w = y1 - y0, x1 - x0
        error = None
        for _ in range(TILE_RETRIES):
            try:
                raw = self._store().getTile(0, c, 0, x0, y0, w, h)
                return np.frombuffer(raw, dtype=self._wire).reshape(h, w).astype(self.dtype)
            except Exception as e:   # transient Ice/server errors: retry with a fresh store
                error = e
                self._drop_store()
        raise RuntimeError(f"Failed to read OMERO tile c={c}, y={y0}:{y1}, x={x0}:{x1} "
                           f"after {TILE_RETRIES} attempts") from error

    def close(self):
        with self._lock:
            stores, self._stores = self._stores, []
        for rps in stores:
            try:
                rps.close()
            except Exception:
                pass


class _LazyTiles:
    """Array-like ``(c, y, x)`` view over a tile-reading function, for ``dask.array.from_array``."""

    def __init__(self, read_tile, shape, dtype):
        self._read_tile = read_tile
        self.shape = tuple(shape)
        self.dtype = np.dtype(dtype)
        self.ndim = 3

    def __getitem__(self, key):
        key = key if isinstance(key, tuple) else (key,)
        key = key + (slice(None),) * (3 - len(key))
        bounds = []
        for k, n in zip(key, self.shape):
            if not isinstance(k, slice):
                k = slice(k, k + 1)
            start, stop, step = k.indices(n)
            if step != 1:
                raise IndexError("strided reads are not supported")
            bounds.append((start, max(start, stop)))
        (c0, c1), (y0, y1), (x0, x1) = bounds
        out = np.empty((c1 - c0, y1 - y0, x1 - x0), dtype=self.dtype)
        if out.size:
            for i, c in enumerate(range(c0, c1)):
                out[i] = self._read_tile(c, y0, y1, x0, x1)
        return out


def lazy_tile_array(read_tile, shape, dtype, name):
    """A dask ``(c, y, x)`` array whose 1 x 1024 x 1024 chunks call ``read_tile(c, y0, y1, x0, x1)``."""
    return da.from_array(_LazyTiles(read_tile, shape, dtype), chunks=(1, TILE, TILE), name=name,
                         meta=np.empty((0, 0, 0), dtype=dtype), lock=False, asarray=False, fancy=False)


def _pixel_size_um(image):
    """In-plane pixel size in µm/px, or None; warns for non-square pixels."""
    def value(length):
        if length is None:
            return None
        v = length.getValue() if hasattr(length, 'getValue') else float(length)
        return float(v) if v and v > 0 else None
    x = value(image.getPixelSizeX(units="MICROMETER"))
    y = value(image.getPixelSizeY(units="MICROMETER"))
    if x is None:
        return y
    if y is not None and abs(y - x) > 0.01 * max(x, y):
        print(f"WARNING: anisotropic pixels (Y={y:.4f}, X={x:.4f} µm/px); using X={x:.4f}")
    return x


def open_omero_image(url, readers=None):
    """Open ``omero://host[:port]/Image:<id>`` as an :class:`~cellsurvey.image_source.ImageSource`.

    Nothing is read until the returned array is computed (stage 1's Zarr
    write). Call :meth:`ImageSource.close` afterwards to close the pixel stores
    and the connection.
    """
    from omero.gateway import OMERO_NUMPY_TYPES
    from cellsurvey.image_source import ImageSource

    host, port, image_id = parse_omero_url(url)
    conn, owns_session, how = connect(host, port)
    reader = None
    try:
        conn.c.enableKeepAlive(60)
        print(f"Connected to OMERO {host}:{port} via {how}")

        # Look the image up across all groups (client-side only; see module docstring).
        conn.SERVICE_OPTS.setOmeroGroup('-1')
        image = conn.getObject('Image', image_id)
        if image is None:
            raise ValueError(f"Image:{image_id} not found on {host}, or you don't have access to it")
        group = image.getDetails().getGroup()
        ctx = conn.SERVICE_OPTS.copy()
        ctx.setOmeroGroup(str(group.getId()))

        _check_fileset(image, image_id, host)

        size_x, size_y, size_c = image.getSizeX(), image.getSizeY(), image.getSizeC()
        size_z, size_t = image.getSizeZ(), image.getSizeT()
        if size_z > 1 or size_t > 1:
            print(f"WARNING: Image:{image_id} has {size_z} Z-planes and {size_t} timepoints; "
                  f"using Z=0, T=0 only")
        pixels = image.getPrimaryPixels()
        pixels_type = pixels.getPixelsType().value
        dtype = np.dtype(OMERO_NUMPY_TYPES[pixels_type])
        if not np.issubdtype(dtype, np.integer):
            raise ValueError(f"Image:{image_id} has {pixels_type} pixels; only integer pixel types "
                             f"are supported for now")
        labels = image.getChannelLabels() or []
        channel_names = [str(labels[i]) if i < len(labels) and labels[i] else str(i) for i in range(size_c)]
        pixel_size_um = _pixel_size_um(image)

        print(f"OMERO Image:{image_id} '{image.getName()}' (group {group.getName()}): "
              f"{size_x} x {size_y} px, {size_c} channels, {pixels_type}, "
              f"{size_x * size_y * size_c * dtype.itemsize / 1e9:.1f} GB")
        print(f"Channels: {channel_names}")

        reader = _TileReader(conn, pixels.getId(), ctx, size_x, size_y, dtype)
        data = lazy_tile_array(reader.read, (size_c, size_y, size_x), dtype,
                               name=f"omero-{host}-{image_id}")
    except BaseException:
        if reader is not None:
            reader.close()
        conn.close(hard=owns_session)
        raise

    def close():
        reader.close()
        conn.close(hard=owns_session)

    return ImageSource(
        data=data,
        channel_names=channel_names,
        name=f"omero_image_{image_id}",
        pixel_size_um=pixel_size_um,
        source={'type': 'omero', 'host': host, 'port': port, 'image_id': image_id,
                'image_name': image.getName(), 'group_id': int(group.getId())},
        readers=readers or DEFAULT_READERS,
        closer=close,
    )
