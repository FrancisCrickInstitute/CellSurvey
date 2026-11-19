# Reading from OMERO

CellSurvey can read its input image straight from an
[OMERO](https://www.openmicroscopy.org/omero/) server instead of a local TIFF.
Pass the image as a URL to `-i`:

```bash
pixi run python run.py \
  -i omero://omero-prod.camp.thecrick.org/Image:602778 \
  -o results/slide1 -p results/slide1_plots
```

Everything else works as with a TIFF, including `--detect-blobs`, either
segmenter and `--resume-from`.

!!! warning "New feature"
    OMERO input is new. Please report anything unexpected.

## How it works

- The image is read **once**, at the start of the run, and copied into
  `<output>.zarr`. Every later stage, including `--resume-from`, works from that
  local copy and never contacts OMERO again.
- Pixels are read at **full resolution** in 1024 × 1024 tiles, with several
  connections in parallel (`--omero-readers`, default 8). On the Crick network
  this ran at roughly 300 MB/s, about 2 minutes for a 37 GB slide.
- Channel names and the physical pixel size come from OMERO, so `--pixel-size`
  is only needed if OMERO has no pixel size for the image.
- OMERO has already decoded the original file (with Bio-Formats), so the
  original format doesn't matter. Only integer pixel types are supported.
- Only the first Z-plane and timepoint are read.
- The server and image ID are recorded in the Zarr (`attrs['cellsurvey']['source']`).

## URL form

`omero://<host>[:<port>]/Image:<id>`. The port defaults to 4064, and
`omero://<host>/<id>` also works. You don't need to say which group the image is
in: CellSurvey searches all your groups.

## Logging in

Your password is never passed on the command line. CellSurvey tries, in order:

1. the `OMERO_SESSION_KEY` environment variable;
2. the session saved by `omero login`, if it is for the same server;
3. the `OMERO_USER` and `OMERO_PASSWORD` environment variables.

For interactive use, log in once with the OMERO command-line client (it is
installed in the CellSurvey pixi environments):

```bash
pixi run omero login -s omero-prod.camp.thecrick.org -u <username>
```

!!! tip "Batch jobs"
    Sessions expire after a period of inactivity (100 minutes on the Crick
    server), so a session created before a job waits in the queue may have
    expired by the time it starts. Either log in inside the job script, or
    create a session and pass its key:

    ```bash
    export OMERO_SESSION_KEY=$(pixi run omero sessions key)
    ```

    Reading only happens in the first few minutes of a run, so the session only
    needs to last until then.

## Which image ID?

One slide often appears as several OMERO images (the scan itself, a thumbnail,
a slide-label photo and an overview photo). Use the **full-resolution scan**.
If you give the ID of a much smaller image from the same file, the run stops
and suggests the right ID. If you give a scan that is smaller but not tiny, for
example a second scan in a multi-scan file, it warns and continues.
