# CellSurvey Documentation Plan

Status: **PLAN** — nothing here is generated yet. This file captures the
intended structure, page inventory, and decisions for a MkDocs Material site
hosted on ReadTheDocs. Delete or archive this file once the real docs exist.

---

## 1. Decisions

| Decision | Choice | Rationale |
|---|---|---|
| Static site generator | **MkDocs** + **Material for MkDocs** | Lighter than Sphinx; better for prose-heavy, mostly-manual docs. No docstrings to render (code is lightly commented). |
| Hosting | **ReadTheDocs** | Versioned, auto-builds from git tags/branches, free for open source. |
| CLI reference | Hand-authored Markdown table (auto-generated later via `mkdocs` + a small script) | argparse already has all flags/help; keep it in sync manually for now, or add a tiny generator. |
| API reference | **Not initially** | 6 modules, minimal public API; documenting internals is lower value than user-facing guides. Revisit if the library surface grows. |
| Python/version pinning | None added | Docs only — no new runtime dependencies. |
| Audience | **Less computationally-literate** than Sopa | Guides written for biologists/researchers, not Python developers (see §4). |

## 2. Repo layout (target)

```
docs/
├── index.md                      # Landing / overview + pipeline diagram + quick-start
├── installation.md               # pixi + Docker
├── getting-started.md            # First run, full example, expected outputs
├── pipeline.md                   # 10-stage walkthrough (user-facing)
├── parameters.md                 # Full CLI reference with tuning guidance
├── outputs.md                    # Every output file/plot explained
├── visualization.md              # Odon / QuPath / TissUUmaps / napari
├── faq.md                        # Question-phrased how-to + troubleshooting (merged)
├── reference/                    # (future) library reference
└── DOCS_OUTLINE.md               # This file (temporary)
mkdocs.yml                        # MkDocs + Material config
.readthedocs.yaml                 # ReadTheDocs build config
```

## 3. Page-by-page content

### 3.1 `index.md`
- One-paragraph summary of what CellSurvey does.
- "Who is this for" (spatial-biology researchers, bioinformaticians).
- **Pipeline overview diagram** (Mermaid or a static image like Sopa's
  `overview_white.png`) showing the main steps. This is the first thing a
  non-technical reader needs.
- Links to Installation and Getting Started.
- Note the entry points: `run.py` (full pipeline) and `python -m cellsurvey.cli`
  (standalone module run). Crash recovery / re-analysis from an existing Zarr is
  handled by the `--resume-from` flag, not a separate script.

### 3.2 `installation.md`
- Pixi install (the only supported path).
- Docker build/run.
- The TF ≥2.18 + `tf_keras` + `TF_USE_LEGACY_KERAS=1` requirement and *why*
  (cuDNN 9 autotuner). Cross-link to `faq.md`.
- Linux-only note (no Windows/macOS via pixi).
- Pure setup only — no analysis walkthrough (mirrors Sopa's Getting Started).

### 3.3 `getting-started.md`
- Minimal command.
- Full example with realistic paths.
- What you should see on a successful run (list of generated files).
- Pointer to `outputs.md` to interpret them.

### 3.4 `pipeline.md` (the big new piece)
Stage-by-stage, each with: *what it does*, *inputs → outputs*, *where in code*.

1. Image loading (`sopa.io.ome_tif`, `BioImage` channel names).
2. Blob/spot detection (only with `--detect-blobs`; LoG via `blob_log`, tiled).
3. Initial Zarr write (checkpoint).
4. Stardist segmentation (`2D_versatile_fluo`, GPU vs CPU).
5. Channel aggregation (`sopa.aggregate`, spots optional).
6. K-means clustering (`StandardScaler` + `KMeans`).
7. Network analysis (`scipy.spatial.Delaunay` + `networkx` Louvain).
8. Spot-to-cell assignment (`gpd.sjoin`, `within`).
9. QuPath GeoJSON export.
10. Spatial neighborhood / UMAP / Leiden / density plots.

Include a simple data-flow diagram (Mermaid, which Material supports natively).

### 3.5 `parameters.md` (replaces the sparse README section)
Full table of every flag + default + type, grouped as in the README, **plus** a
"tuning guidance" subsection per group. Key relationships to document:

- `--channels` / `--thresholds` **must have equal length** (CLI errors otherwise).
- Blob detection: `--threshold` controls sensitivity; `--min/max/num-sigma`
  control the LoG scale range; `--tile-size`/`--overlap`/`--workers` are
  performance knobs (overlap prevents double-counting at tile seams).
- `--n-clusters` sets k-means k (a *required* input, no auto-selection).
- `--community-resolution` = Louvain resolution (higher → more, smaller communities).
- `--max-edge-distance` = Delaunay graph density dial (larger → denser network).
- `--radius-min/max` = spatial-neighbors graph band for the mean-hop heatmap.
- `--resume-from` = crash recovery; two-level checkpoint with `_seg.zarr` reuse.
- `--use-gpu` vs auto-detection.

**Style note:** present this as a *walkthrough* (like Sopa's `cli_usage.md`),
not a bare flag dump. Lead with defaults ("most runs need no changes"), then
explain when to reach for each flag.

### 3.6 `outputs.md`
Every artifact, what it means, and how to open it:

- `*.zarr` (initial) and `*_seg.zarr` (segmented) — SpatialData containers.
- `summary.json` — cell/cluster/community/edge counts.
- `cell_type_to_cell_type.png` — mean hop distance between clusters.
- `umap_kmeans_cluster.png` / `umap_leiden.png` — expression-space embeddings.
- `cell_density.png` — per-cluster/per-community spatial density maps.
- `cluster_intensity_heatmap.png` — mean z-scored intensity per cluster.
- `morphology_by_cluster.png` — cell area by cluster (only if `area` present).
- `qupath_export.geojson` — cells + spots for QuPath.

### 3.7 `visualization.md`
Move/expand the existing README "Visualising Results" section:
Odon (recommended), QuPath, TissUUmaps, napari + napari-spatialdata, with the
performance notes and the ≤24-category color-by caveat.

### 3.8 `faq.md` (question-phrased how-to + troubleshooting, merged)
Model on Sopa's `faq.md`: H2 headings phrased as **user questions**, each with a
plain-language answer first, then a copy-pasteable snippet if needed. Folds in the
former `troubleshooting.md` error cases as "it broke" questions.

Proposed question list (draft — expand during writing):

- What are the inputs and outputs of CellSurvey?
- Do I need a GPU?
- What do I do if "Autotuner could not find any supported configs" appears?
  → `tf_keras` + `TF_USE_LEGACY_KERAS=1` (cuDNN 9 issue).
- Why do I get `ModuleNotFoundError: No module named 'tf_keras'`?
- GPU isn't detected — what now? → `--use-gpu`, or run CPU.
- How do I resume a crashed run? → `--resume-from` + `_seg.zarr` reuse.
- Which parameters should I change from their defaults?
- Can I process formats other than OME-TIFF?
- What does "No cell was returned by the segmentation" mean? (empty/edge patches)

### 3.9 `reference/` (future, deferred)
- Auto-generated CLI reference.
- Module/API reference if the public surface grows.

## 4. Sopa reference patterns (borrow / avoid)

Sopa's docs (MkDocs Material, `prism-oncology/sopa/docs`) are the closest
structural template, but Sopa targets **developers/advanced users**. Borrow its
*structure*, soften its *content*.

### Borrow

- **Home = pitch + pipeline diagram** — a single visual of the main steps, as the
  first thing a non-technical reader sees.
- **Getting Started = pure install/setup**, no analysis walkthrough; links out to
  tutorials/CLI for the actual work.
- **FAQ phrased as user questions** ("How do I…", "What are the inputs…") — the
  de-facto how-to for non-experts, navigated by *their* questions, not module names.
- **Tabbed content** (`=== "Tab name"`) for install variants / per-technology commands.
- **Admonitions** (`!!! note/tip/warning`) — already matches the README style.
- **CLI as a progressive walkthrough** with `--help` pointers, not a bare flag table.
- **Per-data-type tips page** (Sopa's `techno_specific.md`) — collects "for X data,
  use these settings" in one place.

### Avoid (too developer/API-focused for our audience)

- **API-first Getting Started** — no `sopa.settings.*` globals, `SpatialData`/
  `AnnData` objects, dask clients, env vars, or cluster-memory math up front.
- **Auto-generated API pages** (`::: module`) surfaced prominently — defer to a
  secondary `reference/` section or skip.
- **Advanced tutorials at top level** (custom segmentation, alignment,
  parallelization internals) — overwhelm a less-technical audience.
- **Parameter-heavy commands without "just use the default" guidance**.
- **FAQ answers that jump straight into Python** (`sopa.settings.auto_save_on_disk`,
  `logging` levels) — give plain-language answers first.

## 5. MkDocs Material features to enable

- `material` theme with navigation, search, and code-copy.
- Mermaid diagrams (pipeline data flow).
- Admonitions (`!!! note`, `!!! warning`) — already used in README style.
- `tables` for the CLI reference.
- `content.code.annotate` / tabs for platform-specific commands where needed.

## 6. Build & deploy

- `mkdocs.yml` at repo root.
- `.readthedocs.yaml` with `mkdocs` build tool and `python: "3.12"`.
- `requirements-docs.txt` (or pixi `[tasks]`) pinning `mkdocs-material`.
- Local preview: `mkdocs serve`.
- ReadTheDocs auto-builds on push to `main` and on tags.

## 7. Open questions (resolve before generating)

1. Keep the CLI reference **manual** vs. **generated from argparse**?
   (Generated is more maintainable; a 10-line script can dump `--help` to Markdown.)
2. Do we want **versioned** docs (ReadTheDocs per-tag) or a single `latest`?
3. Should the README's Parameters/Visualising sections be **trimmed** and point
   to the docs site, or kept duplicated? (Recommend: keep quick-start, link out.)
4. Does `--resume-from` (including its segmented-Zarr reuse behaviour) get its own
   subsection in `pipeline.md`, or live primarily under `parameters.md` and
   `faq.md`? (Recommend: cover it in `pipeline.md` + a FAQ entry.)
5. **FAQ vs. separate troubleshooting page**: merge (as outlined above, Sopa-style)
   or keep `troubleshooting.md` as a distinct error-reference? (Recommend: merge
   for a single help page; split later only if it grows unwieldy.)

## 8. Order of work (after sign-off)

1. Scaffold `mkdocs.yml`, `.readthedocs.yaml`, `requirements-docs.txt`, `docs/`.
2. Write `index.md`, `installation.md`, `getting-started.md`.
3. Write `pipeline.md` (Mermaid diagram + stage table).
4. Write `parameters.md` (flags + tuning).
5. Write `outputs.md`, `visualization.md`.
6. Write `faq.md` (question-phrased, folding in troubleshooting).
7. Trim README to a quick-start that links to docs.
8. Verify `mkdocs build` and `mkdocs serve` locally.
9. Connect ReadTheDocs and confirm a clean build.
