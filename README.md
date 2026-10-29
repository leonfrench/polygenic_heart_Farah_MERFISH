# Polygenic gene-set tester for the Farah et al. MERFISH developing human heart

A Streamlit app for the MERFISH developing human heart atlas of Farah et al. (2024): three frontal sections of one
13 p.c.w. (post-conception week) heart, 228,635 cells, 238-gene panel. Enter a gene set; every cell gets an exact
set-versus-complement AUROC from its gene expression ranks, shown on the three sections side by side. The app also tests
whether the gene set is specifically expressed in each **cell type** and each **cellular community**, using pseudobulk
profiles as in the [SEA-AD polygenic tester](https://github.com/leonfrench/SEA_AD_multi_polygenic_app). It is a
cell-based relative of the [spatial fetal app](https://github.com/leonfrench/spatial_fetal_app).

## Run locally

```bash
~/miniconda3/bin/conda create -y -p ./.conda-env -c conda-forge --override-channels python=3.12 pip
./.conda-env/bin/python -m pip install -r requirements-prep.txt pytest
./.conda-env/bin/streamlit run app.py
```

Open http://localhost:8501. `requirements.txt` is enough to run the app on the prepared data in
`data/merfish-heart-13pcw/`; `requirements-prep.txt` adds h5py for preprocessing.
`HEART_APP_DATA=/path/to/export` selects a different prepared dataset.

## Deploy on Streamlit Community Cloud

- Main file: `app.py`, with Python 3.12, the development version. The app and its tests pass with a clean install of
  `requirements.txt` on Python 3.12 (numpy 2.5.4, pandas 3.1.0, pyarrow 26.0.0, scipy 1.18.1 and streamlit 1.66.1; the
  preprocessing tests are skipped there). `requirements.txt` is the dependency file Community Cloud installs; do not add
  an `environment.yml`, which it would install instead.
- The prepared data (42 MB in 22 files; the largest is 16 MB) is committed in plain git, so no Git LFS is needed.
- Memory: a local server held about 0.1 GB at rest and peaked at 0.4 GB while a browser session scored every cell,
  switched through the colour modes and cell type levels, and scored a second gene set.
- The first page load sends about 3 MB of plot data to the browser.

## Using the app

- **Sections**: the three replicate experiments (`R77_4C4`, `R78_4C12`, `R78_4C15`) in one row. Each section is centred
  in its panel at a shared scale, oriented as in the paper's Fig. 1d: atria at the top, the heart's right side on the left.
- **Genes**: human symbols by default, matched case-insensitively. Mouse or rhesus macaque symbols are converted to human
  orthologs with the bundled HomoloGene snapshot (`data/homologene_to_human.csv.gz`; case-sensitive source symbols,
  one-to-many orthologs are all kept and unmapped symbols are reported). An optional background is the complete universe,
  including the targets; ranks are recomputed inside it. **Gene matching** lists the targets used and any symbols left
  out, and its download gives every gene of the scoring universe with its role (target or background). The default gene
  set is the His–Purkinje (ventricular conduction system) cardiomyocyte markers named in the paper
  (`IRX1 IRX2 IRX3 TBX3 HCN4`); it ranks `vCM-His-Purkinje` and the `VCS` community first, and that community is in all
  three sections.
- **Colour cells by**: gene-set AUROC, cell type, cellular community, cell type AUROC or cellular community AUROC (the
  last two colour each cell by its label's specificity AUROC from the tables). The cell type level (class / population /
  Leiden cluster) sets both the colours and the cell type table; communities have one level.
- **Plot**: one panel per section on a black background. Dots are sized to the median nearest-cell spacing (~9 µm).
  AUROC ranges differ between sections, so gene-set AUROC colours are relative within each section: each is scaled
  between its own 1st and 99th percentiles, which are shown above its panel, and the legend is in arbitrary units. Hover
  gives the exact AUROC. The label AUROC modes (one pooled value per label) share one scale clipped at the 1st and 99th
  percentiles of the cells shown. The clip is `COLOUR_CLIP` in `app.py`. Colours only; the statistics are unchanged.
  Drag to pan, scroll to zoom, double-click to reset, click a legend entry to highlight that label, and save the current
  view as PNG. Zooming out stops at the fitted view; from there the wheel scrolls the page (a pinch or ctrl+wheel never
  zooms the page). Hover shows a cell's gene-set AUROC, then its cell population, Leiden cluster, and cellular community
  with the section; the class is added when that level is chosen in the sidebar. Abbreviations are spelled out from the
  paper's text (`VEC`, valve endocardial cell; `VCS`, ventricular conduction system).
- **Tables**: cell type and community expression specificity across all sections (AUROC, p-value and the
  Benjamini–Hochberg adjusted pFDR), with CSV downloads, plus a download of every cell's score once scoring completes.
  The tables use pseudobulk profiles, not cell scores, so they appear as soon as the gene set is valid. Tick a row's
  checkbox to colour the plot by that table's specificity AUROC (cell type or community) and highlight the label; untick
  it to return to the previous colour mode. The population table shows the cell class as its second column; the Leiden
  cluster and community tables show each label's three most common populations.

### Annotations

- **Cell types**: the source's 27 cell populations (`populations`, the default level), its 32 Leiden clusters, and 11
  cell classes that group the populations as in Fig. 1c,d (`vCM`, `aCM`, `ncCM`, `Fibro`, `Endocardial`, `BEC`, `LEC`,
  `SMC`, `Epi`, `Neuronal`, `WBC`; pericytes are shown with the smooth muscle cells). The Leiden clusters (resolution 2)
  are not nested in the populations: 10 of the 32 span several, because the authors subclustered the ventricular
  cardiomyocytes. Cluster numbers are zero-padded to two digits so they sort in order (21 and 27 are absent from the
  source).
- **Cellular communities**: the source's 13 communities (`communities`): k-means clusters of the population mix within
  150 µm of each cell. `R78_4C15` has no AVN/AV ring or IFT/SAN cells (and only 16 ncCM cells), and `R77_4C4` has 24
  outflow tract cells.
- **Colours**: the source's scanpy colours for populations and communities. Pericytes are black in the source, the colour
  of the plot background, so they take `#D0D800`: the colour farthest (in OKLab) from the other 26 populations' colours,
  within the palette's lightness (0.6–0.85) and chroma (≤ 0.2). Classes take the colour of their most numerous
  population. Leiden clusters take the hue of their most common population, at evenly spaced OKLCH lightness where
  several clusters share one. Preprocessing checks that every colour has at least 3:1 contrast with the black plot.

## Data

The input is one file of the [Dryad deposit](https://doi.org/10.5061/dryad.w0vt4b8vp), in `data/` (not version
controlled):

| File | Contents |
| --- | --- |
| `overall_merfish.h5ad` | 228,635 cells of the three 13 p.c.w. experiments × 238 genes, log-normalized (dense `float32`); `obs`: `sample_id`, `n_counts`, `leiden`, `zone_cluster`, `communities`, `complexity`, `populations`, `purity`; `obsm/spatial`; colours in `uns` |

What the inspection showed:

- `X` is `log1p(counts / n_counts × 365)`, scanpy's `normalize_total` (to the median count) then `log1p`: `expm1(X)` sums
  to 365.000 in every cell. So `counts = expm1(X) / sum(expm1(X)) × n_counts` recovers the integer counts: every value
  rounds within 0.0002, each cell sums to its `n_counts`, and every cell's smallest detected count is 1. Preprocessing
  checks the rounding and the sums for every cell. Within-cell ranks are the same for normalized values and counts, so
  only the pseudobulks depend on the recovery.
- The paper reports 258,237 cells after segmentation and filtering; the h5ad has 228,635, all with a population and a
  community, and all are kept. Median 386 transcripts and 87 detected genes per cell; 63.4% of cell × gene entries are
  zero. The sequential smFISH genes are not in the h5ad.
- The universe is the 238 panel genes under their panel symbols; all but `FAM213A` (now `PRXL2A`) are current HGNC
  symbols. Eight panel genes cannot be reached from mouse symbols: the HomoloGene snapshot lists `ADGRL1`–`3` under their
  former symbols (`LPHN1`–`3`) and has no mouse entry for `CD24`, `FNDC1`, `PDLIM3`, `PLN` or `PPP1R12B` (`FNDC1` and
  `PDLIM3` can be reached from rhesus).
- `obsm/spatial` is in micrometres with y up, and the sections are set side by side along x (−10.7 to −5.7, −4.6 to 0.3
  and 1.4 to 6.6 mm). The atria have the larger y and the right atrium the smaller x, as in Fig. 1d.
- `zone_cluster` numbers the communities one-to-one, so only the names are kept, and `batch` repeats `sample_id`. Cell
  names are `<cell id>-<sample>`, with the segmentation's cell id as in the per-experiment Dryad files.

## Preprocessing

```bash
./.conda-env/bin/python scripts/prepare_heart_cells.py --force
```

`--raw` (default `data`) holds `overall_merfish.h5ad`. It takes about 15 seconds, with peak memory of 2.8 GB. `--force`
replaces an existing output directory once the new one is complete.

Steps:

1. **Counts**: each section's rows of the dense matrix are read from the h5ad and converted back to integer counts (see
   above).
2. **Cells**: every cell, in h5ad order (grouped by section). Display coordinates centre each section on its bounding
   box and negate y, so the atria are at the top.
3. **Ranks**: average ranks within each cell, with all undetected genes tied at the bottom, stored as doubled `uint16`
   values (format `sparse-doubled-average-ranks-v1`). Thirty-two random cells per section are checked against
   `scipy.stats.rankdata` and against a pairwise AUROC, and every batch is checked after a write/read round trip.
4. **Batches**: each section's cells are shuffled (seed 20261029) and split into equal batches of at most 16,384 cells
   (`--batch-cells`) and 12 million stored cell × gene entries (`--batch-entries`), in `progressive/<SECTION>/NNNNN.npz`.
5. **Annotations**: `annotations.json` holds labels, colours and names for each level (see Annotations above).

Output (`data/merfish-heart-13pcw/`, 42 MB, of which 0.1 MB is pseudobulk):

| File | Contents |
| --- | --- |
| `manifest.json` | Universe size, per-section cell ranges, batch paths, input checksum, count recovery, validation |
| `cells.parquet` | `cell_index`, `slice` (sample), `cell` (source name), `cell_id`, `rx, ry` (µm, display frame), `x, y` (µm, `obsm/spatial`), `cell_class`, `cell_population`, `leiden`, `community`, `complexity`, `purity` (from the source), `draw_order`, `batch`, `batch_row`, `n_counts`, `n_genes_detected` |
| `genes.parquet` | `gene_index`, panel symbol, detection and total counts |
| `gene_detection.parquet` | Detected-cell counts per gene and section (for the sparsity note) |
| `pseudobulk/<SECTION>.npz` | For each annotation level: summed panel-gene counts per (label, gene) for the labels present (`<level>_counts`, `uint32`), their global label indices (`<level>_index`, into `manifest.json` `pseudobulk.labels`) and cell counts (`<level>_cells`) |
| `annotations.json` | Cell class, population, Leiden cluster and community labels and colours, names for the abbreviations, and the populations of each Leiden cluster and community |

| Section | Cells | Size | Batches |
| --- | ---: | ---: | ---: |
| R77_4C4 | 72,962 | 5.0 × 5.4 mm | 5 |
| R78_4C12 | 75,782 | 5.0 × 6.0 mm | 5 |
| R78_4C15 | 79,891 | 5.2 × 6.3 mm | 5 |
| **Total** | **228,635** | | **15** |

## Methods

**Cell scores (plot).** For k target genes in a universe of G genes, a cell's score is
`AUROC = (sum of target ranks − k(k+1)/2) / (k(G − k))`, the probability that a target gene outranks a
non-target gene in that cell, with half credit for ties.

**Expression specificity (tables).** For the chosen level, each label's pseudobulk sums the transcript counts of its
cells across all sections and is converted to log₂(CPM + 1). Each gene is z-scored across the labels, removing its
average expression (genes that do not vary get z = 0). Within each label, genes are ranked by z-score (average ranks for
ties), and the AUROC compares the target genes with the rest of the universe: above 0.5 means the gene set is more
specific to that label than other genes. With a custom background, genes are ranked within that universe. Labels with
fewer than 20 cells are not profiled. Two-sided Mann–Whitney p-values over genes (normal approximation with tie and
continuity correction, matching `scipy.stats.mannwhitneyu`) are adjusted with Benjamini–Hochberg within each table
(pFDR).
Genes are not independent (co-expression), so the p-values are exploratory. This follows the
[SEA-AD polygenic tester](https://github.com/leonfrench/SEA_AD_multi_polygenic_app), which ranks genes by pseudobulk
z-scores in the same way; its p-values omit the tie correction.

## Tests

```bash
./.conda-env/bin/python -m pytest tests -q
node tests/test_cells.cjs
```

These cover rank encoding and scoring against SciPy and pairwise definitions, count recovery, section centring,
pseudobulk sums, annotation colours, specificity z-scores and AUROCs against `mannwhitneyu`, ortholog mapping, the grid
layouts, the score transport protocol, the canvas helpers, and end-to-end app runs on `data/merfish-heart-13pcw`.

## Data credit

Farah et al. *Spatially organized cellular communities form the developing human heart.*
Nature 627, 854–864 (2024). [doi:10.1038/s41586-024-07171-z](https://doi.org/10.1038/s41586-024-07171-z).
Data, cell populations and cellular communities: [Dryad](https://doi.org/10.5061/dryad.w0vt4b8vp) (CC0); also in the
[UCSC Cell Browser](https://cells.ucsc.edu/?ds=hoc). Orthologs: HomoloGene.
