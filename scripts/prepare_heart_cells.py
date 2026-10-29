"""MERFISH developing human heart (Farah et al. 2024) -> exact sparse ranks in per-section batches.

Input in --raw (data/, as downloaded from the Dryad deposit):
  overall_merfish.h5ad   228,635 cells of three frontal sections of one 13 p.c.w. heart x 238 panel genes, stored as
                         log1p of each cell's counts scaled to 365 (scanpy normalize_total and log1p); obs has the
                         section (sample_id), transcript count (n_counts), Leiden cluster, cell population and cellular
                         community, obsm/spatial the cell centres (micrometres, y up, sections offset side by side)

Integer counts are recovered exactly from the normalized values and n_counts. The universe is the 238 panel
genes under their panel (HGNC) symbols, as given.

  .conda-env/bin/python scripts/prepare_heart_cells.py --force
"""
from pathlib import Path
import argparse, datetime, hashlib, json, platform, shutil, time
import h5py
import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[1]
DATASET = 'Farah-heart-13pcw'
EXPRESSION = 'overall_merfish.h5ad'
# The cell classes of the 27 MERFISH populations in Fig. 1c,d (pericytes are shown with the smooth muscle cells).
CLASSES = {**dict.fromkeys(['vCM-LV-Compact', 'vCM-RV-Compact', 'vCM-LV-Trabecular', 'vCM-RV-Trabecular', 'vCM-Proliferating',
                            'vCM-His-Purkinje', 'vCM-LV-AV', 'vCM-RV-AV'], 'vCM'),
           'aCM-LA': 'aCM', 'aCM-RA': 'aCM', 'ncCM-IFT-like': 'ncCM', 'ncCM-AVC-like': 'ncCM',
           **dict.fromkeys(['aFibro', 'vFibro', 'VIC', 'adFibro'], 'Fibro'),
           **dict.fromkeys(['aEndocardial', 'vEndocardial', 'VEC'], 'Endocardial'),
           'BEC': 'BEC', 'LEC': 'LEC', 'VSMC': 'SMC', 'Pericyte': 'SMC', 'Epicardial': 'Epi', 'EPDC': 'Epi',
           'Neuronal': 'Neuronal', 'WBC': 'WBC'}
# Names for the abbreviations, from the paper's text (shown as 'name (label)' in the tooltip).
NAMES = dict(
    cell_class={'vCM': 'Ventricular cardiomyocyte', 'aCM': 'Atrial cardiomyocyte', 'ncCM': 'Non-chamber cardiomyocyte',
                'Fibro': 'Fibroblast-like', 'BEC': 'Blood endothelial cell', 'LEC': 'Lymphatic endothelial cell',
                'SMC': 'Smooth muscle cell or pericyte', 'Epi': 'Epicardial cell or EPDC', 'WBC': 'White blood cell'},
    cell_population={'vCM-LV-Compact': 'Left ventricular compact cardiomyocyte', 'vCM-RV-Compact': 'Right ventricular compact cardiomyocyte',
                     'vCM-LV-Trabecular': 'Left ventricular trabecular cardiomyocyte',
                     'vCM-RV-Trabecular': 'Right ventricular trabecular cardiomyocyte',
                     'vCM-Proliferating': 'Proliferating ventricular cardiomyocyte', 'vCM-His-Purkinje': 'His–Purkinje cardiomyocyte',
                     'vCM-LV-AV': 'Mitral valve leaflet cardiomyocyte', 'vCM-RV-AV': 'Tricuspid valve leaflet cardiomyocyte',
                     'aCM-LA': 'Left atrial cardiomyocyte', 'aCM-RA': 'Right atrial cardiomyocyte',
                     'ncCM-IFT-like': 'Inflow tract / pacemaker cardiomyocyte', 'ncCM-AVC-like': 'Atrioventricular canal / node cardiomyocyte',
                     'aFibro': 'Atrial fibroblast', 'vFibro': 'Ventricular fibroblast', 'VIC': 'Valvular interstitial cell',
                     'adFibro': 'Adventitial fibroblast', 'aEndocardial': 'Atrial endocardial cell', 'vEndocardial': 'Ventricular endocardial cell',
                     'VEC': 'Valve endocardial cell', 'BEC': 'Blood endothelial cell', 'LEC': 'Lymphatic endothelial cell',
                     'VSMC': 'Vascular smooth muscle cell', 'EPDC': 'Epicardium-derived progenitor cell', 'WBC': 'White blood cell'},
    community={'AVN/AV Ring': 'Atrioventricular node / AV ring', 'VCS': 'Ventricular conduction system',
               'IFT/SAN': 'Inflow tract / sinoatrial node', 'Inner-LV': 'Inner left ventricular wall', 'Inner-RV': 'Inner right ventricular wall',
               'Outer-LV': 'Outer left ventricular wall', 'Outer-RV': 'Outer right ventricular wall', 'OFT': 'Outflow tract',
               'Mus. Valve Leaf.': 'Muscular valve leaflet (mitral valve)', 'Subepicardial': 'Subepicardial (outer atrioventricular canal)',
               'Valve': 'Atrioventricular and outflow tract valves'})
# Pericytes are black in the source, the colour of the plot background. They take the colour farthest
# (OKLab) from the other 26 populations' colours, within their lightness (0.6-0.85) and chroma (<= 0.2).
PERICYTE = '#D0D800'
MIN_CONTRAST = 3.0  # against the black plot, as WCAG asks of graphical objects
DARK_BAND = (0.50, 0.665)  # OKLCH lightness of Leiden cluster steps: clear of a black surface, without glare
OBS = {'populations': 'cell_population', 'leiden': 'leiden', 'communities': 'community'}
LEVELS = {'cell_class': 'Cell class', 'cell_population': 'Cell population', 'leiden': 'Leiden cluster', 'community': 'Cellular community'}


def sha256(path):
    with open(path, 'rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def rank_rows(rows, cols, values, n_rows, n_genes):
    """Doubled average ranks per row; positive values only, zeros share (Z+1)/2."""
    rows, cols, values = (np.asarray(a, dtype=np.int64) for a in (rows, cols, values))
    # One integer key per sort; np.lexsort is much slower on unordered rows.
    order = np.argsort(rows*(int(values.max(initial=0))+1) + values)
    rows, cols, values = rows[order], cols[order], values[order]
    counts = np.bincount(rows, minlength=n_rows)
    zeros = n_genes - counts
    start = np.r_[0, np.flatnonzero((np.diff(rows) != 0) | (np.diff(values) != 0)) + 1]
    end = np.r_[start[1:], len(rows)]
    row_start = np.r_[0, np.cumsum(counts)[:-1]]
    # A 0-based [start,end) tie inside the positive block has doubled 1-based rank start+end+1.
    rank2 = np.repeat(start + end + 1 - 2*row_start[rows[start]] + 2*zeros[rows[start]], end - start)
    order = np.argsort(rows*n_genes + cols)
    indptr = np.r_[0, np.cumsum(counts)].astype(np.int64)
    return indptr, cols[order].astype(np.uint16), rank2[order].astype(np.uint16), np.where(zeros > 0, zeros+1, 0).astype(np.uint16)


def pseudobulk(rows, cols, values, codes, n_genes):
    """Summed counts per (label, gene) for the labels present; rows index cells, codes give each cell's label."""
    present, inverse = np.unique(codes, return_inverse=True)
    sums = np.bincount(inverse[rows]*n_genes + cols, weights=values, minlength=len(present)*n_genes)
    if sums.max(initial=0) >= 2**32:
        raise ValueError('Pseudobulk counts exceed uint32')
    return present, np.bincount(inverse, minlength=len(present)), sums.reshape(len(present), n_genes).astype(np.uint32)


def validate_rows(indptr, cols, rank2, zero2, dense_rows, n_genes, rng):
    checked = 0
    for row, dense in dense_rows:
        a, b = indptr[row], indptr[row+1]
        decoded = np.full(n_genes, zero2[row]/2); decoded[cols[a:b]] = rank2[a:b]/2
        np.testing.assert_array_equal(decoded, rankdata(dense, method='average'))
        selected = rng.choice(n_genes, 100, replace=False)
        other = np.setdiff1d(np.arange(n_genes), selected)
        truth = np.mean((dense[selected, None] > dense[other]) + 0.5*(dense[selected, None] == dense[other]))
        auc = (decoded[selected].sum() - 100*101/2) / (100*(n_genes-100))
        np.testing.assert_allclose(auc, truth, atol=1e-12)
        checked += 1
    return checked


def oklab(color):
    rgb = np.array([int(color[i:i+2], 16)/255 for i in (1, 3, 5)])
    rgb = np.where(rgb <= 0.04045, rgb/12.92, ((rgb + 0.055)/1.055)**2.4)
    lms = np.cbrt(np.array([[0.4122214708, 0.5363325363, 0.0514459929], [0.2119034982, 0.6806995451, 0.1073969566],
                            [0.0883024619, 0.2817188376, 0.6299787005]]) @ rgb)
    return np.array([[0.2104542553, 0.7936177850, -0.0040720468], [1.9779984951, -2.4285922050, 0.4505937099],
                     [0.0259040371, 0.7827717662, -0.8086757660]]) @ lms


def from_oklab(lab):
    """sRGB hex, or None when outside the gamut."""
    lms = (np.array([[1, 0.3963377774, 0.2158037573], [1, -0.1055613458, -0.0638541728], [1, -0.0894841775, -1.2914855480]]) @ lab)**3
    rgb = np.array([[4.0767416621, -3.3077115913, 0.2309699292], [-1.2684380046, 2.6097574011, -0.3413193965],
                    [-0.0041960863, -0.7034186147, 1.7076147010]]) @ lms
    if (rgb < -1e-6).any() or (rgb > 1 + 1e-6).any():
        return None
    rgb = np.where(rgb <= 0.0031308, 12.92*rgb, 1.055*np.clip(rgb, 0, 1)**(1/2.4) - 0.055)
    return '#%02X%02X%02X' % tuple(int(round(v*255)) for v in np.clip(rgb, 0, 1))


def subtype_colours(classes, colours, band=DARK_BAND):
    """Each class's hue at evenly spaced OKLCH lightness across the band, one step per subtype (sorted).

    A lone subtype keeps its class colour; chroma is reduced where a step leaves the sRGB gamut.
    """
    out = {}
    for parent in sorted(set(classes.values())):
        members = sorted(k for k, v in classes.items() if v == parent)
        if len(members) == 1:
            out[members[0]] = colours[parent]; continue
        _, a, b = oklab(colours[parent])
        for member, lightness in zip(members, np.linspace(*band, len(members))):
            out[member] = next(c for c in (from_oklab(np.array([lightness, a*f, b*f])) for f in np.linspace(1, 0, 101)) if c)
    return out


def compositions(groups, members, top=3):
    """Each group's most common members and their share of its cells, e.g. a community's populations."""
    shares = pd.crosstab(np.asarray(groups), np.asarray(members), normalize='index')
    return {group: ', '.join(f'{t} {v:.0%}' for t, v in row[row > 0].sort_values(ascending=False, kind='stable').head(top).items())
            for group, row in shares.iterrows()}


# ---------------------------------------------------------------- heart annotations and counts

def categorical(h5, column):
    """An obs column stored as categories and codes, as strings."""
    codes = h5[f'obs/{column}/codes'][:]
    if (codes < 0).any():
        raise ValueError(f'{column} has missing values')
    return h5[f'obs/{column}/categories'].asstr()[:][codes]


def source_colours(h5, column):
    """The scanpy colours of an obs column by label (uns/<column>_colors follows the categories)."""
    labels, colours = h5[f'obs/{column}/categories'].asstr()[:], h5[f'uns/{column}_colors'].asstr()[:]
    if len(labels) != len(colours):
        raise ValueError(f'{column} needs one colour per category')
    return {label: colour.upper() for label, colour in zip(labels, colours)}


def recover_counts(x, n_counts):
    """Integer counts and the largest rounding residual from log1p(counts / n_counts * target) values.

    The target cancels: counts = expm1(x) / sum(expm1(x)) * n_counts. Every value must round
    within 0.01, detected values must stay detected and each cell must sum to its n_counts.
    """
    if (n_counts <= 0).any() or not np.array_equal(n_counts, np.round(n_counts)):
        raise ValueError('Transcript counts must be positive integers')
    scaled = np.expm1(np.asarray(x, dtype=np.float64))
    exact = scaled/scaled.sum(axis=1, keepdims=True)*n_counts[:, None]
    counts = np.rint(exact)
    residual = float(np.abs(exact - counts).max())
    if residual > 0.01 or ((counts == 0) & (scaled > 0)).any() or not np.array_equal(counts.sum(axis=1), n_counts):
        raise ValueError('Expression values are not log1p-normalized integer counts')
    return counts.astype(np.int64), residual


def display_coordinates(slices, x, y):
    """Each section centred on its bounding box, with y negated (display y points down, so the atria are at the top)."""
    frame = pd.DataFrame({'slice': slices, 'x': x, 'y': y})
    centre = frame.groupby('slice')[['x', 'y']].agg(lambda v: (v.min() + v.max())/2)
    return (frame.x - frame.slice.map(centre.x)).to_numpy(), (frame.slice.map(centre.y) - frame.y).to_numpy()


def contrast(colour):
    """Contrast ratio of a colour against black."""
    rgb = np.array([int(colour[i:i+2], 16)/255 for i in (1, 3, 5)])
    luminance = np.where(rgb <= 0.04045, rgb/12.92, ((rgb + 0.055)/1.055)**2.4) @ [0.2126, 0.7152, 0.0722]
    return float((luminance + 0.05)/0.05)


def majority(groups, members):
    """Each group's most common member (ties go to the first name)."""
    return pd.crosstab(np.asarray(groups), np.asarray(members)).idxmax(axis=1).to_dict()


def annotation_levels(cells, populations, communities):
    """Labels, colours, names and compositions of the four levels, from the source colours.

    Classes take the colour of their most numerous population; Leiden clusters take the hue of
    theirs, in lightness steps where several clusters share one.
    """
    unknown = sorted(set(cells.cell_population) - CLASSES.keys())
    if unknown:
        raise ValueError(f'Populations without a class: {unknown}')
    populations = {**populations, 'Pericyte': PERICYTE}
    classes = majority(cells.cell_class, cells.cell_population)
    named = lambda level, labels: {k: NAMES[level][k] for k in labels if k in NAMES.get(level, {})}
    present = {level: sorted(set(cells[level])) for level in LEVELS}
    levels = dict(
        cell_class=dict(kind='cell_type', title=LEVELS['cell_class'], labels={k: populations[classes[k]] for k in present['cell_class']},
                        descriptions=named('cell_class', present['cell_class'])),
        cell_population=dict(kind='cell_type', title=LEVELS['cell_population'], labels={k: populations[k] for k in present['cell_population']},
                             descriptions=named('cell_population', present['cell_population'])),
        # The populations each cluster contains, for the table (descriptions would replace the labels in the tooltip).
        leiden=dict(kind='cell_type', title=LEVELS['leiden'], descriptions={},
                    labels=dict(sorted(subtype_colours(majority(cells.leiden, cells.cell_population), populations).items())),
                    composition=compositions(cells.leiden, cells.cell_population)),
        community=dict(kind='region', title=LEVELS['community'], labels={k: communities[k] for k in present['community']},
                       descriptions=named('community', present['community']), composition=compositions(cells.community, cells.cell_population)))
    dark = {f'{label} ({level})': colour for level, v in levels.items() for label, colour in v['labels'].items() if contrast(colour) < MIN_CONTRAST}
    if dark:
        raise ValueError(f'Colours below {MIN_CONTRAST}:1 contrast with the black plot: {dark}')
    return levels


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--raw', type=Path, default=ROOT/'data')
    parser.add_argument('--out', type=Path, default=ROOT/'data/merfish-heart-13pcw')
    parser.add_argument('--batch-cells', type=int, default=16384, help='Maximum cells per batch.')
    parser.add_argument('--batch-entries', type=int, default=12_000_000,
                        help='Maximum stored (cell, detected gene) entries per batch; keeps files well under 95 MB.')
    parser.add_argument('--seed', type=int, default=20261029)
    parser.add_argument('--force', action='store_true', help='Replace an existing output directory.')
    args = parser.parse_args()
    started = time.perf_counter()
    if args.out.exists() and any(args.out.iterdir()) and not args.force:
        raise SystemExit(f'{args.out} is not empty; pass --force to replace it')
    staging = args.out.with_name(args.out.name + '.staging')
    if staging.exists():
        shutil.rmtree(staging)
    (staging/'progressive').mkdir(parents=True)
    (staging/'pseudobulk').mkdir()

    expression = args.raw/EXPRESSION
    with h5py.File(expression, 'r') as h5:
        if h5['X'].attrs.get('encoding-type') != 'array':
            raise ValueError('Expected a dense expression matrix')
        panel = h5['var/_index'].asstr()[:]
        order = np.argsort(panel, kind='stable')  # gene_index follows the sorted symbols
        symbols = panel[order]
        g = len(symbols)
        if len(set(symbols)) != g or not 1 < g < 32768:
            raise ValueError('uint16 doubled ranks require 2 to 32767 unique genes')
        codes = h5['obs/sample_id/codes'][:]
        if (codes < 0).any() or (np.diff(codes) < 0).any():
            raise ValueError('Cells must be grouped by section')
        sections = h5['obs/sample_id/categories'].asstr()[:]
        bounds = {s: (int(np.searchsorted(codes, k)), int(np.searchsorted(codes, k, side='right'))) for k, s in enumerate(sections)}
        slices = [s for s in sections if bounds[s][1] > bounds[s][0]]
        cells = pd.DataFrame({'cell': h5['obs/_index'].asstr()[:], 'slice': sections[codes],
                              **{column: categorical(h5, source) for source, column in OBS.items()}})
        n_counts = h5['obs/n_counts'][:].astype(np.float64)
        x, y = h5['obsm/spatial'][:].T
        cells['complexity'], cells['purity'] = h5['obs/complexity'][:], h5['obs/purity'][:]
        population_colours, community_colours = source_colours(h5, 'populations'), source_colours(h5, 'communities')
        n = len(cells)
        # Cell names are '<cell id>-<section>'; the id is the segmentation's, as in the per-experiment Dryad files.
        cell_id = cells.cell.str.rsplit('-', n=1, expand=True)
        if not ((cell_id[1] == cells.slice).all() and cell_id[0].str.fullmatch(r'\d+').all()):
            raise ValueError("Cell names must be '<cell id>-<section>'")
        if not cells.leiden.str.fullmatch(r'\d{1,2}').all():
            raise ValueError('Leiden clusters must be numbered below 100')
        cells.insert(0, 'cell_index', np.arange(n, dtype=np.uint32))
        cells.insert(3, 'cell_id', cell_id[0].astype(np.int64).to_numpy())
        rx, ry = display_coordinates(cells.slice, x, y)
        for at, (column, values) in enumerate({'rx': rx, 'ry': ry, 'x': x, 'y': y}.items()):
            cells.insert(4 + at, column, values)
        cells.insert(8, 'cell_class', cells.cell_population.map(CLASSES))
        cells['leiden'] = cells.leiden.str.zfill(2)  # so the clusters sort in order
        cells['draw_order'] = np.random.default_rng(args.seed).permutation(n).astype(np.uint32)
        cells['batch'] = ''; cells['batch_row'] = np.zeros(n, np.int32)
        cells['n_counts'] = np.zeros(n, np.int64); cells['n_genes_detected'] = np.zeros(n, np.int64)
        levels = annotation_levels(cells, population_colours, community_colours)

        # Exact ranks, randomly ordered per-section batches and per-section pseudobulk counts.
        label_index = {level: np.unique(cells[level].to_numpy()) for level in LEVELS}
        rng = np.random.default_rng(args.seed)
        detection = np.zeros(g, np.int64); totals = np.zeros(g, np.float64); per_slice = {}
        manifest_files, slice_entries, validated, residual = [], [], 0, 0.0
        def record(path, **extra):
            item = dict(path=str(path.relative_to(staging)), bytes=path.stat().st_size, sha256=sha256(path), **extra)
            if item['bytes'] >= 95_000_000:
                raise ValueError(f'File exceeds conservative GitHub limit: {path}')
            manifest_files.append(item); return item
        for slice_id in slices:
            t = time.perf_counter()
            start, stop = bounds[slice_id]
            m = stop - start
            counts, section_residual = recover_counts(h5['X'][start:stop][:, order], n_counts[start:stop])
            residual = max(residual, section_residual)
            rows, cols = np.nonzero(counts)
            values = counts[rows, cols]
            part = cells.iloc[start:stop]
            cells.loc[part.index, 'n_counts'] = counts.sum(axis=1)
            cells.loc[part.index, 'n_genes_detected'] = (counts > 0).sum(axis=1)
            per_slice[slice_id] = np.bincount(cols, minlength=g)
            detection += per_slice[slice_id]; totals += np.bincount(cols, weights=values, minlength=g)
            indptr, gene_index, rank2, zero2 = rank_rows(rows, cols, values, m, g)
            dense_rows = [(row, counts[row].astype(np.float64)) for row in rng.choice(m, min(m, 32), replace=False)]
            validated += validate_rows(indptr, gene_index, rank2, zero2, dense_rows, g, rng)
            shuffled = rng.permutation(m)
            # Equal-size batches, small enough for both the cell and stored-entry limits.
            size = max(1, min(args.batch_cells, int(args.batch_entries*m/max(1, len(gene_index)))))
            size = -(-m // -(-m // size))
            batches = []
            for b, first in enumerate(range(0, m, size)):
                chosen = shuffled[first:first+size]
                lengths = indptr[chosen+1] - indptr[chosen]
                ptr = np.r_[0, np.cumsum(lengths)]
                gather = np.repeat(indptr[chosen] - ptr[:-1], lengths) + np.arange(ptr[-1])
                path = staging/'progressive'/slice_id/f'{b:05d}.npz'; path.parent.mkdir(exist_ok=True)
                global_rows = (start + chosen).astype(np.uint32)
                np.savez_compressed(path, cell_index=global_rows, indptr=ptr.astype(np.uint32),
                                    gene_index=gene_index[gather], rank2=rank2[gather], zero_rank2=zero2[chosen])
                with np.load(path, allow_pickle=False) as z:  # round-trip check of the serialized batch
                    np.testing.assert_array_equal(z['gene_index'], gene_index[gather])
                    np.testing.assert_array_equal(z['rank2'], rank2[gather])
                cells.loc[global_rows, 'batch'] = str(path.relative_to(staging))
                cells.loc[global_rows, 'batch_row'] = np.arange(len(chosen), dtype=np.int32)
                batches.append(record(path, slice=slice_id, n_cells=len(chosen))['path'])
            arrays = {}
            for level in LEVELS:
                present, n_cells, sums = pseudobulk(rows, cols, values, np.searchsorted(label_index[level], part[level].to_numpy()), g)
                np.testing.assert_array_equal(sums.sum(axis=0, dtype=np.float64), np.bincount(cols, weights=values, minlength=g))
                arrays.update({f'{level}_index': present.astype(np.uint16), f'{level}_cells': n_cells.astype(np.uint32),
                               f'{level}_counts': sums})
            path = staging/'pseudobulk'/f'{slice_id}.npz'
            np.savez_compressed(path, **arrays)
            record(path, slice=slice_id, layout='pseudobulk')
            slice_entries.append(dict(id=slice_id, n_cells=m, cell_index_start=start, cell_index_stop=stop,
                                      bbox_rxry=[float(part.rx.min()), float(part.ry.min()), float(part.rx.max()), float(part.ry.max())],
                                      batches=batches, pseudobulk=str(path.relative_to(staging)), rank_seconds=round(time.perf_counter()-t, 2)))
            print(f'{slice_id}: {m:,} cells ranked into {len(batches)} batches in {time.perf_counter()-t:.1f}s', flush=True)

    genes = pd.DataFrame({'gene_index': np.arange(g, dtype=np.uint16), 'symbol': symbols,
                          'n_cells_detected': detection, 'total_counts': totals})
    (staging/'annotations.json').write_text(json.dumps(levels, indent=1, ensure_ascii=False)+'\n'); record(staging/'annotations.json')
    cells.to_parquet(staging/'cells.parquet', index=False); record(staging/'cells.parquet')
    genes.to_parquet(staging/'genes.parquet', index=False); record(staging/'genes.parquet')
    pd.DataFrame({'gene_index': np.tile(np.arange(g, dtype=np.uint16), len(slices)), 'slice': np.repeat(slices, g),
                  'n_cells_detected': np.concatenate([per_slice[s] for s in slices])}).to_parquet(
        staging/'gene_detection.parquet', index=False); record(staging/'gene_detection.parquet')
    manifest = dict(
        format='sparse-doubled-average-ranks-v1', layout='per-slice-random-batches-v1', dataset=DATASET,
        created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        n_cells=n, n_genes=g, n_slices=len(slices), seed=args.seed, batch_cells=args.batch_cells, batch_entries=args.batch_entries,
        observations='MERFISH cells of three frontal sections of one 13 p.c.w. human heart, as annotated in overall_merfish.h5ad; '
                     'counts are decoded transcripts per cell',
        coordinates=dict(columns='rx, ry', units='micrometres; each section centred on its bounding box, ry = -(y - centre)', um_per_unit=1.0,
                         raw='x, y are obsm/spatial (micrometres, y up), where the source sets the sections side by side'),
        counts=dict(method='expm1(X) / sum(expm1(X)) * n_counts per cell, rounded; X is log1p(counts / n_counts * 365)',
                    max_rounding_residual=residual),
        region_assignment='cellular community from the source: k-means (k = 13) of the population mix within 150 µm of each cell',
        sources={'expression': dict(path=EXPRESSION, sha256=sha256(expression)),
                 'paper': 'Farah et al., Nature 627, 854-864 (2024), doi:10.1038/s41586-024-07171-z',
                 'data': 'https://doi.org/10.5061/dryad.w0vt4b8vp'},
        ranks=dict(direction='ascending: larger rank means higher expression', ties='average', scale=2,
                   zeros='all unstored genes have rank zero_rank2 / 2; stored ranks are absolute',
                   universe='every MERFISH panel gene, under its panel symbol',
                   auroc='(sum(rank of set) - k*(k+1)/2) / (k*(G-k))'),
        annotation_levels={k: dict(kind=v['kind'], title=v['title'], n_labels=len(v['labels'])) for k, v in levels.items()},
        pseudobulk=dict(description='per section and annotation level: summed panel-gene counts per (label, gene), '
                                    '<level>_counts rows follow <level>_index into labels[level]; <level>_cells counts cells',
                        labels={level: label_index[level].tolist() for level in LEVELS}, cells={level: 'all' for level in LEVELS}),
        validation=dict(exact_rank_and_pairwise_auc_rows=validated, batches_roundtripped=True),
        software=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__,
                      h5py=h5py.__version__),
        slices=slice_entries, files=manifest_files, seconds=round(time.perf_counter()-started, 1))
    (staging/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    if args.out.exists():
        shutil.rmtree(args.out)
    staging.rename(args.out)
    print(json.dumps({k: v for k, v in manifest.items() if k not in ['files', 'slices', 'pseudobulk']}, indent=2))
    print('OUTPUT_BYTES', sum(p.stat().st_size for p in args.out.rglob('*') if p.is_file()), flush=True)


if __name__ == '__main__':
    main()
