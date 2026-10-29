import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from scipy.stats import rankdata

pytest.importorskip('h5py')  # preprocessing needs requirements-prep.txt; the app alone does not
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from prepare_heart_cells import (CLASSES, DARK_BAND, MIN_CONTRAST, NAMES, PERICYTE, annotation_levels, contrast,  # noqa: E402
                                 display_coordinates, oklab, pseudobulk, rank_rows, recover_counts, subtype_colours)


def test_rank_rows_matches_rankdata_with_ties_zeros_and_empty_rows():
    rng = np.random.default_rng(3)
    g, n = 40, 25
    dense = rng.integers(0, 4, (n, g)) * (rng.random((n, g)) < .5)
    dense[0] = 0                    # no detected genes
    dense[1] = rng.integers(1, 3, g)  # no zeros
    rows, cols = np.nonzero(dense)
    shuffle = rng.permutation(len(rows))  # aggregation output is unordered
    indptr, gene_index, rank2, zero2 = rank_rows(rows[shuffle], cols[shuffle], dense[rows, cols][shuffle], n, g)
    for row in range(n):
        decoded = np.full(g, zero2[row]/2)
        a, b = indptr[row], indptr[row+1]
        assert (np.diff(gene_index[a:b].astype(int)) > 0).all()
        decoded[gene_index[a:b]] = rank2[a:b]/2
        np.testing.assert_array_equal(decoded, rankdata(dense[row], method='average'))
    assert zero2[1] == 0 and zero2[0] == g+1


def test_pseudobulk_sums_counts_per_present_label():
    rng = np.random.default_rng(5)
    dense = rng.poisson(1.5, (12, 7))
    codes = np.array([4, 4, 1, 4, 9, 1, 1, 4, 9, 9, 4, 1])  # global label indices; 0, 2, 3 ... absent
    rows, cols = np.nonzero(dense)
    present, n_cells, sums = pseudobulk(rows, cols, dense[rows, cols], codes, 7)
    assert present.tolist() == [1, 4, 9] and n_cells.tolist() == [4, 5, 3] and sums.dtype == np.uint32
    for i, label in enumerate(present):
        np.testing.assert_array_equal(sums[i], dense[codes == label].sum(axis=0))


def test_recover_counts_inverts_scanpy_normalization():
    rng = np.random.default_rng(0)
    counts = rng.poisson(0.8, (50, 30))*rng.integers(1, 4, (50, 30))
    counts[:, 0] += 1  # every cell has counts
    n = counts.sum(axis=1).astype(np.float64)
    x = np.log1p(counts/n[:, None]*np.median(n)).astype(np.float32)  # normalize_total, then log1p
    recovered, residual = recover_counts(x, n)
    np.testing.assert_array_equal(recovered, counts)
    assert residual < 1e-3
    with pytest.raises(ValueError, match='not log1p-normalized'):
        recover_counts(x, n + 1)  # wrong totals give fractional counts
    with pytest.raises(ValueError, match='positive integers'):
        recover_counts(x, n + 0.5)


def test_sections_are_centred_with_the_atria_up():
    slices = np.array(['A', 'A', 'B', 'B'])
    rx, ry = display_coordinates(slices, np.array([0., 10, 100, 120]), np.array([5., 25, -50, -10]))
    np.testing.assert_allclose(rx, [-5, 5, -10, 10])
    np.testing.assert_allclose(ry, [10, -10, 20, -20])  # source y is up, display y down


def test_classes_group_the_27_populations_as_in_fig1():
    assert len(CLASSES) == 27 and len(set(CLASSES.values())) == 11
    assert CLASSES['Pericyte'] == CLASSES['VSMC'] == 'SMC' and CLASSES['VIC'] == 'Fibro'
    assert set(NAMES['cell_population']) <= CLASSES.keys() and set(NAMES['cell_class']) <= set(CLASSES.values())


def test_subtypes_share_their_class_hue_in_rising_lightness_steps():
    colours = subtype_colours({'a': 'X', 'b': 'X', 'c': 'X', 'd': 'Y'}, {'X': '#3987E5', 'Y': '#D95926'})
    assert colours['d'] == '#D95926'  # a lone subtype keeps its class colour
    labs = np.array([oklab(colours[k]) for k in 'abc'])
    np.testing.assert_allclose(labs[:, 0], np.linspace(*DARK_BAND, 3), atol=0.005)
    assert np.ptp(np.degrees(np.arctan2(labs[:, 2], labs[:, 1]))) < 3  # same hue family as the class


def heart_cells():
    populations = ['Pericyte']*3 + ['VSMC']*2 + ['vCM-LV-Compact']*4 + ['vCM-RV-Compact']
    return pd.DataFrame({'cell_population': populations, 'cell_class': [CLASSES[p] for p in populations],
                         'leiden': ['00']*3 + ['01']*2 + ['02']*2 + ['03']*3, 'community': ['Outer-LV']*6 + ['VCS']*4})


SOURCE = {'Pericyte': '#000000', 'VSMC': '#276834', 'vCM-LV-Compact': '#7382BC', 'vCM-RV-Compact': '#BDC1D5'}
COMMUNITIES = {'Outer-LV': '#F37121', 'VCS': '#8B58A4'}


def test_levels_keep_the_source_colours_and_stand_out_on_black():
    levels = annotation_levels(heart_cells(), SOURCE, COMMUNITIES)
    populations = levels['cell_population']
    assert populations['labels'] == {**SOURCE, 'Pericyte': PERICYTE}  # black in the source, like the plot
    assert populations['descriptions'] == {k: NAMES['cell_population'][k] for k in ['VSMC', 'vCM-LV-Compact', 'vCM-RV-Compact']}
    assert levels['cell_class']['labels'] == {'SMC': PERICYTE, 'vCM': '#7382BC'}  # the most numerous population's colour
    # A cluster alone in its population keeps its colour; two that share one step through its hue in lightness.
    leiden = levels['leiden']['labels']
    assert leiden['00'] == PERICYTE and leiden['01'] == '#276834'
    labs = np.array([oklab(leiden[k]) for k in ['02', '03']] + [oklab('#7382BC')])
    hue = np.degrees(np.arctan2(labs[:, 2], labs[:, 1]))
    assert labs[0, 0] < labs[1, 0] and np.ptp(hue) < 3
    assert levels['leiden']['composition']['00'] == 'Pericyte 100%'  # no zero shares listed
    assert levels['leiden']['composition']['03'] == 'vCM-LV-Compact 67%, vCM-RV-Compact 33%'
    assert levels['community']['composition'] == {'Outer-LV': 'Pericyte 50%, VSMC 33%, vCM-LV-Compact 17%',
                                                  'VCS': 'vCM-LV-Compact 75%, vCM-RV-Compact 25%'}
    assert levels['community']['descriptions'] == {'Outer-LV': 'Outer left ventricular wall', 'VCS': 'Ventricular conduction system'}
    assert all(contrast(c) >= MIN_CONTRAST for level in levels.values() for c in level['labels'].values())


def test_levels_reject_dark_colours_and_unknown_populations():
    with pytest.raises(ValueError, match='contrast'):
        annotation_levels(heart_cells(), {**SOURCE, 'VSMC': '#202020'}, COMMUNITIES)
    with pytest.raises(ValueError, match='without a class'):
        annotation_levels(heart_cells().assign(cell_population='Myocyte'), SOURCE, COMMUNITIES)
