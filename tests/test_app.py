from pathlib import Path
import json
import numpy as np
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest
from spatial_app.enrichment import MIN_CELLS

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'data/merfish-heart-13pcw'
pytestmark = pytest.mark.skipif(not (DATA/'manifest.json').exists(), reason='prepared heart dataset not available')


@pytest.fixture
def manifest():
    return json.loads((DATA/'manifest.json').read_text())


def test_three_sections_score_in_one_run_and_fill_both_tables(manifest):
    app = AppTest.from_file(str(ROOT/'app.py'), default_timeout=180).run()
    assert not app.exception
    assert manifest['n_slices'] == 3 and len(app.session_state['scores']) == manifest['n_cells']
    # The first run scores every cell of every batch, with no progress bar or later reruns.
    assert np.isfinite(app.session_state['scores']).all() and not app.get('progress')
    assert app.session_state['done'] == {b for s in manifest['slices'] for b in s['batches']}
    header = next(m.value for m in app.markdown if 'frontal sections' in m.value)
    assert header == '**Farah-heart-13pcw · 3 frontal sections of one 13 p.c.w. heart**'
    panels = app.session_state['cell_plot_transport']['geometry']['panels']
    assert [p['label'] for p in panels] == ['R77_4C4', 'R78_4C12', 'R78_4C15']
    assert app.session_state['job']['query'].universe_size == manifest['n_genes'] == 238
    cell_table, community_table = app.dataframe[0].value, app.dataframe[1].value
    cells = pd.read_parquet(DATA/'cells.parquet')
    for table, column in [(cell_table, 'cell_population'), (community_table, 'community')]:
        counts = cells[column].value_counts()
        assert table.set_index(table.columns[0])['Cells'].to_dict() == counts[counts >= MIN_CELLS].to_dict()
    assert cell_table.columns[:2].tolist() == ['Cell population', 'Cell class']
    assert community_table.columns[:2].tolist() == ['Cellular community', 'Most common populations']
    # The default set marks the His-Purkinje cells of the ventricular conduction system.
    assert app.session_state['job']['query'].target_symbols == ('HCN4', 'IRX1', 'IRX2', 'IRX3', 'TBX3')
    matching = next(e for e in app.expander if e.label.startswith('Gene matching'))
    assert matching.label == 'Gene matching · 5 targets used in a 238-gene universe' and not matching.warning
    assert cell_table.iloc[0, 0] == 'vCM-His-Purkinje' and community_table.iloc[0, 0] == 'VCS'
    app.selectbox(key='cell_level').set_value('leiden').run()
    leiden_table = app.dataframe[0].value
    assert leiden_table.columns[:2].tolist() == ['Leiden cluster', 'Most common populations']
    assert leiden_table.iloc[0, 1] == 'vCM-His-Purkinje 100%'


def test_mouse_input_maps_to_human_and_colour_modes_run(manifest):
    app = AppTest.from_file(str(ROOT/'app.py'), default_timeout=180).run()
    app.selectbox(key='input_species').set_value('Mouse').run()
    before = app.session_state['scores'].copy()
    app.text_area(key='target_input').set_value('Myh6 Nr2f2 Pln NOT_A_GENE').run()
    assert not app.exception
    assert app.session_state['job']['query'].target_symbols == ('MYH6', 'NR2F2')
    # A new gene set is rescored, and its tables filled, in the same run.
    scores = app.session_state['scores']
    assert np.isfinite(scores).all() and not np.array_equal(scores, before) and len(app.dataframe) == 2
    assert app.session_state['mapping_unmapped'] == (['NOT_A_GENE', 'Pln'], [])  # Pln has no human entry in the HomoloGene snapshot
    app.selectbox(key='cell_level').set_value('cell_class').run()
    assert app.dataframe[0].value.iloc[0, 0] == 'aCM'  # atrial markers
    for mode in ['Cell type', 'Region', 'Cell type AUROC', 'Region AUROC']:
        app.radio(key='colour_mode').set_value(mode).run()
        assert not app.exception
