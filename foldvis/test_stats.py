from pathlib import Path

from foldvis.models import Fold
from foldvis.stats import get_alpha_carbon_coords



def test_get_alpha_carbon_coords():
    '''
    https://stackoverflow.com/questions/32527861/python-unit-test-that-uses-an-external-data-file
    '''
    here = Path(__file__).parent

    rel = 'data/1AAY_alphafold/test_676a7_unrelaxed_rank_1_model_2.pdb'
    p = here.parent / rel
    model = Fold(p)
    assert \
    len(model.sequence) == len(list(get_alpha_carbon_coords(model, 'CA')))