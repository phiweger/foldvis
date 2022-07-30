# import shutil
import subprocess
import tempfile
from typing import Union

from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
import numpy as np
import screed


def get_alpha_carbon_atoms(fold, only_coords=False):
    '''
    alpha carbon: https://foldit.fandom.com/wiki/Alpha_carbon
    '''
    for res in fold.structure.get_residues():
        l = []
        for i in res.get_atoms():
            # There is only one alpha carbon per amino acid.
            if i.get_id() == 'CA':
                if only_coords:
                    yield i.coord
                else:
                    yield i


def get_coordinate(x: Union[Atom, Residue]):
    '''
    Calculating center of mass is much slower than looking up the coordinate
    of a carbon atom.
    '''
    if type(x) == Residue:
        return x.center_of_mass()
    elif type(x) == Atom:
        return x.coord
    else:
        raise ValueError('Unsupported type')


def euclidean_distance(a, b):
    '''
    https://stackoverflow.com/questions/1401712/how-can-the-euclidean-distance-be-calculated-with-numpy
    '''
    return np.linalg.norm(a-b)


def is_close(pos, fold, radius, coordinates='alpha_carbons'):
    '''
    fold = Fold('test_676a7_unrelaxed_rank_1_model_2.pdb')
    list(is_close(1, fold, 10))
    # [True, True, True, True, True, False, False, ...]

    coordinates .. center_of_mass, alpha_carbon

    TODO: pseudo single-atom representation of side chains:

    > Specifically, we defined this distance according to the sites' side chain
    center of masses. A consequence of approximating DTL with respect to the
    closest ligand-binding sites is that by definition, any ligand-binding
    residue has a DTL of 0. -- Kiefl et al.,
    https://www.biorxiv.org/content/10.1101/2022.03.02.482602v1.full.pdf

    - https://pymolwiki.org/index.php/Sidechaincenters
    - https://bioinformatics.stackexchange.com/questions/18162/pymol-python-script-for-selecting-a-residues-sidechain-and-calculating-its-cent
    '''

    if coordinates == 'center_of_mass':
        chain = list(fold.structure.get_residues())
    
    elif coordinates == 'alpha_carbons':
        chain = list(get_alpha_carbon_atoms(fold))
    
    a = get_coordinate(chain[pos])

    for i in chain:
        b = get_coordinate(i)
        # https://stackoverflow.com/questions/1401712/how-can-the-euclidean-distance-be-calculated-with-numpy
        # dist = np.linalg.norm(a - b)
        dist = euclidean_distance(a, b)
        if dist < radius:
            yield True
        else:
            yield False


def get_foldseek_vae_states(fold):
    '''
    https://github.com/steineggerlab/foldseek/issues/15
    '''
    tmp = tempfile.TemporaryDirectory()
    p = tmp.name
    fp = fold.path.resolve().__str__()

    steps = [
        f'foldseek createdb {fp} {p}/db',
        f'foldseek lndb {p}/db_h {p}/db_ss_h',
        f'foldseek convert2fasta {p}/db_ss {p}/db_ss.fasta',
    ]
    command = '; '.join(steps)
    log = subprocess.run(command, capture_output=True, shell=True)
    assert log.returncode == 0, log.stderr

    # shutil.copyfile(f'{p}/db_ss.fasta', outfile)
    with screed.open(f'{p}/db_ss.fasta') as file:
       return next(file).sequence


def distance_to_closest_active_site(fold, binding_frequencies, threshold=0.5):
    '''
    Usage:

    f = Fold('serine_hydroxymethyltransferase.pdb')
    b = Binding(f, 'confident')
    b.predict_binding_(pfam)
    bf = b.get_binding('PF00464.18', 'SER')
    distance_to_closest_active_site(f, bf, .5)
    '''
    residues = list(fold.structure.get_residues())
    bf = binding_frequencies

    assert len(residues) == len(bf) 
    active = [r for r, f in zip(residues, bf) if f >= threshold]

    l = []
    for r in residues:
        rm = r.center_of_mass()
        d = np.min([euclidean_distance(rm, a.center_of_mass()) for a in active])
        l.append(float(d))
    
    return l

