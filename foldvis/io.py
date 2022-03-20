import json
from pathlib import Path
from typing import Union, Generator
import re

from Bio import PDB
from Bio.PDB.Structure import Structure

from foldvis.models import Fold


def read_pdb(fp: Union[str, Path], name: str='x') -> Structure:
    fp = Path(fp)
    assert fp.exists()
    
    pdb_parser = PDB.PDBParser(QUIET=True)
    structure = pdb_parser.get_structure(name, str(fp))
    return structure


def read_alphafold(filedir, outdir=None):

    files = Path(filedir).glob('*.pdb')
    d = {}

    # Load scores
    for i in files:
        model = int(re.match(r'.*model_(\d).pdb', i.name).group(1))
        fp = str(i.resolve()).replace(f'{model}.pdb', f'{model}_scores.json')
        
        with open(fp, 'r') as file:
            scores = json.load(file)
       
        fold = Fold(i)
        _ = fold.add_scores(scores)

        v = np.mean(scores['plddt'])
        d[fold] = v

    # Rank models by pLDDT, best is reference
    ref, *queries = [i for i, j in sorted(d.items(), key=lambda x: x[1], reverse=True)]

    # Align into the same space
    rest = []
    for qry, chain in zip(queries, 'BCDE'):
        trx = qry.align_to(ref)
        trx.rename_chains_({'A': chain})
        rest.append(trx)

    if outdir:
        outdir = Path(outdir)
        
        if not outdir.exists():
            outdir.mkdir(parents=True)

        name = ref.path.name.replace('.pdb', '.reference.pdb')
        _ = save_pdb(ref.structure, outdir / name)

        for qry in rest:
            name = qry.path.name.replace('.pdb', '.transform.pdb')
            _ = save_pdb(qry.structure, outdir / name)

    return [ref] + rest


def save_pdb(structure: Structure, out: Union[str, Path]) -> None:
    io = PDBIO()
    io.set_structure(structure)
    
    p = Path(out)
    if not p.parent.exists():
        p.parent.mkdir(parents=True)

    io.save(p.resolve().__str__())
    return None