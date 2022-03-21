from collections import defaultdict
from io import StringIO
from pathlib import Path
from typing import Union

from Bio.PDB import PDBIO, Structure
import screed

from foldvis.utils import entropy, mean_pairwise_similarity


def save_pdb(structure: Structure, out: Union[str, Path, StringIO]) -> None:
    file = PDBIO()
    file.set_structure(structure)
    
    # Save to stream
    if type(out) == StringIO:
        file.save(out)
        return out

    # Save to file
    else:
        p = Path(out)
        if not p.parent.exists():
            p.parent.mkdir(parents=True)
        file.save(p.resolve().__str__())
        return None


def load_conserved(fp, ref, metric=mean_pairwise_similarity):
    variance = defaultdict(list)
    cnt = 0
    with screed.open(fp) as file:
        for i in file:
            if ref in i.name:
                ix = cnt
            
            for pos, aa in enumerate(i.sequence):
                variance[pos].append(aa)
            cnt += 1
        
    l = []
    for pos, residues in variance.items():
        ref_aa = residues[ix]
        if not ref_aa == '-':
            l.append(metric(residues))

    return l