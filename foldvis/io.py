from io import StringIO
from pathlib import Path
from typing import Union

from Bio.PDB import PDBIO, Structure


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