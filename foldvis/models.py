from copy import copy
from pathlib import Path
from typing import Union

from Bio import PDB
from Bio.PDB.Structure import Structure

from foldvis.io import read_alphafold


class Fold():
    def __init__(self, fp):
        self.structure = self.read_pdb(fp)
        self.path = Path(fp)
        self.transformed = False
        self.scores = {}
        return None

    def __repr__(self):
        return 'Fold loaded from ' + self.path.name

    def read_pdb(self, fp: Union[str, Path], name: str='x') -> Structure:
        fp = Path(fp)
        assert fp.exists()
        pdb_parser = PDB.PDBParser(QUIET=True)
        structure = pdb_parser.get_structure(name, str(fp))
        return structure

    def align_to(self, ref):
        rot, tra = align_structures(ref.path, self.path)
        cp = copy(self)
        cp.structure.transform(rot, tra)
        cp.transformed = True
        return cp

    def rename_chains_(self, renames: dict) -> None:
        '''
        - https://stackoverflow.com/questions/70246451/how-do-i-change-the-chain-name-of-a-pdb-file
        - modifies in place, "_" suffix convention like in pytorch
        '''
        for model in self.structure:
            for chain in model:
                old_name = chain.get_id()
                new_name = renames.get(old_name)
                if new_name:
                    print(f'Renaming chain {old_name} to {new_name}')
                    chain.id = new_name
                else:
                    print(f'Keeping chain name {old_name}, no new name found')
        return None

    def add_scores(self, d):
        self.scores = d
        return None


class AlphaFold():
    '''
    fold = AlphaFold(...)
    fold.best -> structure, wrapped in Fold object (keeps track of filepath ...)
    '''
    def __init__(self, indir, workdir=None):
        self.models = {}
        self.workdir = workdir

        for n, i in enumerate(read_alphafold(indir, workdir)):
            self.models[n+1] = i
        
        return None
