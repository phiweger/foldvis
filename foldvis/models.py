from copy import deepcopy
import io
import json
from pathlib import Path
import re
from typing import Union

from Bio import PDB
from Bio.PDB.Structure import Structure
import numpy as np

from foldvis.utils import align_structures
from foldvis.io import save_pdb


class Fold():
    def __init__(self, fp, quiet=True):
        self.path = Path(fp)
        if not quiet:
            print(f'Loading structure in {self.path.name}')
        self.structure = self.read_pdb(self.path)
        self.transformed = False
        self.annotation = {}

        ln = len(list(self.structure.get_residues()))
        self.annotate_('position', [i+1 for i in range(ln)])
        return None

    def __repr__(self):
        return 'Fold loaded from ' + self.path.name

    def __len__(self):
        '''Number of amino acids in the sequence'''
        return len(list(self.structure.get_residues()))

    def read_pdb(self, fp: Union[str, Path], name: str='x') -> Structure:
        fp = Path(fp)
        assert fp.exists()
        pdb_parser = PDB.PDBParser(QUIET=True)
        structure = pdb_parser.get_structure(name, str(fp))
        return structure

    def align_to(self, ref, mode=0, minscore=0.5):
        tmscore, rot, tra = align_structures(ref.path, self.path, mode=mode, minscore=minscore)
        cp = deepcopy(self)
        cp.structure.transform(rot, tra)
        cp.transformed = True
        return tmscore, cp

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

    # def add_scores(self, d):
    #     for k, v in d.items():
    #         self.scores[k] = v        
    #     return None

    def to_stream(self):
        stream = io.StringIO()
        return save_pdb(self.structure, stream).getvalue()

    def annotate_(self, key, values, check=False):
        if check:
            assert len(values) == len(self)
        self.annotation[key] = values
        return None


class AlphaFold():
    '''
    fold = AlphaFold(...)
    fold.best -> structure, wrapped in Fold object (keeps track of filepath ...)
    '''
    def __init__(self, indir, workdir=None):
        self.models = {}
        self.workdir = workdir

        for n, i in enumerate(self.read_alphafold(indir, workdir)):
            self.models[n+1] = i
        
        return None


    def read_alphafold(self, filedir, outdir=None):
        '''
        Load structures and (quality) scores, align them, and put prepare data
        for visualization.
        '''
        files = Path(filedir).glob('*.pdb')
        d = {}
    
        # Load scores
        for i in files:
            model = int(re.match(r'.*model_(\d).pdb', i.name).group(1))
            # print(f'Loading model {model}')
            fp = str(i.resolve())
            fp = fp.replace(f'{model}.pdb', f'{model}_scores.json')
            
            with open(fp, 'r') as file:
                scores = json.load(file)
           
            fold = Fold(i)
            _ = fold.annotate_('plddt', scores['plddt'])
    
            v = np.mean(scores['plddt'])
            d[fold] = v
    
        # Rank models by pLDDT, best is reference
        ref, *queries = [i for i, j in sorted(d.items(), key=lambda x: x[1], reverse=True)]

        print(f'Best model (pLDDT): {ref.path.name}')
        print(f'Align remaining models to best and rename')
        # Align into the same space
        rest = []
        for qry, chain in zip(queries, 'BCDE'):
            _, trx = qry.align_to(ref)
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


