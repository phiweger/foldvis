from collections import defaultdict
from io import StringIO
from pathlib import Path
from typing import Union

from Bio import PDB, SeqUtils
from Bio.PDB import PDBIO, Structure
from Bio.PDB.PDBParser import PDBParser
import pandas as pd
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


def read_pdb(fp: Union[str, Path], name: str='x') -> Structure:
   '''
   # https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ

   p = PDBParser()
   structure = p.get_structure("X", "pdb1fat.ent")
   for model in structure:
       for chain in model:
           for residue in chain:
               for atom in residue:
                   print(atom)
   '''
   fp = Path(fp)
   assert fp.exists()
   # https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
   pdb_parser = PDB.PDBParser(QUIET=True, PERMISSIVE=0)
   structure = pdb_parser.get_structure(name, str(fp))
   return structure


def load_bfactor_column(fp):
    '''
    Load annotation data stored in the bfactor column of a .pdb file
    '''
    parser = PDBParser()
    structure = parser.get_structure('', fp)
    d = {}

    for res in structure.get_residues():
        for atom in res:
            # ALA > Ala
            x = atom.parent.resname[0] + atom.parent.resname[1:].lower()  
            try:
                aa = SeqUtils.IUPACData.protein_letters_3to1[x]  
                # same value for all atoms in a resudue
            except KeyError:
                continue

            num = atom.full_id[3][1]
            d[num] = [atom.bfactor, aa]
    
    return [i for i, j in d.values()]


