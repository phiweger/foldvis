from itertools import combinations
from math import log, e
from pathlib import Path
import subprocess
from typing import Union
import tempfile

from Bio.PDB import PDBIO
from Bio.PDB.Structure import Structure

import numpy as np


def chunks(l, n):
    '''
    Yield successive n-sized chunks from l (stackoverflow, 312443).

    a = [1, 2, 3, 4]
    list(chunks(a, 2))
    # [[1, 2], [3, 4]]

    Returns empty list if list empty.

    For overlapping chunks, see windows()
    '''
    for i in range(0, len(l), n):
        yield l[i:i + n]


def rmdir(directory):
    '''
    DEPRECATED, using tempfile ".cleanup()" fn
    https://stackoverflow.com/questions/13118029/deleting-folders-in-python-recursively/49782093#49782093
    '''
    directory = Path(directory)
    for item in directory.iterdir():
        if item.is_dir():
            try:
                rmdir(item)
            except NotADirectoryError:
                # NotADirectoryError: [Errno 20] Not a directory:
                # '.../tmp/latest' -- a symlink (?) created by foldseek
                item.unlink()
        else:
            item.unlink()
    directory.rmdir()


def align_structures(target: Structure, query: Structure, mode=0, minscore=0.5):
    '''
    PyMOL just wraps programs as well:

    - https://pymolwiki.org/index.php/TMalign
    - https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/tmalign.py

    https://janakiev.com/blog/python-shell-commands/

    https://stackoverflow.com/questions/17742789/running-multiple-bash-commands-with-subprocess

    conda create -n foldseek -c conda-forge -c bioconda foldseek
    conda activate foldseek
    # wget https://mmseqs.com/foldseek/foldseek-osx-universal.tar.gz; tar xvzf foldseek-osx-universal.tar.gz; export PATH=$(pwd)/foldseek/bin/:$PATH

    To superimpose two structures we have to first find out how to transform 
    one into the other:

    https://github.com/steineggerlab/foldseek/issues/13

    mode:

    --cov-mode INT

    0: coverage of query and target
    1: coverage of target
    2: coverage of query
    3: target seq. length has to be at least x% of query length
    4: query seq. length has to be at least x% of target length
    5: short seq. needs to be at least x% of the other seq. length [0]

    --tmscore-threshold [0.5]
    accept alignments with a tmsore > thr [0.0,1.0] [0.500]
    '''
    tmp = tempfile.TemporaryDirectory()
    p = tmp.name

    print(f'Mode {mode}, min. Tm score {minscore}')

    steps = [
        f'foldseek createdb {target} {p}/targetDB',
        f'foldseek createdb {query} {p}/queryDB',
        f'foldseek search {p}/queryDB {p}/targetDB {p}/aln {p}/tmp -a --cov-mode {mode} --tmscore-threshold {minscore}',
        f'foldseek aln2tmscore {p}/queryDB {p}/targetDB {p}/aln {p}/aln_tmscore',
        f'foldseek createtsv {p}/queryDB {p}/targetDB {p}/aln_tmscore {p}/aln_tmscore.tsv'
    ]

    command = '; '.join(steps)
    log = subprocess.run(command, capture_output=True, shell=True)

    with open(f'{p}/aln_tmscore.tsv', 'r') as file:
        qry, rest = next(file).strip().split('\t')
        ref, score, *rest = rest.split(' ')
        print(f'Tm-score: {score}')
        rest = [float(i) for i in rest]
        translation = rest[:3]
        rotation = list(chunks(rest[3:], 3))

    tmp.cleanup()

    return np.array(rotation), np.array(translation)


def transform_(structure: Structure, translation: np.ndarray, rotation: np.ndarray):
    '''
    DEPRECATED, can use "structure.transform(rotation, translation)"
    '''
    for i in structure.get_atoms():
        i.transform(rotation, translation)
    return None


def entropy(labels, base=None):
    '''
    Computes entropy of label distribution.

    https://stackoverflow.com/questions/15450192/fastest-way-to-compute-entropy-in-python

    nextstrain uses entropy:

    https://docs.nextstrain.org/en/latest/guides/share/download-data.html#diversity-entropy-data
    '''

    n_labels = len(labels)

    if n_labels <= 1:
      return 0

    value,counts = np.unique(labels, return_counts=True)
    probs = counts / n_labels
    n_classes = np.count_nonzero(probs)

    if n_classes <= 1:
      return 0

    ent = 0.

    # Compute entropy
    base = e if base is None else base
    for i in probs:
      ent -= i * log(i, base)

    return ent


def mean_pairwise_similarity(labels):
    '''
    Geneious calculates "Identity" as "Mean pairwise identity over all pairs
    in the column." (hover over "Identity" after loading a MSA).

    See also discussions here:

    - https://www.biostars.org/p/3856/
    - https://www.biostars.org/p/5067/
    - https://www.biostars.org/p/223824/

    For en extensive review:

    - https://onlinelibrary.wiley.com/doi/abs/10.1002/prot.10146
    '''
    l = [1 if i == j else 0 for i, j in combinations(labels, 2)]
    return round(sum(l) / len(l), 4)


