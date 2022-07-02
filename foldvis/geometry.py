# import shutil
import subprocess
import tempfile

import numpy as np
import screed


def is_close(pos, fold, radius):
    '''
    fold = Fold('test_676a7_unrelaxed_rank_1_model_2.pdb')
    list(is_close(1, fold, 10))
    # [True, True, True, True, True, False, False, ...]
    '''
    residues = list(fold.structure.get_residues())
    ref = residues[pos]
    a = ref.center_of_mass()

    for res in residues:
        b = res.center_of_mass()
        # https://stackoverflow.com/questions/1401712/how-can-the-euclidean-distance-be-calculated-with-numpy
        dist = np.linalg.norm(a - b)
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
    print(fp)

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
