import numpy as np


def is_close(pos, fold, radius):
    '''
    fold = Fold('test_676a7_unrelaxed_rank_1_model_2.pdb')
    list(is_close(1, fold, 10))
    # [True, True, True, True, True, False, False, ...]
    '''
    residues = list(fold.structure.get_residues())
    ref = residues[pos]

    for res in residues:
        a = ref.center_of_mass()
        b = res.center_of_mass()
        # https://stackoverflow.com/questions/1401712/how-can-the-euclidean-distance-be-calculated-with-numpy
        dist = np.linalg.norm(a - b)
        if dist < radius:
            yield True
        else:
            yield False






'''
As input we now need the site-wise number of synonymous and non-syn. substitutions of all residues in the protein. We should not write a routine but assume this has been come by on another manner.

theory:

- https://academic.oup.com/mbe/article/32/4/1097/1077799
- https://pubmed.ncbi.nlm.nih.gov/16239014/
- https://www.hyphy.org/resources/slides-selection-2016.pdf
- https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000304
- BUSTED gene-wide selection https://academic.oup.com/mbe/article/32/5/1365/1134918?login=false

tutorials: 

- https://pubmed.ncbi.nlm.nih.gov/25388108/

tools:

- https://www.datamonkey.org/
- https://github.com/rjorton/vnvs
- https://github.com/im3sanger/dndscv

'''

