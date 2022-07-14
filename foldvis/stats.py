from hdbscan import HDBSCAN
import numpy as np

from foldvis.geometry import get_alpha_carbon_atoms 


def spatial_association(fold, features, statistic='Gi', radius=8, coordinates='alpha_carbons'):
    '''
    Calculates statistic(s) that can be used to test for spatial
    association.

    TODO:

    - multiple comparison correction (control FDR)
    '''
    assert len(fold) == len(features), 'Residue to feature mapping not unique'

    l = []
    for pos in range(len(fold)):
        sphere = list(is_close(pos, fold, radius, coordinates))
    
        fn = globals()[statistic]
        # https://stackoverflow.com/questions/3061/calling-a-function-of-a-module-by-using-its-name-a-string

        _, z = fn(sphere, features, pos)
        l.append(z)
    
    return l


def Gi(w_ij, features, pos):
    '''
    Getis-Ord statistic for spatial association.

    "The analysis of Spatial Association by Use of Distance Statistics",
    Getis & Ord, Geographical Analysis, 1992

    Gi, other than Gi_star (see below) will exclude the position in the
    "center" of the current fn call (i != j).
    '''
    n = len(w_ij)

    # exclude the query position
    w_ij = w_ij[:pos] + [False] + w_ij[pos+1:]
    features = features[:pos] + [0] + features[pos+1:]

    scores = sum([score for include, score in zip(w_ij, features) if include])
    total = sum(features)
    G = scores / total
    
    Wi = sum(w_ij)
    Yi1 = total / (n-1)
    
    x = sum(i**2 for i in features)
    Yi2 = (x / (n-1)) - (Yi1**2)
    
    E = Wi / (n-1)
    Var = (Wi * (n-1-Wi) * Yi2) / ((n-1)**2 * (n-2) * Yi1**2)
    Z = (G-E) / np.sqrt(Var)
    
    return G, Z


def Gi_star(w_ij, features, pos):
    '''
    "The analysis of Spatial Association by Use of Distance Statistics",
    Getis & Ord, Geographical Analysis, 1992

    Gi_star, other than Gi (see above) will include the position in the
    "center" of the current fn call (i ==j).
    '''
    n = len(w_ij)

    scores = sum([score for include, score in zip(w_ij, features) if include])
    total = sum(features)
    G = scores / total

    Wi = sum(w_ij)
    Yi1 = total / n

    x = sum((features[pos] * i)**2 for i in features)
    Yi2 = (x / n) - (Yi1**2)
    # The double sum in table 1 of the original paper just means add x_i to
    # the sum of all the other features.

    E = Wi / n
    Var = (Wi * (n-Wi) * Yi2) / (n**2 * (n-1) * Yi1**2)
    Z = (G-E) / np.sqrt(Var)
    
    return G, Z


def cluster(fold, mask, *args, **kwargs):
    '''
    from foldvis.stats import cluster
    from foldvis.geometry import get_alpha_carbon_atoms

    mask = [1 if i < 0.05 else 0 for i in d['meme']['positive']['scores']]
    cluster(model, mask, min_cluster_size=2)
    '''
    points = list(get_alpha_carbon_atoms(fold, only_coords=True))
    X = [i for i, j in zip(points, mask) if j]
    clusterer = HDBSCAN(*args, **kwargs)
    return clusterer.fit_predict(X)



'''
TODO:

random forest, predict clusters from features, then variable importance to see
which features matter most.

==?

enrichment: given spatial clusters, are they enriched in any feature?

- distance to ligand
- solvent accessibility
- interface
'''


'''
TODO: We could feed the clusters to the Getis-Ord statistic or calculate eg
the (adjusted) Rand score.
'''

