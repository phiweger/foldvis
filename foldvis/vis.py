import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
# plt.get_cmap('viridis')

import py3Dmol


def plot_alphafold(af, width=400, height=300):

    styles = {
        'A': ('black', 1),
        'B': ('grey', 0.50),
        'C': ('grey', 0.50),
        'D': ('grey', 0.50),
        'E': ('grey', 0.50),
    }
    view=py3Dmol.view(width=width, height=height)
    for m in af.models.values():
        view.addModel(m.to_stream(),'pdb')

    view.setBackgroundColor('white')

    for k, (color, opacity) in styles.items():
        view.setStyle({'chain': k},{'cartoon': {'color': color, 'opacity': opacity}})

    view.zoomTo()
    return view


def map_colors(v, palette='viridis'):
    norm = matplotlib.colors.Normalize(vmin=min(v), vmax=max(v), clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=plt.get_cmap(palette))
    
    cols = [matplotlib.colors.to_hex(mapper.to_rgba(i)) for i in v]
    return cols


def plot_annotation(model, label, palette='viridis', surface=False, opacity=1., width=400, height=300):
    '''
    Available color maps:

    https://matplotlib.org/stable/tutorials/colors/colormaps.html
    '''
    assert label in model.annotation.keys(), 'Label not found in annotation'
    stream = model.to_stream()
    
    # https://stackoverflow.com/questions/28752727/map-values-to-colors-in-matplotlib
    anno = model.annotation[label]
    mn = min(anno)
    mx = max(anno)
    
    norm = matplotlib.colors.Normalize(vmin=mn, vmax=mx, clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=plt.get_cmap(palette))
    
    # Map colors to residues
    d = {k+1: matplotlib.colors.to_hex(mapper.to_rgba(v)) for k, v in enumerate(anno)}
    
    view = py3Dmol.view(width=width, height=height)
    # view.addModelsAsFrames(stream)
    view.addModel(stream, 'pdb')

    i = 0
    l = []
    for line in stream.split("\n"):
        split = line.split()
        if len(split) == 0 or split[0] != "ATOM":
            continue

        color = d[int(split[5])]
        l.append(color)
        idx = int(split[1])
        view.setStyle({'model': -1, 'serial': i+1}, {'cartoon': {'color': color}})
        i += 1
        
    if surface:
        map_ = {(i+1): j for i, j in zip(range(len(model)), map_colors(anno, palette))}
        view.addSurface(py3Dmol.VDW, {'opacity': opacity, 'colorscheme': {'prop': 'resi', 'map': map_}})

    view.zoomTo()
    return view


def plot_superposition(models, colors, width=400, height=300):
    
    view=py3Dmol.view(width=width, height=height)
    for m in models:
        view.addModel(m.to_stream(),'pdb')

    view.setBackgroundColor('white')

    for chain, color in colors.items():
        view.setStyle({'chain': chain}, {'cartoon': {'color': color}})

    view.zoomTo()
    return view

