import nglview
from ipywidgets import HBox, VBox, Label

def show_ngl_row(mols, show_indices=True, captions=None):
    full_width = 1500
    w = full_width//len(mols)
    views = [nglview.show_ase(mol) for mol in mols]
    for view in views:
        if show_indices:
            view.add_label(labelType='atomindex', color='black')
        view._remote_call('setSize', target='Widget', args=[f'{w}px',f'400px'])
    
    result = HBox(views)
    
    if captions:
        if len(captions) == len(mols):
            result = HBox([VBox([v,Label(c)]) for v,c in zip(views,captions)])
            for view in views:
                view.center()
    
    return result