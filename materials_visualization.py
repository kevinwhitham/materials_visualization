import nglview
import math
from ipywidgets import HBox, VBox, Label

def show_ngl_row(mols, show_indices=False, captions=None, trajectories=False):
    full_width = 1500
    w = full_width//len(mols)
    if trajectories:
        views = [nglview.show_asetraj(mol) for mol in mols]
    else:
        views = [nglview.show_ase(mol) for mol in mols]

    for view in views:

        # Add indices to all atoms
        if show_indices:
            view.add_label(labelType='atomindex', color='black')

        # Set width of each view
        view._remote_call('setSize', target='Widget', args=[f'{w}px',f'400px'])

        # We usually want to view the xz plane of a unit cell
        view.control.spin([1, 0, 0], math.pi/2)

        # Add axes
        view.add_representation(repr_type='unitcell')

        # Set view to orthographic
        view.camera = 'orthographic'
    
    result = HBox(views)
    
    if captions:
        if len(captions) == len(mols):
            result = HBox([VBox([v,Label(c)]) for v,c in zip(views,captions)])
            for view in views:
                view.center()
    
    return result