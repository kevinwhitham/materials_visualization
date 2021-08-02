import nglview
import math
from ipywidgets import HBox, VBox, Label

def show_ngl_row(mols, show_indices=False, captions=None, trajectories=False, view_axis='y'):

    mols = make_list(mols)

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

        # The default view axis is z
        if view_axis == 'x':
            view.control.spin([1, 0, 0], math.pi/2)
            view.control.spin([0, 1, 0], math.pi/2)
        elif view_axis == 'z':
            continue
        else:
            # view along the y (normal to xz plane) by default
            view.control.spin([1, 0, 0], math.pi / 2)

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

import pandas as pd
import matplotlib.pyplot as plt
from datetime import timedelta

def plot_fmax_vs_time(timing_filenames, labels=None):
    '''
    Plot data from an optimization log file.

    :param timing_filenames: [List] paths to log files
    :param labels: [List] optional labels for plot legend, otherwise will use filenames
    :return: None
    '''

    timing_filenames = make_list(timing_filenames)
    if labels:
        labels = make_list(labels)
    
    # Create a new figure
    fig = plt.figure()

    summary_col_names = ['File', 'iterations']
    timing_summary = pd.DataFrame(data=None, columns=summary_col_names)
    for i, timing_file in enumerate(timing_filenames):
      
        # if labels are supplied, use them in the plot legend
        label = None
        if labels != None:
            if len(labels) == len(timing_filenames):
                label = labels[i]

        # otherwise use filenames in the plot legend
        if label == None:
            label = timing_file

        # Get the type of optimization so we can interpret the log data
        # Sometimes there is no header, sometimes 1 line, sometimes 2 lines, so we skip the first 2
        # Hopefully there is more than two rows of data in the file
        logfile_contents = pd.read_table(timing_file, skiprows=[0,1], sep=r'\s+', error_bad_lines=False)

        algo_name = logfile_contents.iloc[-1,0]
        #print('Algorithm name:', algo_name)

        # Default formatting options
        header=0
        cols = [1,2,3,4]
        skiprows = None
        col_names = ['Step', 'Time', 'Energy', 'fmax']

        # Choose formatting based on algorithm name
        if 'bfgslinesearch' in str(algo_name).lower():
            cols = [1,3,4,5]
            col_names = ['Step', 'FC', 'Time', 'Energy', 'fmax']
        elif 'precon' in str(algo_name).lower():
            header=None

        timing_data = pd.read_table(timing_file, header=header, index_col=1, names=col_names, parse_dates=['Time'], infer_datetime_format=True, sep=r'\[*\s+|\]\s+', engine='python', error_bad_lines=False)

        # Correct for change of day in elapsed time
        dt = timing_data['Time'].diff()
        dt[0] = timedelta(0)

        for i, d in enumerate(dt):
            if d < timedelta(0):
                dt.iloc[i] += timedelta(days=1)
            # if logfiles have been concatenated, the step numbers may not be consecutive
            # also there will be a gap in the time between runs
            # We do not know how long the first iteration of a run takes
            # we could set time of the first step equal to the time of the previous or next step
            # or just set the time difference to zero
            if i > 0:
                if (timing_data.index[i] - timing_data.index[i-1]) < 1:
                    dt.iloc[i] = timedelta(0)
        
        # Plot time in units of hours
        plt.plot([td.days*24+td.seconds/3600 for td in dt.cumsum()], timing_data['fmax'], '-o', label=label)
        timing_summary = timing_summary.append(pd.DataFrame(data=[[timing_file,len(timing_data)]], columns=summary_col_names), ignore_index=True)

    display(timing_summary)  
    plt.tight_layout()
    plt.yscale('log')
    plt.xlabel('Time (hours)')
    plt.ylabel('fmax (eV/Ang)')
    # Put a legend to the right of the current axis
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    return fig

from ase.io import read
from ase.io import Trajectory
import numpy as np
import matplotlib.pyplot as plt

def plot_total_displacement(starting_structure_filename, trajectory_filenames, labels):
    starting_structure = read(starting_structure_filename)

    trajectory_filenames = make_list(trajectory_filenames)
    labels = make_list(labels)

    for file, label in zip(trajectory_filenames,labels):
            traj = Trajectory(file)
            disp = []
            for atoms in traj:
                disp.append(np.sum(np.sqrt(np.sum((starting_structure.positions - atoms.positions)**2, axis=1))))

            plt.plot(disp, label=label)

    plt.xlabel('Iteration')        
    plt.ylabel('Total Displacement (Angstrom)')
    plt.legend()        
    plt.show()


def plot_unit_cell_volume_change(trajectory_filenames, labels):

    trajectory_filenames = make_list(trajectory_filenames)
    labels = make_list(labels)

    fig = plt.figure()

    for filename, label in zip(trajectory_filenames, labels):
        traj = Trajectory(filename)
        plt.plot(list(range(len(traj))),
                 [round((atoms.get_volume() - traj[0].get_volume()) / traj[0].get_volume() * 100.0, 2) for atoms in
                  traj], label=label + ': $\Delta$V=' + str(
                round((traj[-1].get_volume() - traj[0].get_volume()) / traj[0].get_volume() * 100.0, 2)) + '%')

    plt.plot(plt.xlim(), [0, 0], '--', color='0.5')
    plt.xlabel('Step')
    plt.ylabel('% Volume Change')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    return fig

def make_list(obj):
    if type(obj) is not list:
        obj = [obj]
    return obj