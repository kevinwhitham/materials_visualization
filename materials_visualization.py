import nglview
from os.path import splitext
import pickle
from ase.io.vasp import read_vasp_out
import pandas as pd
import matplotlib.pyplot as plt
from datetime import timedelta
from ase.io import read
from ase.io import Trajectory
import math
from ipywidgets import HBox, VBox, Label
import numpy as np
import ase.units


def show_ngl_row(mols, show_indices=False, captions=None, trajectories=False, view_axis='y', show_cell=True):
    mols = make_list(mols)

    full_width = 1500
    w = full_width // len(mols)
    if trajectories:
        views = [nglview.show_asetraj(mol) for mol in mols]
    else:
        views = [nglview.show_ase(mol) for mol in mols]

    for view in views:

        # Add indices to all atoms
        if show_indices:
            view.add_label(labelType='atomindex', color='black')

        # Set width of each view
        view._remote_call('setSize', target='Widget', args=[f'{w}px', f'400px'])

        # The default view axis is z
        if view_axis == 'x':
            view.control.spin([1, 0, 0], math.pi / 2)
            view.control.spin([0, 1, 0], math.pi / 2)
        elif view_axis == 'z':
            continue
        else:
            # view along the y (normal to xz plane) by default
            view.control.spin([1, 0, 0], math.pi / 2)

        # Add axes
        if show_cell:
            view.add_representation(repr_type='unitcell')

        # Set view to orthographic
        view.camera = 'orthographic'

    result = HBox(views)

    if captions:
        if len(captions) == len(mols):
            result = HBox([VBox([v, Label(c)]) for v, c in zip(views, captions)])
            for view in views:
                view.center()

    return result



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
        logfile_contents = pd.read_table(timing_file, skiprows=[0, 1], sep=r'\s+', error_bad_lines=False)

        algo_name = logfile_contents.iloc[-1, 0]
        # print('Algorithm name:', algo_name)

        # Default formatting options
        header = 0
        cols = [1, 2, 3, 4]
        skiprows = None
        col_names = ['Step', 'Time', 'Energy', 'fmax']

        # Choose formatting based on algorithm name
        if 'bfgslinesearch' in str(algo_name).lower():
            cols = [1, 3, 4, 5]
            col_names = ['Step', 'FC', 'Time', 'Energy', 'fmax']
        elif 'precon' in str(algo_name).lower():
            header = None

        timing_data = pd.read_table(timing_file, header=header, index_col=1, names=col_names, parse_dates=['Time'],
                                    infer_datetime_format=True, sep=r'\[*\s+|\]\s+', engine='python',
                                    error_bad_lines=False)

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
                if (timing_data.index[i] - timing_data.index[i - 1]) < 1:
                    dt.iloc[i] = timedelta(0)

        # Plot time in units of hours
        plt.plot([td.days * 24 + td.seconds / 3600 for td in dt.cumsum()], timing_data['fmax'], '-', label=label)
        timing_summary = timing_summary.append(
            pd.DataFrame(data=[[timing_file, len(timing_data)]], columns=summary_col_names), ignore_index=True)

    display(timing_summary)
    plt.tight_layout()
    plt.yscale('log')
    plt.xlabel('Time (hours)')
    plt.ylabel('fmax (eV/Ang)')
    # Put a legend to the right of the current axis
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    return fig


def plot_total_displacement(trajectory_filenames, labels):
    trajectory_filenames = make_list(trajectory_filenames)
    labels = make_list(labels)

    fig = plt.figure()

    for file, label in zip(trajectory_filenames, labels):
        traj = Trajectory(file)
        disp = []
        for atoms in traj:
            disp.append(np.sum(np.sqrt(np.sum((traj[0].positions - atoms.positions) ** 2, axis=1))))

        plt.plot(disp, label=label)

    plt.xlabel('Iteration')
    plt.ylabel('Total Displacement (Angstrom)')
    plt.legend()
    return fig


def plot_unit_cell_volume_change(trajectories, labels):
    '''
    Plot relative change in unit cell volume.
    :param trajectories: .traj filenames or Trajectory objects
    :param labels: strings for plot legend
    :return: matplotlib figure
    '''

    trajectories = make_list(trajectories)
    labels = make_list(labels)

    fig = plt.gcf()
    if fig is None:
        fig = plt.figure()

    for traj, label in zip(trajectories, labels):
        if type(traj) == str:
            traj = Trajectory(traj)

        plt.plot(list(range(len(traj))),
                 [round((atoms.get_volume() - traj[0].get_volume()) / traj[0].get_volume() * 100.0, 2) for atoms in
                  traj], label=label + ': $\Delta$V=' + str(
                round((traj[-1].get_volume() - traj[0].get_volume()) / traj[0].get_volume() * 100.0, 2)) + '%')

    plt.plot(plt.xlim(), [0, 0], '--', color='0.5')
    plt.xlabel('Step')
    plt.ylabel('% Volume Change')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    return fig

def get_octahedral_angles_and_distances(center_atom_symbol, vertex_atom_symbol, trajectory, plane_of_interest='xy'):
    '''
    Plots the angle connecting octahedral centers over a trajectory.
    :param center_atom_symbol: name of atom at octahedral centers (e.g. 'Pb')
    :type center_atom_symbol: str
    :param vertex_atom_symbol: name of atom at octahedral vertices (e.g. 'I')
    :type vertex_atom_symbol: str
    :param trajectory: ASE trajectory
    :param plane_of_interest: way to specify equitorial plane (e.g. 'xy')
    :type plane_of_interest: str
    :return: tuple of DataFrame (angles, distances)
    :rtype: tuple
    '''

    if 'x' not in plane_of_interest:
        apical_axis = 0
    elif 'y' not in plane_of_interest:
        apical_axis = 1
    else:
        apical_axis = 2

    # DataFrame to hold the angle data
    angle_data = pd.DataFrame()
    distance_data = pd.DataFrame()

    # Step through the trajectory
    for step,atoms in enumerate(trajectory):
        # For each center atom, find the nearest 6 atoms of vertex type
        all_center_atom_indices = np.array([a.index for a in atoms if a.symbol == center_atom_symbol])
        all_vertex_atom_indices = np.array([a.index for a in atoms if a.symbol == vertex_atom_symbol])
        all_distances = atoms.get_all_distances(mic=True)
        second_center_atom_indices = all_center_atom_indices

        for center_atom_index in all_center_atom_indices:

            # Remove the center atom index to avoid double counting
            second_center_atom_indices = np.delete(second_center_atom_indices,
                                                    np.argwhere(second_center_atom_indices == center_atom_index))

            vertex_atom_indices = all_vertex_atom_indices[np.argsort(all_distances[center_atom_index][all_vertex_atom_indices])][:6]
            #print('vertex_atom_indices', vertex_atom_indices)
            # Choose four vertices closest to the center in the apical direction
            equitorial_vertex_atom_indices = vertex_atom_indices[np.argsort([np.abs(atoms.get_distance(center_atom_index, vertex_index, vector=True, mic=True)[apical_axis]) for vertex_index in vertex_atom_indices])][:4]
            #print('equitorial_vertex_atom_indices', equitorial_vertex_atom_indices)
            for vertex_atom_index in equitorial_vertex_atom_indices:
                if len(second_center_atom_indices):
                    # Get nearest atom of type center_atom_symbol that is not center_atom_index
                    #print('center_atom_index', center_atom_index)
                    #print('all_center_atom_indices', all_center_atom_indices)
                    #print('vertex_atom_index', vertex_atom_index)
                    distance_sorted_center_atom_indices = all_center_atom_indices[np.argsort(all_distances[vertex_atom_index][all_center_atom_indices])]
                    distance_sorted_center_atom_indices = np.delete(distance_sorted_center_atom_indices, np.argwhere(distance_sorted_center_atom_indices == center_atom_index))
                    #print('distance_sorted_center_atom_indices', distance_sorted_center_atom_indices)
                    nearest_center_atom_index = distance_sorted_center_atom_indices[0]
                    if any(nearest_center_atom_index == second_center_atom_indices):
                        angle_data = angle_data.append(pd.DataFrame(dict(step=step,
                                                            angle=atoms.get_angle(center_atom_index,
                                                                                  vertex_atom_index,
                                                                                  nearest_center_atom_index, mic=True),
                                                            atoms=','.join(map(str,[center_atom_index,
                                                                            vertex_atom_index,
                                                                            nearest_center_atom_index]))),
                                                       index=[0]),
                                          ignore_index=True)

                distance_data = distance_data.append(pd.DataFrame(dict(step=step,
                                                                       distance=all_distances[center_atom_index][vertex_atom_index],
                                                                       atoms=','.join(map(str,[center_atom_index, vertex_atom_index]))),
                                                                  index=[0]),
                                                     ignore_index=True)




    return angle_data, distance_data



def vasp_to_trajectory(outcar_filenames, trajectory_filename):
    '''
    Convert and concatenate VASP OUTCAR files to ASE Trajectory

    :param [list,str] outcar_filenames:  paths to OUTCAR files to convert
    :param str trajectory_filename: path to save trajectory file
    :return Trajectory:
    '''

    outcar_filenames = make_list(outcar_filenames)

    atoms_list = []
    for f in outcar_filenames:
        before_len = len(atoms_list)
        atoms_list += read_vasp_out(f, index=':')

    traj = Trajectory(trajectory_filename, mode='w')
    for atoms in atoms_list:
        traj.write(atoms)

    return Trajectory(trajectory_filename)


def make_list(obj):
    if type(obj) is not list:
        obj = [obj]
    return obj

def load_bands(filename):
    '''
    Load a 2D numpy array file of eigenvalues vs k-points.
    :param filename: path to npy file
    :type filename: basestring
    :return: 2D numpy array
    :rtype: array
    '''

    e_mk = np.load(filename)
    emax_n = np.max(e_mk, axis=1) # greatest eigenvalue per band
    soc_vb_n = max(np.argwhere(emax_n < 0))[0]  # Fermi level should be at zero energy

    print('Bands:', e_mk.shape[0])
    print('K-points:', e_mk.shape[1])
    print('Valence band index: ',soc_vb_n)

    return e_mk


def get_band_orbital_weights(bs_calc, species, n, orbital, f_kmsi=None):
    '''
    Get the atomic orbital character for every k-point and band.

    :param bs_calc: band structure calculator
    :param species: string of atomic species e.g. 'Pb'
    :param n: band index
    :param orbital: string of orbital e.g. 's', 'p'
    :param f_kmsi: for spin-orbit bands, provide the projections
    :returns array (k-points x bands)
    '''

    if f_kmsi is not None:
        f_kni = abs(f_kmsi[:, :, 0, :]) + abs(f_kmsi[:, :, 1, :])
    else:
        # projectors method works for LCAO calculations
        f_kni = bs_calc.get_projections(locfun='projectors')

    wfs = bs_calc.wfs

    anl_ki = []

    for kpt in wfs.kpt_u:
        if kpt.s == 0:
            anl_i = []
            for a, P_ni in kpt.P_ani.items():
                i = 0
                setup = wfs.setups[a]
                for lj, nj in zip(setup.l_j, setup.n_j):
                    if nj >= 0:
                        for j in range(i, i + 2 * lj + 1):
                            anl_i.append([a, nj, lj])
                    i += 2 * lj + 1

            anl_ki.append(anl_i)
    anl_ki = np.array(anl_ki)

    letter_to_angular = dict(s=0, p=1, d=2, f=3)
    l = letter_to_angular[orbital]

    w_kn = np.zeros(f_kni.shape[:2])

    for k in range(f_kni.shape[0]):
        for a in [a.index for a in bs_calc.atoms if a.symbol == species]:
            # get a weight from [0,1] for the contribution of a,n,l for all bands at this k point
            anl_index = np.argwhere(np.all(anl_ki[k] == [a, n, l], axis=1)).flatten()
            w_kn[k] += (np.sum((abs(f_kni[k, :, anl_index]) ** 2).T, axis=1) / np.sum(abs(f_kni[k]) ** 2,
                                                                                      axis=1)).flatten()

    return w_kn.T



def plot_bands(e_mk, path_data,
               energy_limits,
               bands_to_highlight=None, band_labels=None,
               title=None,
               weight_nk=None,
               weight_color=(1,0,0),
               weight_label=None,
               thickness=None):
    '''
    Plot a band structure diagram from 2D array of E vs k
    :param e_mk: 2D array of eigenvalues vs k-points
    :type e_mk: numpy array
    :param energy_limits: list of min,max energies to plot
    :type energy_limits: list
    :param path_data: band path and k-points tuple of (x, X, labels)
    :type path_data: tuple
    :param bands_to_highlight: list of band indices e.g. [660, 662]
    :type bands_to_highlight: list
    :param band_labels: list of strings to show in legend e.g. ['Valence', 'Conduction']
    :type band_labels: list
    :param weight_nk: array with same shape as e_mk giving a weight for each point
    :type weight_nk: numpy array
    :param weight_color: RGB value to color points of high weight
    :type weight_color: tuple
    :param thickness: optional size of symbols
    :type thickness: int
    :return: None
    :rtype: None
    '''

    band_max = np.max(e_mk, axis=1)
    if bands_to_highlight is None:
        valence_band_index = np.max(np.argwhere(band_max < 0))
        bands_to_highlight = [valence_band_index, valence_band_index + 1]
        band_labels = ['Valence', 'Conduction']

    if band_labels is None:
        band_labels = bands_to_highlight

    bands_to_highlight = make_list(bands_to_highlight)
    band_labels = make_list(band_labels)

    min_plot_energy = min(energy_limits)
    max_plot_energy = max(energy_limits)

    for b in bands_to_highlight:
        band_min = np.min(e_mk[b])
        band_max = np.max(e_mk[b])
        print(
            f'Width of band {b}: {np.round(band_max - band_min, 4)} ({np.round(band_min, 4)} to {np.round(band_max, 4)})')

    def pretty_label(label):
        if label == 'G':
            label = r'$\Gamma$'
        elif len(label) > 1:
            # Assume extra chars are part of a subscript, e.g. M1 becomes $M_{1}$
            label = '$' + label[0] + '_{' + str(label[1:]) + '}$'
        return label

    # Get band path
    # x are the bandpath points in k-space
    # X are the symmetry point locations in k-space

    # for backward compatability, check if path_data is a path to a pyc file
    if type(path_data) is str:
        with open(path_data, 'rb') as file:
            x, X, orig_labels = pickle.load(file)
    else:
        x, X, orig_labels = path_data

    labels = [pretty_label(l) for l in orig_labels]

    # Band structure diagrams
    if plt.gca() is None:
        plt.figure(figsize=(4, 3), dpi=128)

    plt.xticks(X, labels, size=10)
    plt.yticks(size=10)

    # Plot vertical grey lines at each symmetry point label
    for i in range(len(X))[1:-1]:
        plt.plot(2 * [X[i]], [min_plot_energy, max_plot_energy],
                 c='0.5', linewidth=0.5)

    # Plot the spin-orbit corrected bands in grey
    # If weights are given, use them to color the data
    from matplotlib.colors import LinearSegmentedColormap
    weight_cmap = LinearSegmentedColormap.from_list('weight_cmap', [(1,1,1),weight_color], N=256)

    band_max = np.max(e_mk, axis=1)
    band_min = np.min(e_mk, axis=1)
    bands_to_plot = np.argwhere((band_max > min_plot_energy) & (band_min < max_plot_energy))
    for band in bands_to_plot:
        if weight_nk is None:
            if thickness is None:
                thickness = 0.5

            plt.plot(x, e_mk[band].flatten(), c='0.5', linewidth=thickness)

        else:
            if thickness is None:
                thickness = 200

            # Vary color and thickness by weight
            # plt.scatter(x, e_mk[band], c=weight_nk[band], cmap=weight_cmap, vmin=0, vmax=1, marker='.',
            #             s=thickness * weight_nk[band], alpha=0.5, edgecolors='none')

            # Vary just thickness by weight
            plt.scatter(x, e_mk[band], color=weight_color, marker='.',
                        s=thickness * weight_nk[band], alpha=0.5, edgecolors='none')

    # for the legend
    if weight_nk is not None:
        plt.scatter(x[0], e_mk[0,0], color=weight_color, label=weight_label)

    # Plot the bands of interest in colors
    # Plot in descending order so that the legend shows higher energy bands at the top
    bands_of_interest = np.array(bands_to_highlight)
    band_labels = np.array(band_labels)
    band_order = list(np.argsort(bands_of_interest)[::-1])
    for boi, label in zip(bands_of_interest[band_order], band_labels[band_order]):
        plt.plot(x, e_mk[boi], lw=1, label=label)

    # Plot a horizontal dotted grey line at zero energy
    plt.plot([0.0, x[-1]], 2 * [0.0], c='0.5', linestyle=':')
    plt.ylabel(r'$\varepsilon_n(k)$ [eV]', size=10)
    plt.axis([0, x[-1], min_plot_energy, max_plot_energy])

    if len(bands_to_highlight):
        plt.legend()

    if title:
        plt.title(title)
    else:
        plt.title(f'Band {bands_to_highlight}')

    plt.tight_layout()



def plot_band_path(structure_file, band_path_str):
    '''
    Plot a path in the Brillouin zone as a png file and output a text file
    with the special k-points.
    :param structure_file: path to structure file readable by ASE
    :type structure_file: str
    :param path: special points, e.g. 'XGY'
    :type path: str
    :return: None
    :rtype:
    '''

    atoms = read(structure_file)

    basename, _ = splitext(structure_file)

    lat = atoms.cell.get_bravais_lattice()

    bp = atoms.cell.bandpath(band_path_str, 48)
    reduced_bp = lat.bandpath(band_path_str, 48)

    plt.figure(figsize=(8, 8), dpi=128)
    bp.plot()
    plt.savefig(f'{basename}_band_path.png')

    with open(f'{basename}_band_path.log', 'w') as file:
        file.write('Reduced Bravais Lattice:\n')
        file.write(lat.description())
        file.write(f'Path: {band_path_str}\n')
        file.write('K-points for reduced bravais lattice:\n')
        for p in list(band_path_str):
            file.write(f'{p}: {reduced_bp.special_points[p]}\n')
        file.write('K-points for structure as input:\n')
        file.write(f'a={atoms.cell.cellpar()[0]:.4f}, b={atoms.cell.cellpar()[1]:.4f}, c={atoms.cell.cellpar()[2]:.4f}\n')
        for p in list(band_path_str):
            file.write(f'{p}: {bp.special_points[p]}\n')