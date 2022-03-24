import glob
import re
import nglview
import os.path
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
from gpaw import GPAW

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
            view.control.spin([0, 1, 0], math.pi / 2)
            view.control.spin([1, 0, 0], math.pi / 2)
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

    return result, views



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
    for step, atoms in enumerate(trajectory):

        # Structure must be periodic to find all angles and distances
        atoms.set_pbc(True)

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
        atoms_list += read_vasp_out(f, index=':')

    traj = Trajectory(trajectory_filename, mode='w')
    for atoms in atoms_list:
        traj.write(atoms)

    return Trajectory(trajectory_filename)


def plot_relaxation(traj, label, fmax_target=0.01, incar_files=None):
    '''
    Plot progress of a relaxation (pressure, fmax, energy, volume)
    :param traj: the steps of the relaxation
    :type traj: ASE trajectory
    :param label: name for this relaxation
    :type label: str
    :param fmax_target: attempt to estimate at which step fmax will be achieved
    :type fmax_target: float
    :param incar_files: If this is a VASP calculation, you may plot SMASS and POTIM
    :type incar_files: list of str
    :return: None
    :rtype:
    '''

    if incar_files:
        # Get SMASS and POTIM values
        potim = []
        smass = []
        for file_name in incar_files:
            # Get number of steps in this run
            outcar_file_name = re.sub('INCAR', 'OUTCAR', file_name)
            total_time = 0
            steps = 0
            with open(outcar_file_name) as outcar_file:
                search_result = re.findall(r'LOOP\+:\s+cpu time\s+(\d+)\.', outcar_file.read())

            if search_result:
                for loop_time in search_result:
                    total_time += int(loop_time)
                    steps += 1

            with open(file_name) as incar_file:
                txt = incar_file.read()
                potim_search = re.search(r'POTIM\s+=\s+(\d*\.?\d*)', txt)
                smass_search = re.search(r'SMASS\s+=\s+(\d*\.?\d*)', txt)
                if potim_search:
                    potim_value = float(potim_search.group(1))
                    for step in range(steps):
                        potim.append(potim_value)
                if smass_search:
                    smass_value = float(smass_search.group(1))
                    for step in range(steps):
                        smass.append(smass_value)

    cols = 1
    rows = 4
    if len(potim) and len(smass):
        rows = 5

    fig = plt.gcf()
    if fig is None:
        plt.subplots(rows, cols, figsize=(3.25, rows*3.25))

    if len(potim) and len(smass):
        ax = plt.subplot(rows, cols, 5)
        plt.plot(list(range(len(potim))), potim, 'rx-', label='POTIM')
        plt.ylabel('POTIM', color='r')
        plt.yticks(color='r')
        plt.gca().twinx()
        plt.plot(list(range(len(smass))), smass, 'bo-', label='SMASS')
        plt.ylabel('SMASS', color='b')
        plt.yticks(color='b')

        plt.subplot(rows, cols, 4, sharex=ax)
        plot_unit_cell_volume_change(traj, labels=[label])
        plt.tick_params('x', labelbottom=False, bottom=False)
        plt.xlabel('')

        plt.sca(ax)
        plt.xlabel('Step')
    else:
        ax = plt.subplot(rows, cols, 4)
        plt.subplot(rows, cols, 4)
        plot_unit_cell_volume_change(traj, labels=[label])

    plt.sca(plt.subplot(rows, cols, 1, sharex=ax))
    # Sign convention is opposite for VASP (<0 is tension) vs ASE (<0 is compression)
    # Default pressure units in ASE are eV/Angstrom^3
    pressure_kbar = [np.trace(atoms.get_stress(voigt=False)) / 3 / ase.units.GPa * -10 for atoms in traj]
    print(f'Final pressure: {pressure_kbar[-1]:.4f} kBar')
    plt.plot(list(range(len(traj))), pressure_kbar, 'o-', label=label)
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    plt.ylabel('Pressure (kBar)')
    plt.tick_params('x', labelbottom=False, bottom=False)

    # We can estimate the accuracy of the pressure by comparing repeated calculations on the same structure
    # with PREC = Normal, (EDIFF = 1E-4) the accuracy of the pressure is approximately 0.1 kBar

    plt.subplot(rows, cols, 2, sharex=ax)
    max_forces = [np.max(np.linalg.norm(atoms.get_forces(), axis=1)) for atoms in traj]
    print(f'Minimum fmax:\t{min(max_forces):.4e} eV/Ang.')
    print(f'Final fmax:\t{max_forces[-1]:.4e} ev/Ang.')
    plt.plot(list(range(len(traj))), max_forces, 'o-', label=label)
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))

    # Predict how many steps until fmax < threshold
    if len(max_forces) > 4:
        fmax_fit_start = len(max_forces) - len(max_forces)//4
        fmax_fit_end = len(max_forces)
        m, b = \
        np.linalg.lstsq(np.vstack([np.arange(fmax_fit_start, fmax_fit_end), np.ones(fmax_fit_end - fmax_fit_start)]).T,
                        np.log(max_forces[fmax_fit_start:fmax_fit_end]), rcond=None)[0]
        print(f'fmax will be < {fmax_target} at step {np.ceil((np.log(fmax_target) - b) / m)}')
        plt.plot(np.arange(fmax_fit_start, fmax_fit_end), np.exp(m * np.arange(fmax_fit_start, fmax_fit_end) + b), '-',
                 color='red')

    plt.yscale('log')
    plt.ylabel('fmax (eV/$\AA$)')
    plt.tick_params('x', labelbottom=False, bottom=False)

    plt.subplot(rows, cols, 3, sharex=ax)
    plt.plot(list(range(len(traj))), [a.get_potential_energy() for a in traj], label=label)
    plt.ylabel('Energy (eV)')
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    plt.tick_params('x', labelbottom=False, bottom=False)



def plot_trajectory_angles_and_distances(traj, atom1, atom2, label):
    '''
    Plot angles and distances between two atom types.

    :param traj: sequence of structures
    :type traj: ASE trajectory
    :param atom1: symbol of the center atom
    :type atom1: str
    :param atom2: symbol of the vertex atom
    :type atom2: str
    :param label: label for plots
    :type label: str
    :return: None
    :rtype:
    '''

    # Plot Pb-I-Pb angles and distances
    angle_data, distance_data = get_octahedral_angles_and_distances(atom1, atom2, traj)

    fig = plt.gcf()
    if fig is None:
        fig, axes = plt.subplots(2,1)

    ax = plt.subplot(2,1,1)
    angle_data.pivot(index='step', columns='atoms', values='angle').plot(ax=ax)
    plt.title(label)
    #plt.legend(title='Angle No.', loc='center left', bbox_to_anchor=(1, 0.5))
    plt.gca().get_legend().remove()
    plt.ylabel(f'{atom1}-{atom2}-{atom1} Angle (deg)')
    plt.tick_params('x', labelbottom=False, bottom=False)
    print(f'Final {atom1}-{atom2}-{atom1} Angle', angle_data.query(f'step=={angle_data["step"].max()}')['angle'].mean(),
          '+/-',
          angle_data.query(f'step=={angle_data["step"].max()}')['angle'].std())

    ax = plt.subplot(2,1,2, sharex=ax)
    distance_data.pivot(index='step', columns='atoms', values='distance').plot(ax=ax)
    #plt.legend(title='Bond No.', loc='center left', bbox_to_anchor=(1, 0.5))
    plt.gca().get_legend().remove()
    plt.xlabel('Step')
    plt.ylabel(f'{atom1}-{atom2} Distance ($\AA$)')

    print(f'Final {atom1}-{atom2} distance',
          distance_data.query(f'step=={distance_data["step"].max()}')['distance'].mean(), '+/-',
          distance_data.query(f'step=={distance_data["step"].max()}')['distance'].std())

def plot_trajectory_structure_params(traj):
    '''
    Plot lattice vector lengths and angles over a trajectory
    :param traj: trajectory
    :type traj: ASE Trajectory
    :return:
    :rtype:
    '''

    fig = plt.gcf()
    if fig is None:
        fig = plt.figure()

    # Get structure lengths and angles
    structure = np.empty((len(traj), 6))
    for step, atoms in enumerate(traj):
        structure[step] = atoms.cell.cellpar()
        structure[step] = (structure[step] - traj[0].cell.cellpar())/traj[0].cell.cellpar()

    labels = ['a', 'b', 'c', 'alpha', 'beta', 'gamma']
    for trace in range(6):
        plt.plot(structure[:,trace], label=labels[trace])

    plt.legend(loc='center left', bbox_to_anchor=(1,0.5))
    plt.ylabel('$\Delta$ (%)')
    plt.xlabel('Step')

def plot_vasp_relaxations(exclude_keywords=[], convergence_steps=10, fmax_target=0.01):
    '''
    Finds all INCAR, OUTCAR files in all sub directories and plots relaxations.
    :param exclude_keywords: exclude paths that include these string
    :type exclude_keywords: list of str
    :param convergence_steps: number of steps to calculate standard deviation in structural parameters
    :type convergence_steps: int
    :param fmax_target: estimate the number of steps to reach this fmax
    :type fmax_target: float
    :return: volume change, fmax, functional, Pb-I-Pb angle, unit cell vector lengths
    :rtype: DataFrame
    '''
    import glob
    import os.path
    import materials_visualization as mv
    import pandas as pd
    from ase.io import Trajectory
    import matplotlib.pyplot as plt
    import numpy as np

    plt.rcdefaults()

    # get directories, assume all directories correspond to a relaxation
    dirs = glob.glob('**/', recursive=True)

    data = []
    paths = []
    relaxation_summary = pd.DataFrame()
    for i, path in enumerate(dirs):

        # Exclude some data sets
        if not any([word in path for word in exclude_keywords]):
            # Get a list of OUTCAR files
            files = glob.glob(path + 'OUTCAR_*')

            if len(files):
                files.sort(key=lambda x: int(re.search(r'OUTCAR_(\d+)',x).group(1)))
                data.append(files)
                paths.append(path)

    print('Plotting data from:', paths)
    print('OUTCAR files:')
    for files in data:
        base_path = os.path.commonpath(files)
        print(f'{base_path}:\t' + ', '.join([os.path.basename(f) for f in files]))

    for i, files in enumerate(data):
        path = paths[i]
        label = path[:-1]  # remove trailing /


        incar_files = [re.sub(r'OUTCAR', 'INCAR', filename) for filename in files]
        print('\n'+incar_files[-1]+':')
        with open(incar_files[-1]) as file:
            incar = file.read()
            print(incar)

        # Concatenate sorted OUTCAR files into one trajectory
        traj = mv.vasp_to_trajectory(files, path + 'vasp_relaxation.traj')

        angle_data, distance_data = mv.get_octahedral_angles_and_distances('Pb', 'I', traj[-1:])
        pb_i_pb_angle = np.mean(angle_data['angle'])
        relaxation_summary = relaxation_summary.append(pd.DataFrame(
            dict(delta_volume_pct=(traj[-1].cell.volume - traj[0].cell.volume) / traj[0].cell.volume * 100.0,
                 fmax_final=np.max(np.linalg.norm(traj[-1].get_forces(), axis=1)),
                 functional=label,
                 pb_i_pb_angle=pb_i_pb_angle,
                 a_vector_delta_pct=(traj[-1].cell.cellpar()[0] - traj[0].cell.cellpar()[0]) / traj[0].cell.cellpar()[
                     0] * 100.0,
                 b_vector_delta_pct=(traj[-1].cell.cellpar()[1] - traj[0].cell.cellpar()[1]) / traj[0].cell.cellpar()[
                     1] * 100.0,
                 c_vector_delta_pct=(traj[-1].cell.cellpar()[2] - traj[0].cell.cellpar()[2]) / traj[0].cell.cellpar()[
                     2] * 100.0,
                 volume_delta_pct_std=np.std([(atoms.cell.volume - traj[0].cell.volume)/traj[0].cell.volume * 100.0
                                    for atoms in traj[-convergence_steps:]]),
                 a_vector_std=np.std([atoms.cell.cellpar()[0] for atoms in traj[-convergence_steps:]]),
                 b_vector_std=np.std([atoms.cell.cellpar()[1] for atoms in traj[-convergence_steps:]]),
                 c_vector_std=np.std([atoms.cell.cellpar()[2] for atoms in traj[-convergence_steps:]]),
                 ),
            index=[i]
            )
                                                       )
        plt.figure(dpi=128, figsize=(6,6))
        plot_relaxation(traj, label=label, incar_files=incar_files, fmax_target=fmax_target)

        # Show now so that it shows below printed info
        plt.show()

        plt.figure(dpi=128, figsize=(6, 6))
        plot_trajectory_angles_and_distances(traj, 'Pb', 'I', label)
        plt.show()

        plot_trajectory_structure_params(traj)
        plt.show()

    return relaxation_summary

def compare_relaxations(relaxation_summary):
    '''
    Bar plot comparison of change in volume after relaxation by different functionals.
    :param relaxation_summary: contains delta_volume_pct, functional
    :type relaxation_summary: DataFrame
    :return:
    :rtype:
    '''
    plt.style.reload_library()
    plt.rcdefaults()
    plt.style.use(['acs'])

    # Plot relative change in volume
    plt.subplots(1, 2, dpi=300, figsize=(3.25, 3.25 / 2))
    plt.subplot(1, 2, 2)
    plt.gca().grid(axis='y', linewidth=0.5)
    plt.gca().set_axisbelow(True)
    for i in relaxation_summary.index:
        plt.bar(x=i, height=relaxation_summary.iloc[i]['delta_volume_pct'],
                label=relaxation_summary.iloc[i]['functional'],
                yerr=relaxation_summary.iloc[i]['volume_delta_pct_std'])
    plt.axhline(y=0, color='grey', linewidth=0.5)
    plt.ylabel('$\Delta V_0$ (%)')
    plt.xlabel(None)
    plt.xticks([])

    # Plot absolute change in volume
    plt.subplot(1, 2, 1)
    plt.gca().grid(axis='y', linewidth=0.5)
    plt.gca().set_axisbelow(True)
    for i in relaxation_summary.index:
        plt.bar(x=i, height=abs(relaxation_summary.iloc[i]['delta_volume_pct']),
                label=relaxation_summary.iloc[i]['functional'],
                yerr=relaxation_summary.iloc[i]['volume_delta_pct_std'])
    plt.axhline(y=0, color='grey', linewidth=0.5)
    plt.ylabel('$| \Delta V_0 |$ (%)')
    plt.xlabel(None)
    plt.xticks([])
    plt.legend(loc='best', fontsize=4, edgecolor='w', borderpad=0)

    # Add panel labels
    for i, label in enumerate(('a', 'b')):
        ax = plt.subplot(1, 2, i + 1)
        ax.text(-0.1, 1.15, label, transform=ax.transAxes,
                fontsize=6, fontweight='bold', va='top', ha='right')

    plt.tight_layout()
    plt.subplots_adjust(left=0)

def get_vasp_runtimes(exclude_keyword=None):
    '''

    :param exclude_keyword: exclude data in paths with this string
    :type exclude_keyword: str
    :return: run time data
    :rtype: DataFrame
    '''
    # get directories, assume all directories correspond to a relaxation
    dirs = glob.glob('**/', recursive=True)

    data = []
    labels = []
    for i, path in enumerate(dirs):
        # Exclude some data sets
        if (exclude_keyword is None) or (exclude_keyword not in path):
            # Get a list of OUTCAR files
            files = glob.glob(path + 'OUTCAR_*')

            if len(files):
                files.sort(key=lambda x: int(re.search(r'OUTCAR_(\d+)',x).group(1)))
                data.append(files)
                labels.append(path[:-1])


    run_times = []
    iterations = []
    for outcar_files in data:
        total_time = 0
        steps = 0
        for outcar in outcar_files:
            with open(outcar) as txt:
                res = re.findall(r'LOOP\+:\s+cpu time\s+(\d+)\.', txt.read())

            if res:
                for loop_time in res:
                    total_time += int(loop_time)
                    steps += 1

        run_times.append(total_time / 3600)
        iterations.append(steps)

    return pd.DataFrame(dict(functional=labels, run_time_hrs=run_times, iterations=iterations))

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


def get_band_orbital_weights(bs_calc, species, n, orbital, M=None, atoms=None, f_kmsi=None):
    '''
    Get the atomic orbital character for every k-point and band.

    :param bs_calc: band structure calculator
    :param species: string of atomic species e.g. 'Pb'
    :param n: principal quantum number
    :param orbital: string of orbital e.g. 's', 'p'
    :param M: list of total angular momentum to plot. Can be any integers from -L to L (i.e. -1,0,1 for orbital="s")
    :param atoms: limit contribution to specific atom indices
    :param f_kmsi: for spin-orbit bands, provide the projections
    :returns array (k-points x bands)
    '''

    if type(bs_calc) is str:
        bs_calc = GPAW(bs_calc)

    if f_kmsi is not None:
        # add spin up and spin down contributions to each band
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

    if M is None:
        M = np.arange(-l, l + 1)

    if atoms is None:
        atoms = [a.index for a in bs_calc.atoms if a.symbol == species]

    w_kn = np.zeros(f_kni.shape[:2])

    for k in range(f_kni.shape[0]):
        for a in atoms:
            # get a weight from [0,1] for the contribution of a,n,l,m for all bands at this k point
            anl_index = np.argwhere(np.all(anl_ki[k] == [a, n, l], axis=1)).flatten()
            for im in np.argwhere(np.array(M) == np.arange(-l, l + 1)):
                w_kn[k] += (np.sum((abs(f_kni[k, :, anl_index[im]]) ** 2).T, axis=1) / np.sum(abs(f_kni[k]) ** 2,
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

    plt.xticks(X, labels)

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
    plt.ylabel(r'$\varepsilon_n(k)$ [eV]')
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
    :return: figure
    :rtype: Matplotlib figure
    '''

    atoms = read(structure_file)

    basename, _ = os.path.splitext(structure_file)

    lat = atoms.cell.get_bravais_lattice()

    bp = atoms.cell.bandpath(band_path_str, 48)
    reduced_bp = lat.bandpath(band_path_str, 48)

    fig = plt.gcf()
    if fig is None:
        fig = plt.figure(figsize=(8, 8), dpi=128)

    # Increase the size of the special point labels
    import matplotlib as mpl
    with mpl.rc_context(rc={'font.size': 12}):
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

    return fig