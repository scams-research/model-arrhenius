# Copyright (c) Andrew R. McCluskey
# Distributed under the terms of the MIT License
# author: Andrew R. McCluskey (arm61)

from typing import List
import numpy as np
import MDAnalysis as mda
from tqdm import tqdm


def find_step_skip(first_file: str) -> int:
    """
    Determine the step skip from the file.

    :param first_file: Path to first data file.
    :return: step skip value.
    """
    f_open = open(first_file, 'r')
    i = 0
    for line in f_open.readlines():
        if i == 2:
            f_open.close()
            return int(line)
        elif 'TIMESTEP' in line:
            i += 1


def flatten_universe_list(u_list: List[mda.Universe], data_file: str) -> mda.Universe:
    """
    Concatenate a list of Universes.

    :param u_list: List of Universe objects.
    :param data_file: File path for data.
    :return: A single flattened Universe.
    """
    j = 0
    trajectory = np.zeros((sum([len(i.trajectory) for i in u_list]), len(u_list[0].atoms), 3))
    for i in u_list:
        for t in tqdm(i.trajectory):
            trajectory[j] = t.positions
            j += 1
    flat_universe = mda.Universe(data_file, trajectory)
    for t in flat_universe.trajectory:
        t.dimensions = u_list[0].trajectory[0].dimensions + np.array([1, 1, 1, 0, 0, 0])
    return flat_universe


def universe_slice(u: mda.Universe, data_file: str, slice_end: int) -> mda.Universe:
    """
    Concatenate a list of Universes.

    :param u_list: List of Universe objects.
    :param data_file: File path for data.
    :return: A single flattened Universe.
    """
    j = 0
    trajectory = np.zeros((len(u.trajectory[:slice_end]), len(u.atoms), 3))
    for t in tqdm(u.trajectory[:slice_end]):
        trajectory[j] = t.positions
        j += 1
    shorter_universe = mda.Universe(data_file, trajectory)
    for t in shorter_universe.trajectory:
        t.dimensions = u.trajectory[0].dimensions + np.array([1, 1, 1, 0, 0, 0])
    return shorter_universe