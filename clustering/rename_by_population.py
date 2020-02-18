# -*- coding: utf-8 -*-
"""
Rename microstate trajectory by population.

BSD 3-Clause License
Copyright (c) 2020, Daniel Nagel
All rights reserved.

"""
import argparse
import numpy as np
import pandas as pd


def get_args():
    """Parse all input."""
    parser = argparse.ArgumentParser(
        description='Calculate secondary structure with DSSP.')
    parser.add_argument(
        '-f', '--file', type=str, required=True,
        help='Path to microstate trajectory file.' +
             'Should be single column ascii file.')

    return parser.parse_args().file


def main():
    """Run main routine."""
    # parse arguments
    file_name = get_args()

    # load file
    data = pd.read_csv(file_name,
                       sep=r'\s+',
                       header=None,
                       comment='#').values.flatten()

    # rename states
    data = rename_by_population(data)

    # save renamed traj
    np.savetxt(r'{}Sorted'.format(file_name), data, fmt='%.0f',
               header='States sorted by population.')


def shift_data(data, val_old, val_new):
    """
    Shift data from old to new values.

    The basic function is based on Ashwini_Chaudhary solution:
    https://stackoverflow.com/a/29408060

    Parameters
    ----------
    data : ndarray, list, list of ndarrays
        1D data or a list of data.

    val_old : ndarray or list
        Values in data which should be replaced. All values needs to be within
        the range of `[data.min(), data.max()]`

    val_new : ndarray or list
        Values which will be used instead of old ones.

    Returns
    -------
    data : ndarray
        Shifted data in same shape as input.

    """
    dtype = np.uint32
    # check data-type
    if not np.issubdtype(dtype, np.integer):
        raise TypeError('An unsigned integer type is needed.')

    # offset data and val_old to allow negative values
    offset = np.min(data)

    # convert to np.array
    val_old = (np.asarray(val_old) - offset).astype(dtype)
    val_new = (np.asarray(val_new) - offset).astype(dtype)

    # convert data and shift
    data = (data - offset).astype(dtype)

    # shift data
    conv = np.arange(data.max() + 1, dtype=dtype)
    conv[val_old] = val_new
    data_shifted = conv[data]

    # shift data back
    data_shifted = data_shifted.astype(np.int32) + offset
    return data_shifted


def rename_by_population(traj):
    r"""
    Rename states sorted by their population starting from 1.

    Parameters
    ----------
    traj : ndarray, list of ndarrays
        State trajectory or list of state trajectories.

    Returns
    -------
    traj : ndarray
        Renamed data.

    """
    # get unique states with population
    states, pop = np.unique(traj, return_counts=True)

    # get decreasing order
    idx_sort = np.argsort(pop)[::-1]

    # rename states
    traj_renamed = shift_data(traj,
                              val_old=states[idx_sort],
                              val_new=np.arange(len(states)) + 1)
    return traj_renamed


if __name__ == '__main__':
    main()
