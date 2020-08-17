# -*- coding: utf-8 -*-
"""Rename microstate trajectory by population.

BSD 3-Clause License
Copyright (c) 2020, Daniel Nagel
All rights reserved.

"""
import click

import msmhelper as mh


@click.command(no_args_is_help='-h')
@click.option(
    '--filename',
    '-f',
    required=True,
    type=click.Path(exists=True),
    help='Path to microstate trajectory file (single column ascii file).',
)
@click.option(
    '--output',
    '-o',
    type=click.Path(),
    help='Path to output. Default append "Sorted" to input name.',
)
def rename_by_population(filename, output):
    """Rename microstates by decreasing population."""
    # load file
    traj = mh.openmicrostates(filename)[0]

    # rename states
    traj = mh.rename_by_population(traj)

    # output name
    if not output:
        output = r'{0}Sorted'.format(filename)

    # save renamed traj
    mh.savetxt(
        output,
        traj,
        fmt='%.0f',  # noqa: WPS323
        header='States sorted by population.',
    )


if __name__ == '__main__':
    rename_by_population()
