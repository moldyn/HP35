#!/usr/bin/env python3
"""Generate Dendrogram from MPP and Lump."""
# ~~~ IMPORT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import sys
from functools import lru_cache

import click
import msmhelper as mh
import numpy as np
import prettypyplot as pplt
from matplotlib import pyplot as plt
from matplotlib import ticker
from matplotlib.cm import ScalarMappable
from matplotlib.collections import LineCollection
from matplotlib.colors import to_hex, Normalize, LinearSegmentedColormap
from scipy.cluster.hierarchy import dendrogram
from tqdm import tqdm

sys.setrecursionlimit(20000)
MIN_EDGE_WEIGHT = 6 * 0.025


@click.command(no_args_is_help='-h')
@click.option(
    '--linkage',
    'linkage_file',
    required=True,
    type=click.Path(exists=True),
    help=(
        'Path to linkage matrix created by MPP named '
        '"{basename}_transitions.dat".'
    ),
)
@click.option(
    '--color-threshold',
    'threshold',
    default=1.0,
    type=click.FloatRange(min=0, max=1),
    help='All links above threshold are colored in dark gray.',
)
@click.option(
    '--tlag',
    required=True,
    type=click.IntRange(min=1),
    help='Lagtime in frames',
)
@click.option(
    '--cut-params',
    default=(0.005, 0.2),
    type=click.FloatRange(min=0, max=0.99),
    nargs=2,
    help='Tuple defining the cut parameters, (pop, q) both in [0, 1).',
)
@click.option(
    '--state-traj',
    'state_traj',
    required=True,
    type=click.Path(exists=True),
    help='Used to apply MPP+ automated lumping.',
)
@click.option(
    '--fraction-of-native-contacts',
    'qtraj',
    required=True,
    type=click.Path(exists=True),
    help='File of holding fraction of native contacts Q.',
)
@click.option(
    '--hide-labels',
    is_flag=True,
    help='If set, no xlabels are plotted.',
)
def plot_dendrogram(
    linkage_file,
    threshold,
    tlag,
    cut_params,
    state_traj,
    qtraj,
    hide_labels,
):
    """Plot MPP result."""
    # parse input and create output basename
    pop_thr, qmin_thr = cut_params
    output_file = (
        f'{linkage_file}.renamed_by_q.pop{pop_thr:.3f}_qmin{qmin_thr:.2f}'
    )

    # setup matplotlib
    pplt.use_style(figsize=2.6, figratio='golden', true_black=True)

    # load transitions and sort them
    transitions = np.loadtxt(linkage_file)
    (
        linkage_mat,
        states_idx_to_microstates,
        states_idx_to_rootstates,
        labels,
    ) = _transitions_to_linkage(transitions, qmin=0)

    # get states
    nstates = len(linkage_mat) + 1
    states = np.unique(linkage_mat[:, :2].astype(int))

    # replace state names by their indices
    transitions_idx, states_idx = mh.rename_by_index(
        transitions[:, :2].astype(int),
        return_permutation=True,
    )
    transitions[:, :2] = transitions_idx

    # estimate population of states
    traj = mh.opentxt(state_traj)
    microstates, counts = np.unique(traj, return_counts=True)
    pops = counts / len(traj)
    pops = {
        idx_state: np.sum([
            pops[microstates == state]
            for state in states_idx_to_microstates[idx_state]
        ])
        for idx_state in states
    }
    pops[2 * (nstates - 1)] = 1.0

    # use population as edge widths
    edge_widths = {
        state: 6 * pops[state] for state in range(2 * nstates - 1)
    }

    # find optimal cut
    macrostates, macrostates_assignment = mpp_plus_cut(
        states_idx_to_rootstates=states_idx_to_rootstates,
        states_idx_to_microstates=states_idx_to_microstates,
        linkage_mat=linkage_mat,
        microstates=microstates,
        pops=pops,
        pop_thr=pop_thr,
        qmin_thr=qmin_thr,
    )
    n_macrostates = len(macrostates_assignment)

    # estimate Q(state)
    q_of_t = mh.opentxt(qtraj, dtype=np.float32)
    q_state = {
        idx_state: _fraction_of_native_contacts(
            states_idx_to_microstates[idx_state],
            q_of_t,
            traj,
        )
        for idx_state in tqdm(states)
    }

    # define colors
    colors = {
        idx_state: _color_by_q(q_state[idx_state])
        for idx_state in states
    }
    # add global value
    colors[2 * (nstates - 1)] = _color_by_q(1.0)

    fig, (ax, ax_mat) = plt.subplots(
        2,
        1,
        gridspec_kw={
            'hspace': 0.05 if hide_labels else 0.3,
            'height_ratios': [9, 1],
        },
    )
    # hide spines of lower mat
    for key, spine in ax_mat.spines.items():
        spine.set_visible(False)

    dendrogram_dict = _dendrogram(
        ax=ax,
        linkage_mat=linkage_mat,
        colors=colors,
        threshold=threshold,
        labels=labels,
        qmin=0,
        edge_widths=edge_widths,
    )

    # plot legend
    cmap, bins = _color_by_q(None)
    norm = Normalize(bins[0], bins[-1])
    label = r'$\langle Q \rangle_\text{state} $'

    cmappable = ScalarMappable(norm, cmap)
    plt.sca(ax)
    pplt.colorbar(cmappable, width='5%', label=label, position='top')

    yticks = np.arange(0.5, 1.5 + n_macrostates)
    xticks = 10 * np.arange(0, nstates + 1)
    cmap = LinearSegmentedColormap.from_list(
        'binary', [(0, 0, 0, 0), (0, 0, 0, 1)],
    )

    # permute macrostate assignment and label them
    macrostates_assignment = macrostates_assignment.T[
        dendrogram_dict['leaves']
    ].T
    macrostates = macrostates[dendrogram_dict['leaves']]
    microstates = microstates[dendrogram_dict['leaves']]

    # apply dynamical correction of minor branches
    dyn_corr_macrostates = mpp_plus_dyn_cor(
        macrostates=macrostates,
        microstates=microstates,
        n_macrostates=n_macrostates,
        pops=pops,
        traj=traj,
        tlag=tlag,
    )

    # rename macrostates by fraction of native contacts
    macrotraj = mh.shift_data(traj, microstates, dyn_corr_macrostates)
    macrostates_q = [
        1 - _fraction_of_native_contacts([state], q_of_t, macrotraj)
        for state in np.unique(macrostates)
    ]
    macroperm = np.unique(dyn_corr_macrostates)[np.argsort(macrostates_q)]
    dyn_corr_macrostates = mh.shift_data(
        dyn_corr_macrostates, macroperm, np.unique(dyn_corr_macrostates),
    )

    mh.savetxt(
        f'{output_file}.macrostates',
        np.array([microstates, dyn_corr_macrostates]).T,
        header='microstates macrostates',
        fmt='%.0f',
    )

    macrotraj = mh.shift_data(traj, microstates, dyn_corr_macrostates)

    # print final sorting order
    macrostates_q = [
        _fraction_of_native_contacts([state], q_of_t, macrotraj)
        for state in np.unique(macrostates)
    ]
    mh.savetxt(
        f'{output_file}.macrotraj',
        macrotraj,
        header='macrostates',
        fmt='%.0f',
    )

    # recalculate macrostates_assignment)
    for idx, mstate in enumerate(np.unique(dyn_corr_macrostates)):
        macrostates_assignment[idx] = dyn_corr_macrostates == mstate

    xvals = 0.5 * (xticks[:-1] + xticks[1:])
    for idx, assignment in enumerate(macrostates_assignment):
        xmean = np.median(xvals[assignment == 1])

        pplt.text(
            xmean,
            yticks[idx] - (yticks[1] - yticks[0]),
            f'{idx + 1:.0f}',
            ax=ax_mat,
            va='top',
            contour=True,
            size='small',
        )

    ax_mat.pcolormesh(
        xticks,
        yticks,
        macrostates_assignment,
        snap=True,
        cmap=cmap,
        vmin=0,
        vmax=1,
    )
    # set x-labels
    ax_mat.set_yticks(yticks)
    ax_mat.set_yticklabels([])
    ax_mat.grid(visible=True, axis='y', ls='-', lw=0.5)
    ax_mat.tick_params(axis='y', length=0, width=0)
    ax_mat.set_xlim(ax.get_xlim())
    ax.set_xlabel('')
    ax_mat.set_xlabel('macrostates')
    ax_mat.set_ylabel('')
    fig.align_ylabels([ax, ax_mat])

    ax_mat.set_xticks(np.arange(0.5, 0.5 + len(states)))

    if hide_labels:
        for axes in (ax, ax_mat):  # if statemat_file else [ax]:
            axes.set_xticks([])
            axes.set_xticks([], minor=True)
            axes.set_xticklabels([])
            axes.set_xticklabels([], minor=True)

    pplt.savefig(f'{output_file}.pdf')


def _color_by_q(q, qmax=1, qmin=0, steps=10):
    cmap = plt.get_cmap('plasma_r', steps)
    colors = [cmap(idx) for idx in range(cmap.N)]

    bins = np.linspace(
        qmin, qmax, steps + 1,
    )

    if q is None:
        return cmap, bins

    for color, rlower, rhigher in zip(colors, bins[:-1], bins[1:]):
        if rlower <= q <= rhigher:
            return color
    return 'k'


def _fraction_of_native_contacts(states, q_of_t, traj):
    """Mean fraction of native contacts per state."""
    if len(states):
        mask = np.full(q_of_t.shape[0], False)
        for state in states:
            mask = np.logical_or(
                mask,
                traj == state,
            )
        cs = q_of_t[mask]
    else:
        cs = q_of_t
    return np.mean(cs)


def _transitions_to_linkage(trans, *, qmin=0.0):
    """Convert transition matrix to linkage matrix.

    Parameters
    ----------
    transitions: ndarray of shape (nstates - 1, 3)
        Three column: merged state, remaining state, qmin lebel.

    qmin: float [0, 1]
        Qmin cut-off. Returns only sublinkage-matrix.

    """
    transitions = np.copy(trans)
    states = np.unique(transitions[:, :2].astype(int))

    # sort by merging qmin level
    transitions = transitions[
        np.argsort(transitions[:, 2])
    ]

    # create linkage matrix
    mask_qmin = transitions[:, 2] > qmin
    nstates_qmin = np.count_nonzero(mask_qmin) + 1
    linkage_mat = np.zeros((nstates_qmin - 1, 4))

    # replace state names by their indices
    transitions_idx, states_idx = mh.rename_by_index(
        transitions[:, :2][mask_qmin].astype(int),
        return_permutation=True,
    )
    transitions[:, :2][mask_qmin] = transitions_idx
    linkage_mat[:, :3] = transitions[mask_qmin]

    # holds for each state (index) a list corresponding to the microstates
    # it consist of.
    states_idx_to_microstates = {
        idx: [
            state,
            *transitions[~mask_qmin][:, 0][
                transitions[~mask_qmin][:, 1] == state
            ].astype(int),
        ]
        for idx, state in enumerate(states_idx)
    }
    states_idx_to_rootstates = {
        idx: [idx]
        for idx, _ in enumerate(states_idx)
    }

    for idx, nextstate in enumerate(
        range(nstates_qmin, 2 * nstates_qmin - 1),
    ):
        statefrom, stateto = linkage_mat[idx, :2].astype(int)
        states_idx_to_microstates[nextstate] = [
            *states_idx_to_microstates[stateto],
            *states_idx_to_microstates[statefrom],
        ]
        states_idx_to_rootstates[nextstate] = [
            *states_idx_to_rootstates[stateto],
            *states_idx_to_rootstates[statefrom],
        ]

        states = linkage_mat[idx, :2].astype(int)
        for state in states:
            linkage_mat[idx + 1:, :2][
                linkage_mat[idx + 1:, :2] == state
            ] = nextstate

    labels = [
        states_idx_to_microstates[idx][0]
        for idx in range(nstates_qmin)
    ]

    return (
        linkage_mat,
        states_idx_to_microstates,
        states_idx_to_rootstates,
        labels,
    )


def _dendrogram(
    *, ax, linkage_mat, colors, threshold, labels, qmin, edge_widths,
):
    nstates = len(linkage_mat) + 1
    # convert color dictionary to array
    colors_arr = np.array(
        [
            to_hex(colors[state]) for state in range(2 * nstates - 1)
        ],
        dtype='<U7',
    )

    dendrogram_dict = dendrogram(
        linkage_mat,
        leaf_rotation=90,
        get_leaves=True,
        color_threshold=1,
        link_color_func=lambda state_idx: colors_arr[state_idx],
        no_plot=True,
    )
    _plot_dendrogram(
        icoords=dendrogram_dict['icoord'],
        dcoords=dendrogram_dict['dcoord'],
        ivl=dendrogram_dict['ivl'],
        color_list=dendrogram_dict['color_list'],
        threshold=threshold,
        ax=ax,
        colors=colors_arr,
        labels=labels,
        qmin=qmin,
        edge_widths=edge_widths,
    )

    ax.set_ylabel(r'metastability $Q_\text{min}$')
    ax.set_xlabel('microstates')
    ax.grid(visible=False, axis='x')

    return dendrogram_dict


def _show_xlabels(*, ax, states_perm):
    """Show the xticks together with the corresponding state names."""
    # undo changes of scipy dendrogram
    xticks = ax.get_xticks()
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    for line in ax.get_xticklines():
        line.set_visible(True)

    for is_major, length_scale in ((True, 4), (False, 1)):
        ax.tick_params(
            axis='x',
            length=length_scale * plt.rcParams['xtick.major.size'],
            labelrotation=90,
            pad=2,
            labelsize='xx-small',
            width=plt.rcParams['xtick.major.width'],
            which='major' if is_major else 'minor',
            top=False,
        )
        offset = 0 if is_major else 1
        ax.set_xticks(xticks[offset::2], minor=not is_major)
        ax.set_xticklabels(states_perm[offset::2], minor=not is_major)


def _plot_dendrogram(
    *,
    icoords,
    dcoords,
    ivl,
    color_list,
    threshold,
    ax,
    colors,
    labels,
    qmin,
    edge_widths,
):
    """Plot dendrogram with colors at merging points."""
    threshold_color = to_hex('pplt:grey')
    # Independent variable plot width
    ivw = len(ivl) * 10
    # Dependent variable plot height
    dvw = 1.05

    iv_ticks = np.arange(5, len(ivl) * 10 + 5, 10)

    ax.set_ylim([qmin, dvw])
    ax.set_xlim([-0.005 * ivw, 1.005 * ivw])
    ax.set_xticks(iv_ticks)

    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticklabels(np.asarray(labels)[np.asarray(ivl).astype(int)])

    get_ancestor = _get_ancestor_func(
        icoords, dcoords, ivl,
    )

    # Let's use collections instead. This way there is a separate legend item
    # for each tree grouping, rather than stupidly one for each line segment.
    colors_used = np.unique(colors)
    color_to_lines = {color: [] for color in (*colors_used, threshold_color)}
    width_to_lines = {color: [] for color in (*colors_used, threshold_color)}
    for xline, yline, color in zip(icoords, dcoords, color_list):
        if np.max(yline) <= threshold:
            # split into left and right to color separately
            xline_l = [*xline[:2], np.mean(xline[1:3])]
            xline_r = [np.mean(xline[1:3]), *xline[2:]]

            color_l = _get_ancestor_color(
                icoords,
                dcoords,
                xline[0],
                yline[0],
                color_list,
                ivl,
                colors,
            )
            ancestors_l = get_ancestor(xline[0], yline[0])
            weight_l = np.sum([
                edge_widths[ancestor] for ancestor in ancestors_l
            ])
            color_r = _get_ancestor_color(
                icoords,
                dcoords,
                xline[3],
                yline[3],
                color_list,
                ivl,
                colors,
            )
            ancestors_r = get_ancestor(xline[3], yline[3])
            weight_r = np.sum([
                edge_widths[ancestor] for ancestor in ancestors_r
            ])
            color_to_lines[color_l].append(list(zip(xline_l, yline[:3])))
            width_to_lines[color_l].append(
                max(weight_l, MIN_EDGE_WEIGHT),
            )
            color_to_lines[color_r].append(list(zip(xline_r, yline[1:])))
            width_to_lines[color_r].append(
                max(weight_r, MIN_EDGE_WEIGHT),
            )

        elif np.min(yline) >= threshold:
            color_to_lines[threshold_color].append(list(zip(xline, yline)))
        else:
            yline_bl = [yline[0], np.max([threshold, yline[1]])]
            yline_br = [np.max([threshold, yline[2]]), yline[3]]
            color_to_lines[color].append(list(zip(xline[:2], yline_bl)))
            color_to_lines[color].append(list(zip(xline[2:], yline_br)))

            yline_thr = np.where(np.array(yline) < threshold, threshold, yline)
            color_to_lines[threshold_color].append(list(zip(xline, yline_thr)))

    # Construct the collections.
    colors_to_collections = {
        color: LineCollection(
            color_to_lines[color], colors=(color,),
            linewidths=width_to_lines[color],
        )
        for color in (*colors_used, threshold_color)
    }

    # Add all the groupings below the color threshold.
    for color in colors_used:
        ax.add_collection(colors_to_collections[color])
    # If there's a grouping of links above the color threshold, it goes last.
    ax.add_collection(colors_to_collections[threshold_color])


def _get_ancestor_color(
    xlines, ylines, xval, yval, color_list, ivl, colors,
):
    """Get the color of the ancestors."""
    # if ancestor is root
    if not yval:
        ancestor = int(ivl[int((xval - 5) // 10)])
        return colors[ancestor]

    # find ancestor color
    xy_idx = np.argwhere(
        np.logical_and(
            np.array(ylines)[:, 1] == yval,
            np.array(xlines)[:, 1:3].mean(axis=1) == xval,
        ),
    )[0][0]
    return color_list[xy_idx]


def _get_ancestor_func(
    xlines, ylines, ivl,
):
    """Get the color of the ancestors."""
    @lru_cache(maxsize=1024)
    def _get_ancestor_rec(xval, yval):
        # if ancestor is root
        if not yval:
            ancestor = int(ivl[int((xval - 5) // 10)])
            return (ancestor, )

        # find ancestor color
        xy_idx = np.argwhere(
            np.logical_and(
                np.array(ylines)[:, 1] == yval,
                np.array(xlines)[:, 1:3].mean(axis=1) == xval,
            ),
        )[0][0]
        xleft, yleft = xlines[xy_idx][0], ylines[xy_idx][0]
        xright, yright = xlines[xy_idx][3], ylines[xy_idx][3]

        return (
            *_get_ancestor_rec(xleft, yleft),
            *_get_ancestor_rec(xright, yright),
        )

    return _get_ancestor_rec


def state_sequences(macrostates, state):
    """Get continuous index sequences of macrostate in mstate assignment."""
    state_idx = np.where(macrostates == state)[0]
    idx_jump = state_idx[1:] - state_idx[:-1] != 1
    return np.array_split(
        state_idx,
        np.nonzero(idx_jump)[0] + 1,
    )


def mpp_plus_cut(
    *,
    states_idx_to_rootstates,
    states_idx_to_microstates,
    linkage_mat,
    microstates,
    pops,
    pop_thr,
    qmin_thr,
):
    """Apply MPP+ step1: Identify branches."""
    nstates = len(linkage_mat) + 1

    macrostates_set = [
        set(states_idx_to_rootstates[2 * (nstates - 1)]),
    ]
    macrostates_leaf_set = [
        set(states_idx_to_microstates[2 * (nstates - 1)]),
    ]
    for state_i, state_j, qmin in reversed(linkage_mat[:, :3]):
        if pops[state_i] > pop_thr and qmin > qmin_thr:
            mstate_i = set(states_idx_to_rootstates[state_i])
            macrostates_set = [
                mstate - mstate_i
                for mstate in macrostates_set
            ]
            macrostates_set.append(mstate_i)

            mstate_leaf_i = set(states_idx_to_microstates[state_i])
            macrostates_leaf_set = [
                mstate - mstate_leaf_i
                for mstate in macrostates_leaf_set
            ]
            macrostates_leaf_set.append(mstate_leaf_i)

    n_macrostates = len(macrostates_set)
    macrostates_assignment = np.zeros((n_macrostates, nstates))
    for idx, mstate in enumerate(macrostates_set):
        macrostates_assignment[idx][list(mstate)] = 1

    macrostates = np.empty(len(microstates), dtype=np.int64)
    for idx, microstate in enumerate(microstates):
        for idx_m, macroset in enumerate(macrostates_leaf_set):
            if microstate in macroset:
                macrostates[idx] = idx_m + 1
                break
        else:
            print(f'{microstate} not in macrostate')

    return macrostates, macrostates_assignment


def mpp_plus_dyn_cor(
    *,
    macrostates,
    microstates,
    n_macrostates,
    pops,
    traj,
    tlag,
):
    """Apply MPP+ step2: Dynamically correct minor branches."""
    # fix dynamically missassigned single-state branches
    # identify them
    dyn_corr_macrostates = macrostates[:]
    for mstate in np.unique(macrostates):
        idx_sequences = state_sequences(macrostates, mstate)
        if len(idx_sequences) > 1:
            highest_pop_sequence = np.argmax([
                np.sum([
                    pops[s] for s in microstates[seq]
                ]) for seq in idx_sequences
            ])
            idx_sequences = [
                seq for idx, seq in enumerate(idx_sequences)
                if idx != highest_pop_sequence
            ]
            for seq in idx_sequences:
                largest_state = np.max(dyn_corr_macrostates)
                for newstate, seq_idx in enumerate(
                    seq,
                    largest_state + 1,
                ):
                    dyn_corr_macrostates[seq_idx] = newstate

    # dynamically reassign all new state to previous macrostates
    mstates = np.unique(dyn_corr_macrostates)
    while len(mstates) > n_macrostates:
        tmat, mstates = mh.msm.estimate_markov_model(
            mh.shift_data(traj, microstates, dyn_corr_macrostates),
            lagtime=tlag,
        )

        # sort new states by increasing metastability
        qs = np.diag(tmat)[n_macrostates:]
        idx_sort = np.argsort(qs)
        newstates = mstates[n_macrostates:][idx_sort]

        deletestate = newstates[0]

        # reassign them
        idx = np.where(mstates == deletestate)[0][0]
        idxs_to = np.argsort(tmat[idx])[::-1]
        for idx_to in idxs_to:
            if idx_to == idx:
                continue
            dyn_corr_macrostates[
                dyn_corr_macrostates == deletestate
            ] = mstates[idx_to]
            break

        mstates = np.unique(dyn_corr_macrostates)

    return dyn_corr_macrostates


if __name__ == '__main__':
    plot_dendrogram()
