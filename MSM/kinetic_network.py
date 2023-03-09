#!/usr/bin/env python3
import msmhelper as mh
import networkx as nx
import numpy as np
import math
import prettypyplot as pplt
import click
from curved_edges import curved_edges
from fa2 import ForceAtlas2
from scipy.spatial import distance_matrix
from matplotlib import pyplot as plt
from matplotlib.colors import to_hex
from matplotlib.collections import LineCollection
from pathlib import Path

pplt.use_style(figsize=1.8, figratio=1)

USE_FA2 = True
DRAW_FLUX = True

def calc_dist(coords2D):
    nodes_dists = distance_matrix(coords2D, coords2D)
    return(nodes_dists)

def check_superposition(node_i, node_j, nodes_dists, node_size):
    if nodes_dists[node_i, node_j] <= node_size[node_i]+node_size[node_j]:
        return(True)

def check_and_shift(old_coords2D, node_i, node_j, nodes_dists, node_size):
    while check_superposition(node_i, node_j, nodes_dists, node_size):
        old_coords2D[node_i, :] += np.sign(old_coords2D[node_i, :]-old_coords2D[node_j, :])*10e-03/nodes_dists[node_i, node_j]
        old_coords2D[node_j, :] -= np.sign(old_coords2D[node_i, :]-old_coords2D[node_j, :])*10e-03/nodes_dists[node_i, node_j]
        nodes_dists = distance_matrix(old_coords2D, old_coords2D)
    new_coords2D = old_coords2D
    return(new_coords2D)

def assign_color(qoft_file, states, traj, levels):
    qoft = mh.opentxt(qoft_file)
    states_qoft = np.array([
        1 - np.mean(qoft[traj == state])
        for state in states
    ])
    states_bin = np.array([_bin(q, levels)
        for q in states_qoft
    ])
    colors_list = [_color(q_bin, levels) for q_bin in states_bin]
    return(colors_list)

def _color(val, levels):
    cmap = plt.get_cmap('plasma', levels)
    return to_hex(
        cmap(val),
    )

def _bin(val, levels):
    # get bin
    bins = np.linspace(0, 1, levels + 1)

    for rlower, rhigher in zip(bins[:-1], bins[1:]):
        if rlower <= val <= rhigher:
            return rlower

    return bins[-1]

def get_luminance(hex_color):
    color = hex_color[1:]
    hex_red = int(color[0:2], base=16)
    hex_green = int(color[2:4], base=16)
    hex_blue = int(color[4:6], base=16)
    return(hex_red*0.2126 + hex_green*0.7152 + hex_blue*0.0722)


@click.command()
@click.option(
    '--states_traj',
    type=str,
    help='path to macrostates trajectory',
)
@click.option(
    '-u',
    multiple=True,
    type=int,
    help='list of states in unfolded basin',
)
@click.option(
    '-f',
    multiple=True,
    type=int,
    help='list of states in folded basin',
)
@click.option(
    '--qoft',
    required=True,
    type=click.Path(exists=True),
    help='Path to fraction of native contacts time evolution',
)
@click.option(
    '--set_min_node_size',
    type=bool,
    default=True,
    required=False,
    help='whether to set a lower threshold for the nodes size',
)
@click.option(
    '--tlag',
    type=int,
    help='lag time in frames for respective state traj',
)
def draw_knetwork(states_traj, tlag, qoft, set_min_node_size, u, f):
    
    path = Path(f'{states_traj}')
    _, ax = plt.subplots()
    traj = mh.opentxt(states_traj).astype(int)
    tmat, states = mh.estimate_markov_model(traj, tlag)
    n_nodes = len(np.unique(states))
    color_list = assign_color(qoft, states, traj, levels=10)

    # get detailed balance
    pop_eq = mh.equilibrium_population(tmat)
    mat = tmat * pop_eq[:, np.newaxis]
    mat = 0.5 * (mat + mat.T)

    # prepare mats for networkx
    mat[np.diag_indices_from(mat)] = 0
    mat[mat < 2e-5] = 0

    # node size
    node_size = 1000 * np.log(pop_eq + 1)

    # set minimum node size
    if set_min_node_size:
        node_size = np.where(node_size<(np.min(node_size)+np.max(node_size))/2,
                            .7*(np.min(node_size)+np.max(node_size))/2, node_size)

    graph = nx.from_numpy_array(mat, create_using=nx.Graph)

    # get position
    # initial guess of simple spring model
    pos = nx.spring_layout(
        graph, fixed=None, iterations=1000, threshold=1e-4, scale=0.1, weight='weight',
    )
    if USE_FA2:
        # improve pos by forceatlas2
        forceatlas2 = ForceAtlas2(
            adjustSizes=False, verbose=False, strongGravityMode=True,
        )

        pos = forceatlas2.forceatlas2_networkx_layout(
            graph, pos=pos, iterations=1000,
        )
        coords2D = np.asarray(list(pos.values()))
        nodes_dists = calc_dist(np.asarray(coords2D))

        for i in range(n_nodes):
            for j in range(i+1, n_nodes):
                new_coords2D = check_and_shift(coords2D, i, j, nodes_dists, node_size)
        # rotate network so that the folded basin - native basin axis is parallel to the x axis
        coords_u = new_coords2D[u, :]
        coords_f = new_coords2D[f, :]
        a = np.mean(coords_u[:, 0])-np.mean(coords_f[:, 0])
        b = np.mean(coords_u[:, 1])-np.mean(coords_f[:, 1])
        theta = math.atan2(b, a)
        rotated_coords = []
        for i in range(n_nodes):     
            x = new_coords2D[i, 0] * math.cos(-theta) - new_coords2D[i, 1] * math.sin(-theta)
            y = new_coords2D[i, 0] * math.sin(-theta) + new_coords2D[i, 1] * math.cos(-theta)
            rotated_coords.append((x, y))
        new_coords2D = rotated_coords
        keys = list(pos.keys())
        pos = dict(zip(keys, new_coords2D))

    if DRAW_FLUX:
        edge_width = 0.1 + 300 * np.array(
            [graph[i][j]['weight'] for i, j in graph.edges],
        )
        curves = curved_edges(graph, pos)
        lc = LineCollection(
            curves, color='black', linewidth=edge_width, alpha=1,
        )
        ax.add_collection(lc)

    if not DRAW_FLUX:
        # create directed graph to draw edges
        digraph = nx.from_numpy_array(tmat, create_using=nx.DiGraph)
        edge_width = 0.2 + 5 * np.array(
            [digraph[i][j]['weight'] for i, j in digraph.edges],
        )
        nx.draw_networkx_edges(
            digraph,
            arrowstyle="-",
            pos=pos,
            connectionstyle='arc3,rad=0.4',
            width=edge_width,
            edge_color='black',
            node_size=node_size,
            arrowsize=3,
        )

    nx.draw_networkx_nodes(
        graph,
        pos=pos,
        node_color=color_list,
        node_size=node_size,
        linewidths=0.55,
        edgecolors='black'
    )
    # write node labels
    for node_idx, (x, y) in pos.items():
        luminance = get_luminance(color_list[node_idx])
        if luminance < 140 and set_min_node_size:
            c_text = 'white'
        else:
            c_text = 'black'
        pplt.text(x, y, states[node_idx], contour=False, fontsize='medium', color=c_text)
    # calc limits
    lims = np.array([
        (
            x - max(node_size),
            x + max(node_size),
            y - max(node_size),
            y + max(node_size),
        )
        for n, (x, y) in pos.items()
    ])
    ax.set_xlim(lims[:, 0].min(), lims[:, 1].max())
    ax.set_ylim(lims[:, 2].min(), lims[:, 3].max())

    ax.set_axis_off()
    pplt.savefig(f'images/{path.name}.dynamical_network.{tlag}f.pdf')




if __name__=='__main__':
    draw_knetwork()
