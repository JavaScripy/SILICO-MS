"""Visualization utilities for mass spectrometry data and fragment networks.

This module provides plotting functions to visualize:
1. MS/MS (tandem mass spectrometry) spectra with proper formatting
2. Empty spectrum placeholders when no MS2 data is available
3. Lipid fragment networks with node/edge highlighting for structural annotations

All functions return matplotlib Axes objects for integration into subplot layouts.

Functions:
    spectrum_plot: Plots a matchms Spectrum object with custom title.
    none_spectrum_plot: Creates an empty plot indicating no MS/MS spectrum available.
    fragment_network_plot: Draws a NetworkX fragment network with highlighted structural matches.
"""

import networkx as nx
from matplotlib.axes import Axes
from matchms import Spectrum
from matchms.plotting import plot_spectrum


def spectrum_plot(
    spectrum: Spectrum,
    title: str = "",
    grid: bool = False,
    ax: Axes = None
) -> Axes:
    """Plots a formatted MS/MS spectrum using matchms plotting backend.

    Args:
        spectrum: Input matchms Spectrum object containing MS2 data.
        title: Custom title displayed above the spectrum plot.
        grid: If True, shows grid lines on the plot.
        ax: Matplotlib Axes object to draw the plot on. Creates new Axes if None.

    Returns:
        Axes: Matplotlib Axes with the rendered spectrum plot.
    """
    ax = plot_spectrum(
        spectrum, 
        grid=grid,
        ax=ax
    )
    ax.set_title(title)

    return ax


def none_spectrum_plot(
    title: str = "",
    ax: Axes = None
) -> Axes:
    """Creates a placeholder plot indicating no MS/MS spectrum exists.

    Displays centered text "No MS/MS Spectrum" with labeled axes for consistency.

    Args:
        title: Title displayed above the placeholder plot.
        ax: Matplotlib Axes object to draw on. Creates new Axes if None.

    Returns:
        Axes: Matplotlib Axes with the empty spectrum placeholder.
    """
    ax.set_title(title)
    ax.set_xlabel("m/z")
    ax.set_ylabel("intensity")
    ax.plot([], [])
    ax.figure.text(
        x=0.5,
        y=0.5,
        s="No MS/MS Spectrum",
        ha="center",
        va="center",
        fontsize=14,
        color="black"
    )
    
    return ax

# TODO
'''
def fragment_network_plot(
    fragment_network: nx.Graph,
    structure_annotated_name: str,
    ax: Axes
) -> Axes:
    """Draws a fragment network graph with nodes and edges colored by feature type.

    Highlights nodes that are candidate features (in-degree = 0) and edges
    corresponding to a specific lipid structure annotation.

    Args:
        fragment_network: NetworkX Graph object representing fragment relationships.
        structure_annotated_name: Target lipid structure name to highlight edges.
        ax: Matplotlib Axes object for rendering the network graph.

    Returns:
        Axes: Matplotlib Axes with the rendered fragment network plot.
    """
    labels = [
        f"F{node.get_feature_id()}"
        for node in fragment_network.nodes()
    ]
    candidate_feature_node = [
        node 
        for node in fragment_network.nodes()
        if fragment_network.in_degree(node) == 0
    ]
    highlight_edges = [
        (u, v)
        for u, v, attrs in fragment_network.edges(data=True)
        if (attrs.get("fragment_attr")) and 
        (attrs.get("fragment_attr").structure_annotated_name 
            == structure_annotated_name)
    ]

    node_color = [
        "#C0A060" if node in candidate_feature_node else    
        "#507050"
        for node in fragment_network.nodes()
    ]
    edge_color = [
        "#620404" if (u, v) in highlight_edges else
        "#333333" 
        for u, v in fragment_network.edges()
    ]
    pos = nx.spring_layout(fragment_network)

    nx.draw_networkx(
        G=fragment_network,
        pos=pos,
        node_color=node_color,
        edge_color=edge_color,
        labels=labels,
        ax=ax
    )

    return ax
'''
