from typing import List


import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matchms import Spectrum
from matchms.plotting import plot_spectrum


from silico_ms.algorithm import FragmentNetowrk



def spectrum_plot(
    spectrum: Spectrum,
    title: str = "",
    grid: bool = False,
    ax: Axes = None
) -> Axes:
    """
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
    """
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


def fragment_network_plot(
    fragment_network: FragmentNetowrk,
    structure_annotated_name: str,
    ax: Axes
) -> Axes:
    """
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