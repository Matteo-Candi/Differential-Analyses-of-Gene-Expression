import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx
from scipy import sparse
import numpy as np
import random

def plot_network(name, hub_color, node_color):
    # Load the CSV file into a pandas DataFrame
    df = pd.read_csv(f'networks/{name}_adj.csv')
    hub_names = []
    with open(f'networks/{name}_names.txt', 'r') as file:
        for line in file:
            hub_name = line.strip()  
            hub_names.append(hub_name)

    random.seed(23)

    # Extract node names from the first row
    nodes = list(df.columns)

    # Create a graph
    G = nx.Graph()

    # Add nodes to the graph with color and name attributes
    for node in nodes:
        color = hub_color if node in hub_names else node_color
        node_size = 150 if node in hub_names else 10  # Larger size for hubs
        G.add_node(node, color=color, name=node, size=node_size)


    overlapping_gene = 'ENSG00000159714.12'
    G.nodes[overlapping_gene]['color'] = 'green'

    if name == 'hub_diff':
        overlapping_genes = ["ENSG00000043143.20", "ENSG00000077782.21", "ENSG00000100311.17", "ENSG00000105556.11", "ENSG00000114554.11", "ENSG00000162972.10", "ENSG00000186205.13", "ENSG00000189241.8"]
        for gene in overlapping_genes:
            G.nodes[gene]['color'] = '#F26419'


    for col, row in df.iterrows():
        for index, value in row.items():
            if value:  
                G.add_edge(nodes[col], index)

    pos = nx.spring_layout(G)  

    radius = .3
    theta = 2 * 3.14159 / len(hub_names)
    h = .1
    for i, hub in enumerate(hub_names):
        x = radius * np.cos(i * theta) + random.uniform(-h, h)
        y = radius * np.sin(i * theta) + random.uniform(-h, h)
        pos[hub] = (x, y)

    other_nodes = [node for node in nodes if node not in hub_names]
    theta = 2 * 3.14159 / len(other_nodes)
    radius = .7
    h = .05
    for i, node in enumerate(other_nodes):
        x = radius * np.cos(i * theta) + random.uniform(-h, h)
        y = radius * np.sin(i * theta) + random.uniform(-h, h)
        pos[node] = (x, y)

    # Draw the graph
    plt.figure(figsize=(5, 4))
    node_colors = [G.nodes[node]['color'] for node in nodes]
    node_sizes = [G.nodes[node]['size'] for node in nodes]

    nx.draw(G, pos, with_labels=False,  node_color=node_colors, node_size=node_sizes, edge_color='whitesmoke', linewidths=.01, arrows=False)



    hub_labels = {overlapping_gene: G.nodes[overlapping_gene]['name'][9:]}
    if name == 'hub_diff':
        hub_labels = {node: G.nodes[node]['name'][9:] for node in overlapping_genes + [overlapping_gene]}

    hub_pos = {node: pos[node] for node in hub_labels.keys()}
    nx.draw_networkx_labels(G, hub_pos, labels=hub_labels, font_color='black', font_size=7)

    plt.title("Network Graph with Hubs Highlighted")
    plt.show()




plot_network("hub_n", hub_color='#FDCA40', node_color='#2EBFA5')
plot_network("hub_c", hub_color='#FDCA40', node_color='#F26419')
plot_network("hub_diff", hub_color='dimgray', node_color='#FDCA40')
