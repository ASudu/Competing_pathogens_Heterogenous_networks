# In this file, we generate the 3 types of networks:
# a) Homogeneous networks
# b) Heterogeneous networks with different degree distribution
# c) Heterogeneous networks with different clustering coeffecient

import networkx as nx
import numpy as np
import math
from typing import Optional, Literal
import os

def save_edgelist(G: nx.Graph, save_file: str)-> None:
    """Saves the edgelist of the graph G.

    Args:
        G (nx.Graph): Graph to save
        save_file (str): File name for the edgelist file is saved.
    """
    str_edges = [f"{e[0]} {e[1]}" for e in list(G.edges)]

    save_path = "../data/networks/" + save_file

    with open(save_path, "w") as f:
        f.writelines(str_edges)

# Homogeneous network
def gen_homog(n: int, k: Optional[int] = None) -> nx.Graph:
    """
    Generates a k-regular graph on n nodes or a uniform random graph on n nodes.

    Args:
        n (int): Number of nodes required in the network
        k (Optional[int], optional): Degree of every node in  the regular graph (needed when kind=0). Defaults to None.

    Returns:
        nx.Graph: _description_
    """

    # If k is None or n*k is odd or k>=n
    if not k or ((n*k)%2) or k>=n:
        return nx.complete_graph(n)
    
    # If k is valid
    else:
        return nx.random_regular_graph(k,n)
 
# Heterogeneous network - degree distribution
def gen_hetero_deg(n:int, mean: float, var: float) -> nx.Graph:
    """Generates a network with skewed degree distribution dictated by an
    average degree and variance. The degrees are sampled from 

    Args:
        n (int): _description_
        mean (float): _description_
        var (float): _description_

    Returns:
        nx.Graph: _description_
    """
    scale_param = var/(mean-1) # k of gamma dist
    shape_param = (mean-1)/scale_param # theta of gamma dist

    samples = np.random.gamma(scale_param - 1,shape_param,n)
    deg = [int(x) + 1 for x in samples]

    return nx.random_degree_sequence_graph(deg)

def gen_small_world(n:int, k:int, p:float) -> nx.Graph:
    """Generates a small world network on n nodes with average degree k and rewiring probability p.

    Args:
        n (int): Number of nodes in the network
        k (int): Average degree of the network
        p (float): Rewiring probability

    Returns:
        nx.Graph: Required small world network
    """
    return nx.watts_strogatz_graph(n,k,p)

# Heterogeneous network - clustering coeffecient
def gen_hetero_cc(n:int, k: float, cc:float, tol:float=1e-3) -> nx.Graph:
    """Generates a network with target clustering coeffecient.

    Args:
        n (int): Number of nodes in the network
        k (float): average degree of the network
        cc (float): Target clustering coeffecient

    Returns:
        nx.Graph: Required network on n nodes with average degree k and clustering coeffecient cc
    """

    # Take a small world network of same size and average degree
    G = gen_small_world(n,math.floor(k),0.2)

    gcc = nx.transitivity(G) # Global clustering coeffecient
    nodes = list(G.nodes)
    seed = np.random.default_rng() # Random number generator
    while abs(gcc - cc) > tol:
        print(f"GCC: {gcc}")
        # rewire edges from each node
        # loop over all nodes in order (label) and neighbors in order (distance)
        # no self loops or multiple edges allowed
        for j in range(1, k // 2 + 1):  # outer loop is neighbors
            targets = nodes[j:] + nodes[0:j]  # first j nodes are now last in list
            # inner loop in node order
            for u, v in zip(nodes, targets):
                if seed.random() < p:
                    w = seed.choice(nodes)
                    # Enforce no self-loops or multiple edges
                    while w == u or G.has_edge(u, w):
                        w = seed.choice(nodes)
                        if G.degree(u) >= n - 1:
                            break  # skip this rewiring
                    else:
                        G.remove_edge(u, v)
                        G.add_edge(u, w)
        return G



# Make data folder
path = "../data"

if not os.path.exists(path):
    os.mkdir(path)

if not os.path.exists(os.path.join(path,"networks")):
    os.mkdir(os.path.join(path,"networks"))

# Generate Homogeneous graph
n = 1000
# A = gen_homog(n)
# save_edgelist(A,"homog_comp.txt")
# print("Generated complete graph")

# B = gen_homog(n,4) # same as mean=4 and var=0
# save_edgelist(B,"homog_regular.txt")
# print("Generated regular graph")

# # Generate hetero graph - degree distribution
# C = gen_hetero_deg(n,4,1)
# save_edgelist(C,"hetero_deg_var_1.txt")
# print("Generated for mu=4 and var=1")

# D = gen_hetero_deg(n,4,3)
# save_edgelist(D,"hetero_deg_var_3.txt")
# print("Generated for mu=4 and var=3")

# E = gen_hetero_deg(n,4,6)
# save_edgelist(E,"hetero_deg_var_6.txt")
# print("Generated for mu=4 and var=6")

# Generate hetero graph - clustering coeffecient
F = gen_hetero_cc(n,4,0.24)
save_edgelist(F,"hetero_cc_0.24.txt")