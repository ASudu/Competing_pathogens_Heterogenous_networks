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
    str_edges = [f"{e[0]} {e[1]}\n" for e in list(G.edges)]

    save_path = "../data/networks/" + save_file

    with open(save_path, "w") as f:
        f.writelines(str_edges)

# Havel-Hakimi algorithm
def havel_hakimi(deg_seq: list) -> nx.Graph:
    """Generates a graph from a given degree sequence using the Havel-Hakimi algorithm.

    Args:
        deg_seq (list): Degree sequence of the graph

    Returns:
        nx.Graph: Graph generated from the degree sequence
    """
    G = nx.Graph()
    G.add_nodes_from(range(len(deg_seq)))
    deg_seq = sorted(deg_seq, reverse=True)
    deg_seq = [list(x) for x in zip(range(len(deg_seq)), deg_seq)]
    edges = []

    while len(deg_seq) > 0:
        # Pick the vertex with lowest degree
        v = deg_seq.pop()
        # If the degree is greater than the number of remaining vertices, return None
        if v[1] > len(deg_seq):
            return None
        # Connect the vertex to the v vertices having the highest remaining degree
        for i in range(v[1]):
            deg_seq[i][1] -= 1
            edges.append((v[0], deg_seq[i][0]))
        # Remove the vertex from the list if its degree is 0
        deg_seq = [d for d in deg_seq if d[1] != 0]
        deg_seq = sorted(deg_seq, key=lambda x: x[1], reverse=True)
    
    G.add_edges_from(edges)
    
    return G

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
        n (int): Number of nodes in the network
        mean (float): Average degree of the network
        var (float): Variance of the degree distribution

    Returns:
        nx.Graph: Final network with the required degree distribution
    """
    s = []
    theta = var[0] / (mean - 1)
    k = (mean - 1) / theta
    samples = np.random.gamma(k - 1, theta, len(var) * n)
    
    for i in range(len(var)):
        s.append([int(x) + 1 for x in samples[i*n:(i+1)*n]])
    
    act_means = [np.mean(sample) for sample in s]
    
    for i in range(1, len(var)):
        for j in range(len(s[i]) - 1):
            p = np.random.random() * var[i]
            s[i][j] += (act_means[0] - act_means[i]) - p
            s[i][j + 1] += (act_means[0] - act_means[i]) + p
        
        s[i] = [int(x)+1 for x in s[i]]
    
    print(f"Mean={[np.mean([int(x) + 1 for x in sample]) for sample in s]} Variance={[np.var([int(x) + 1 for x in sample]) for sample in s]}")

    # Generate graph from degree sequence
    G = []
    for i in range(len(var)):
        g = havel_hakimi(s[i])
        G.append(g)
        deg1 = [g.degree(n) for n in g.nodes]
        print(f"Mean: {np.mean(deg1)}, Variance: {np.var(deg1)}")
        print(f"Number of connected components: {nx.number_connected_components(g)}")

    return G

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
def gen_hetero_cc(n:int, k: float, cc:float, p:float=0.1, tol:float=1e-3) -> nx.Graph:
    """Generates a network with target clustering coeffecient.

    Args:
        n (int): Number of nodes in the network
        k (float): average degree of the network
        cc (float): Target clustering coeffecient

    Returns:
        nx.Graph: Required network on n nodes with average degree k and clustering coeffecient cc
    """

    # Take a small world network of same size and average degree
    G = gen_small_world(n,math.floor(k),p)

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

# Generate hetero graph - degree distribution
C, D, E = gen_hetero_deg(n,mean=4,var=[1,3,6])
save_edgelist(C,"hetero_deg_var_1.txt")
save_edgelist(D,"hetero_deg_var_3.txt")
save_edgelist(E,"hetero_deg_var_6.txt")

# # Generate hetero graph - clustering coeffecient
# F = gen_hetero_cc(n,k=4,cc=0.24,p=0.1)
# save_edgelist(F,"hetero_cc_0.24.txt")
# Generate 3 samples from gamma distribution each of size n and each with same mean equal to 4 but different variances 1, 3, and 6 respectively
