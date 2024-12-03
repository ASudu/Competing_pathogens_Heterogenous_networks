import networkx as nx
import random
import math
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import os

def network_SISD(G, N, beta_1, beta_2, beta_3, gamma, delta, mu, S0, I0, max_time):
    g=G.copy()
    pop_size = S0 + I0
    
    i_nodes = list(range(N))
    for i in range(N):
        i_node = random.choice(i_nodes)
        i_nodes.remove(i_node)

        if (i < math.floor(S0/pop_size * N)):
            g.nodes[i_node]['state'] = 'S'
        else:
            g.nodes[i_node]['state'] = 'I1'

    S = len([x for x,y in g.nodes(data=True) if y['state']=='S'])
    I1 = len([x for x,y in g.nodes(data=True) if y['state']=='I1'])
    I2 = len([x for x,y in g.nodes(data=True) if y['state']=='I2'])
    I3 = len([x for x,y in g.nodes(data=True) if y['state']=='I3'])
    D = len([x for x,y in g.nodes(data=True) if y['state']=='D'])
    degrees = [d for n, d in G.degree()]
    first_moment = sum(degrees) / len(degrees)
    second_moment = sum(d**2 for d in degrees) / len(degrees)
    percolation_threshold = first_moment / (second_moment - first_moment)


    states = [(0, S/N, I1/N, I2/N, I3/N, D/N, percolation_threshold)]

    strain_2_time = None
    strain_3_time = None
    S_at_2, I1_at_2, = None, None
    S_at_3, I1_at_3, I2_at_3 = None, None, None

    for t in range(max_time):
        nodes_to_remove = []
        for node in list(g.nodes()):
            if node not in g:
                continue

            all_neighbors = list(g.neighbors(node))
            neighbors = [n for n in all_neighbors if g.nodes[n]['state'] != 'D']
            if neighbors:
                contact = random.choice(neighbors)
            else:
                contact = None
            if g.nodes[node]['state'] == 'S':
                if contact is not None:
                    if (g.nodes[contact]['state'] == 'I1'):
                        if random.random() < beta_1:
                            g.nodes[node]['state'] = 'I1'
                            S -= 1
                            I1 += 1
                    elif (g.nodes[contact]['state'] == 'I2'):
                        if random.random() < beta_2:
                            g.nodes[node]['state'] = 'I2'
                            S -= 1
                            I2 += 1
                    elif (g.nodes[contact]['state'] == 'I3'):
                        if random.random() < beta_3:
                            g.nodes[node]['state'] = 'I3'
                            S -= 1
                            I3 += 1
            elif g.nodes[node]['state'] == 'I1':
                if contact is not None:
                    if (g.nodes[contact]['state'] == 'I2'):
                        if random.random() < beta_2:
                            g.nodes[node]['state'] = 'I2'
                            I1 -= 1
                            I2 += 1
                    elif (g.nodes[contact]['state'] == 'I3'):
                        if random.random() < beta_3:
                            g.nodes[node]['state'] = 'I3'
                            I1 -= 1
                            I3 += 1
                    elif random.random() < gamma:
                        g.nodes[node]['state'] = 'S'
                        I1 -= 1
                        S += 1
                    elif random.random() < delta:
                        nodes_to_remove.append(node)
                        I1 -= 1
                        D += 1
                    elif random.random() < mu:
                        g.nodes[node]['state'] = 'I2'
                        I1 -= 1
                        I2 += 1
                else:
                    if random.random() < gamma:
                        g.nodes[node]['state'] = 'S'
                        I1 -= 1
                        S += 1
                    elif random.random() < delta:
                        nodes_to_remove.append(node)
                        I1 -= 1
                        D += 1
                    elif random.random() < mu:
                        g.nodes[node]['state'] = 'I2'
                        I1 -= 1
                        I2 += 1
            elif g.nodes[node]['state'] == 'I2':
                if contact is not None:
                    if (g.nodes[contact]['state'] == 'I1'):
                        if random.random() < beta_1:
                            g.nodes[node]['state'] = 'I1'
                            I2 -= 1
                            I1 += 1
                    elif (g.nodes[contact]['state'] == 'I3'):
                        if random.random() < beta_3:
                            g.nodes[node]['state'] = 'I3'
                            I2 -= 1
                            I3 += 1
                    elif random.random() < gamma:
                        g.nodes[node]['state'] = 'S'
                        I2 -= 1
                        S += 1
                    elif random.random() < delta:
                        nodes_to_remove.append(node)
                        I2 -= 1
                        D += 1
                    elif random.random() < mu:
                        g.nodes[node]['state'] = 'I3'
                        I2 -= 1
                        I3 += 1
                else:
                    if random.random() < gamma:
                        g.nodes[node]['state'] = 'S'
                        I2 -= 1
                        S += 1
                    elif random.random() < delta:
                        nodes_to_remove.append(node)
                        I2 -= 1
                        D += 1
                    elif random.random() < mu:
                        g.nodes[node]['state'] = 'I3'
                        I2 -= 1
                        I3 += 1
            elif g.nodes[node]['state'] == 'I3':
                if contact is not None:
                    if (g.nodes[contact]['state'] == 'I1'):
                        if random.random() < beta_1:
                            g.nodes[node]['state'] = 'I1'
                            I3 -= 1
                            I1 += 1
                    elif (g.nodes[contact]['state'] == 'I2'):
                        if random.random() < beta_2:
                            g.nodes[node]['state'] = 'I2'
                            I3 -= 1
                            I2 += 1
                    elif random.random() < gamma:
                        g.nodes[node]['state'] = 'S'
                        I3 -= 1
                        S += 1
                    elif random.random() < delta:
                        nodes_to_remove.append(node)
                        I3 -= 1
                        D += 1
                else:
                    if random.random() < gamma:
                        g.nodes[node]['state'] = 'S'
                        I3 -= 1
                        S += 1
                    elif random.random() < delta:
                        nodes_to_remove.append(node)
                        I3 -= 1
                        D += 1

        # Remove nodes marked as dead
        g.remove_nodes_from(nodes_to_remove)
        
        # Calculate percolation threshold
        if len(g) > 0:
            degrees = [d for n, d in g.degree()]
            first_moment = sum(degrees) / len(degrees)
            second_moment = sum(d**2 for d in degrees) / len(degrees)
            percolation_threshold = first_moment / (second_moment - first_moment)
        else:
            largest_cc = []
            percolation_threshold = 0
        
        states.append((t+1, S/N, I1/N, I2/N, I3/N, D/N, percolation_threshold))
    
        if strain_2_time is None and I2 > 0:
            strain_2_time = t
            S_at_2 = S
            I1_at_2 = I1
        if strain_3_time is None and I3 > 0:
            strain_3_time = t
            S_at_3 = S
            I1_at_3 = I1
            I2_at_3 = I2

    # Calculate R0 values
    R0_2 = (mu * I1_at_2 + beta_2 * (S_at_2 + I1_at_2) / N) / ((beta_1 * I1_at_2 / N) + delta + gamma + mu) if I1_at_2 is not None else 0
    R0_3 = (mu * I2_at_3 + beta_3 * (S_at_3 + I2_at_3 + I1_at_3) / N) / ((beta_2 * I2_at_3 / N) + (beta_1 * I1_at_3 / N) + delta + gamma) if I2_at_3 is not None else 0

    return states, R0_2, R0_3, strain_2_time, strain_3_time

networks = ["hetero_cc_0.05","hetero_cc_0.5","hetero_cc_0.34",
            "hetero_deg_var_4.46","hetero_deg_var_8.76","hetero_deg_var_30.35",
            "homog_regular"]

for network in networks:
    path = "./"
    header = os.path.join(path, f"{network}.txt")
    G = nx.read_edgelist(header)

    mapping = {node: i for i, node in enumerate(G.nodes())}
    G = nx.relabel_nodes(G, mapping)

    n_nodes = G.number_of_nodes()

    beta_1 = 0.115
    beta_2 = beta_1 * 9/8
    beta_3 = beta_2 * 10/9
    gamma = 0.05
    delta = 0.0005
    mu = 0.0005
    I0 = 0.05
    S0 = 1 - I0
    max_time = 1000
    runs = 1

    cumulative_S = np.zeros(max_time + 1)
    cumulative_I1 = np.zeros(max_time + 1)
    cumulative_I2 = np.zeros(max_time + 1)
    cumulative_I3 = np.zeros(max_time + 1)
    cumulative_D = np.zeros(max_time + 1)
    cumulative_perc_threshold = np.zeros(max_time + 1)
    cumulative_R0_2 = np.zeros(max_time + 1)
    cumulative_R0_3 = np.zeros(max_time + 1)

    fixed_2 = 0
    fixed_3 = 0

    for i in tqdm(range(runs)):
        results, R0_2, R0_3, strain_2_time, strain_3_time = network_SISD(G, n_nodes, beta_1, beta_2, beta_3, gamma, delta, mu, S0, I0, max_time)
        t, S, I1, I2, I3, D, percolation_threshold = zip(*results)

        cumulative_S = np.add(cumulative_S, S)
        cumulative_I1 = np.add(cumulative_I1, I1)
        cumulative_I2 = np.add(cumulative_I2, I2)
        cumulative_I3 = np.add(cumulative_I3, I3)
        cumulative_D = np.add(cumulative_D, D)
        cumulative_perc_threshold = np.add(cumulative_perc_threshold, percolation_threshold)
        cumulative_R0_2 = np.add(cumulative_R0_2, R0_2)
        cumulative_R0_3 = np.add(cumulative_R0_3, R0_3)

        if I1[max_time] == 0:
            if I2[max_time] == 0:
                fixed_3 += 1
            else:
                fixed_2 += 1

    average_S = cumulative_S / runs
    average_I1 = cumulative_I1 / runs
    average_I2 = cumulative_I2 / runs
    average_I3 = cumulative_I3 / runs
    average_D = cumulative_D / runs
    average_perc_threshold = cumulative_perc_threshold / runs
    average_R0_2 = cumulative_R0_2 / runs
    average_R0_3 = cumulative_R0_3 / runs

    fixed_2 = fixed_2 / runs
    fixed_3 = fixed_3 / runs

    degrees = [d for n, d in G.degree()]
    first_moment = sum(degrees) / len(degrees)
    second_moment = sum(d**2 for d in degrees) / len(degrees)
    
    # R0_1 = beta_1 / (delta + gamma + mu)
    # per_1 = beta_1 / R0_1
    # per_2 = beta_2 / R0_2
    # per_3 = beta_3 / R0_3

    # if ("homog" in network):
    #     lambda_c = 1 / first_moment
    # else:
    #     lambda_c = first_moment / (second_moment - first_moment)

    # r0_eff_1 = R0_1 / lambda_c
    # r0_eff_2 = R0_2 / lambda_c
    # r0_eff_3 = R0_3 / lambda_c

    # per_eff_1 = per_1 * lambda_c
    # per_eff_2 = per_2 * lambda_c
    # per_eff_3 = per_3 * lambda_c

    # epidemic_size = 1 - average_S[max_time]

    # peak_infection = max(average_I1 + average_I2 + average_I3)
    # peak_time = np.argmax(average_I1 + average_I2 + average_I3)

    metrics = f"Betas: {beta_1}, {beta_2}, {beta_3}\n\
    # \n\
    # R0 effective 1: {r0_eff_1}\n\
    # R0 effective 2: {r0_eff_2}\n\
    # R0 effective 3: {r0_eff_3}\n\
    # \n\
    # Epidemic size: {epidemic_size}\n\
    # \n\
    # Peak infection: {peak_infection}\n\
    # Peak time: {peak_time}\n\
    # \n\
    # Percolation threshold 1: {per_eff_1}\n\
    # Percolation threshold 2: {per_eff_2}\n\
    # Percolation threshold 3: {per_eff_3}"

    plt.figure(figsize=(5,4))
    plt.plot(t, average_perc_threshold, label="Percolation Threshold")
    plt.axhline(y=average_perc_threshold[0], color='r', linestyle='--', label='Initial Percolation Threshold')
    if strain_2_time is not None:
        plt.axvline(x=strain_2_time, color='b', linestyle=':', label='Strain 2 Arrival Time')
    if strain_3_time is not None:
        plt.axvline(x=strain_3_time, color='g', linestyle=':', label='Strain 3 Arrival Time')
    plt.xlabel("Time")
    plt.ylabel("Percolation Threshold")
    plt.legend()
    plt.title(f"Network: {network}")
    plt.tight_layout()
    
    output_dir = "results_final"
    output_path = os.path.join(output_dir, network)

    with open(f"{output_path}.txt", "w") as f:
        f.write(metrics)

    plt.savefig(f"{output_path}.png")
