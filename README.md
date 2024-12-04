# Competing_pathogens_Heterogenous_networks

## Introduction

In this project, we study the following:
- How sequential mutant strains of a virus, each mutnt more transmissible than its parent, compete in population
- How the heterogeneity in human contact patterns affect the competition dynamics and infection spread

## Contributors
| Team member name |
| ---------------- |
| [Sudarshan Anand](https://github.com/ASudu) |
| Rohini Janivara |
| Alejandro Danies Lopez |


## Setup

We recommend using a Linux/Unix(Mac) environment. For windows users, we recommend using [WSL2](https://learn.microsoft.com/en-us/windows/wsl/install) to create a linux environment.

1. Install `rye` using the following command:
```bash
curl -sSf https://rye.astral.sh/get | bash
source "$HOME/.rye/env"
````
2. `cd` into the project directory and run the following command to create a virtual environment:
```bash
rye init .
```
3. Add dependencies as follows:
```bash
rye add jupyter networkx numpy matplotlib seaborn scipy tqdm sympy pyvis plotly 
```
4. Activate the virtual environment using the following command:
```bash
source .venv/bin/activate
```

## Files

### Synthetic networks generation
- Move into ``Network_generation``
- Run the ``network_gen.py`` to generate the synthetic networks
```bash
cd Network_generation
python network_gen.py
```

### Real world data calibration
- The folder ``Real_world_data`` contains all the data we have tried on for the calibration
- We finally calibrated on Influenza data (``Real_world_data/calibration/Fludata_US_2016_2024.csv``) and on the US High school contact network (``Real_world_data/high_school_nw.txt``)
- The calibration is done in ``calibration.ipynb``

### Simulations
- The ``simulations_and_analysis`` folder has all the notebooks and script files used for the simulation of the infection spread on various contact networks (synthetic and real-world), and the analysis part done for obtaining metrics for epidemic analysis and parameter sensitivity

- At the end of running the simulations, the results are stored in a ``results_final`` folder
