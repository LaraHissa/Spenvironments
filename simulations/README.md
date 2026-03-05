# Speciation Dynamics Model

![C++ Version](https://img.shields.io/badge/C++-17-blue.svg)

This repository contains a C++ simulator for an individual-based model (IBM) designed to investigate speciation in sexually reproducing populations. The model is implemented on a discrete lattice with reflective (non-periodic) boundaries, where population dynamics are self-regulated by local and global carrying capacities.

🧬 Model Framework

The simulation explores how individuals adapt and diversify across distinct environments under natural selection and a hard genetic compatibility constraint.

### Implementation Details:

### **node.h**
The `node.h` class represents each individual $i$ within the population, acting as the fundamental structure that stores its specific attributes during the simulation. It maintains the unique `NodeId` and the spatial coordinates $(x, y)$ on the lattice through `NodePosition`. Furthermore, it manages the local connectivity of the individual by storing the indices and positions of its neighbors within the mating radius $S$, providing the necessary data to determine the node's degree and the network of reproductive compatibility used to identify species as connected components.

* **`get_nodeid() / get_position()`**: Returns the individual's unique identifier and its discrete spatial coordinates on the $L_1 \times L_2$ lattice.
* **`degree()`**: Provides the number of neighbors within the mating radius $S$, representing the individual's local connectivity and potential for gene flow.
* **`add_neighbor_position() / add_neighbor_id()`**: Dynamically updates the individual's spatial and topological neighbor lists during the graph reconstruction phase of the simulation.
* **`get_neighbors() / get_neighbors_id()`**: Retrieves the full set of coordinates or IDs of individuals located within the interaction neighborhood $S$.
* **`clear()`**: Resets the individual's attributes, ensuring a clean state for the subsequent generation in the discrete-time cycle of the IBM.



### **graph.h**

The `graph.h` class manages the spatial environment and the positioning of individuals within the lattice. It is responsible for the initial random distribution of nodes and the implementation of the circular neighborhood template based on the mating radius $S$. By iterating through the population, it establishes spatial connectivity through the `add_neighbors` function, identifying all individuals within the specified distance while accounting for reflective boundary conditions. Additionally, it provides essential metrics for the model, such as the calculation of local density and the mean degree of the network, which are necessary for the self-regulation of the population dynamics.

* **`generate_random_positions()`**: Distributes the initial $M_0$ individuals across the lattice using a stochastic process, assigning random $(x, y)$ coordinates to each node while accounting for site occupancy.
* **`region_S()`**: Pre-calculates the circular neighborhood template by identifying all relative coordinates that fall within the Euclidean distance defined by the mating radius $S$.
* **`add_neighbors()`**: Constructs the population network by mapping each individual to its spatial neighbors based on the `region_S` template, respecting the reflective (non-periodic) boundary conditions.
* **`local_density()`**: Measures the population density surrounding a specific site, used to evaluate local carrying capacity and reproductive inhibition based on the threshold $K_{local}$.
* **`mean_degree()`**: Calculates the average connectivity ($\langle k \rangle$) of the graph, serving as a global parameter for characterizing the gene flow potential when the population is at carrying capacity.
* **`clear_node_list() / add_node()`**: Manages the population transition between generations by resetting the current node list and populating the graph with the newly generated offspring.


### **speciation_model.h**

The `speciation_model.h` class serves as the princiapal part of the simulation, orchestrating the evolutionary processes and the genetic architecture of the IBM. It manages the multi-chromosomal system consisting of reproductive ($\gamma$) and phenotypic ($\rho$) binary sequences, while implementing the life cycle dynamics. This includes the calculation of Gaussian fitness based on environmental adaptation, the execution of assortative mating governed by the genetic compatibility threshold $G$, and the spatial dispersal of offspring within the radius $S$. By integrating these mechanisms, the class drives the adaptation of individuals to local optima and the emergence of reproductive isolation through the identification of species as connected components.



* **`dynamics()`**: Orchestrates the full generational cycle, including mating attempts, the transition to the new graph structure, and the update of phenotypic and fitness values for the entire population.
* **`try_reproduce()`**: Manages the probabilistic mating attempts and the compensatory recruitment mechanism ($\eta$), which regulates per-capita fecundity to maintain population stability near the global carrying capacity $M_0$.
* **`reproduce()`**: Implements sexual reproduction with free recombination and applies genetic mutations with probability $\mu$ per locus to generate the genomes of the offspring.
* **`offspring_region()`**: Handles the spatial dispersal of offspring by randomly selecting a site within the reproduction radius $S$ of the parent, respecting the reflective lattice boundaries.
* **`phenotype()`**: Calculates the ecological phenotype $P_i$ through the additive effects of the phenotypic chromosome, providing the quantitative trait required for natural selection.
* **`compatibility_environmental()`**: Computes individual fitness $w_i$ using a Gaussian selection function based on the distance to the local environmental optimum ($P_E$), supporting both homogeneous and heterogeneous niche scenarios.
* **`check_species()`**: Identifies species as the connected components of the population network where edges are defined by the genetic compatibility threshold $G$.
* **`calculate_genetic_neighbors()`**: Determines the set of potential mating partners for each individual by filtering spatial neighbors based on the hard genetic compatibility constraint $D_{ij} < G$.


### **auxiliar_functions.h**

The `auxiliar_functions.h` class provides a specialized diagnostic and data-export interface for the simulation. It acts as the bridge between the numerical simulation and post-simulation analysis, facilitating the systematic recording of the population's evolutionary trajectory. This class handles the generation of structured data files that capture spatial, genomic, and phenotypic distributions, allowing for the characterization of species richness and the quantification of phenotypic divergence across environmental boundaries.

[Image of a data logging and scientific visualization flowchart for an individual-based simulation]

* **`save_diff_phenotype_analysis()`**: Exports the phenotypic differences ($DIFF$) recorded at the environment boundary compared to the global population variance, essential for tracking divergence under heterogeneous selection.
* **`save_node_positions_and_species()`**: Generates spatial snapshots of the lattice, recording each individual's coordinates, species ID, local degree, phenotype, and fitness for distribution analysis.
* **`save_tudo()`**: Performs a comprehensive export of the population network, saving pairs of genetically compatible individuals along with their spatial and phenotypic attributes for network-based gene flow analysis.
* **`save_rep_genomes()`**: Records the complete binary sequences of the reproductive chromosomes ($\gamma$) for the entire population to track genetic drift and molecular divergence.
* **`save_density()`**: Maps the spatial distribution of offspring production across the lattice, identifying reproductive hotspots and the impact of local carrying capacity on recruitment.
* **`save_genetic_neighbors()`**: Exports the adjacency list representing reproductive compatibility, where edges are defined by the genetic threshold $G$.
* **`print_offspring_position() / print_density_offspring_position()`**: Provides real-time console feedback on the spatial dispersal of the new generation during the discrete-time transition.


### **main.cpp**

The `main.cpp` file serves as the primary orchestrator for the speciation model, managing the simulation lifecycle from parameter initialization to the execution of the generational loop. It handles the interface between the user-defined parameters and the internal classes, ensuring that the ecological and genetic constraints are correctly applied. The script is designed to support high-throughput numerical experiments by utilizing command-line arguments and a robust seeding mechanism to guarantee scientific reproducibility.

* **Parameter Handling and Reproducibility**: Parses input arguments such as initial population size ($M_0$), genetic threshold ($G$), and selection width ($\sigma$). It implements a deterministic seeding strategy using `std::mt19937` combined with a hash of the simulation parameters, allowing specific evolutionary trajectories to be reconstructed and audited.
* **Spatial and Biological Initialization**: Orchestrates the setup of the spatial lattice through the `graph` class and initializes the genomic state of the population. This phase includes the pre-calculation of the neighborhood template ($S$), the definition of environmental optima ($P_E$), and the initial calculation of genetic distances and species clusters.
* **Generational Execution Loop**: Executes the discrete-time simulation over $t = 5000$ generations. In each iteration, it invokes the `dynamics` method to process reproduction, mutation, and selection, followed by a re-evaluation of species clusters based on the potential gene flow network.
* **Data Management**: Manages the automated creation of the directory structure (`/All`, `/diff`) and ensures the periodic logging of global population metrics and individual-level snapshots. It maintains the `evolution_global.dat` file, providing a continuous record of species richness and population stability throughout the evolutionary run.



## **Data Outputs**

The simulation generates structured data files for post-simulation analysis, providing a multi-scale view of the evolutionary process. These outputs are organized within the specified directory:

* **`evolution_global[run_id].dat`**: A temporal log updated at each generation $t$. This file records the primary macroscopic variables of the system: the time step, the current species richness (number of connected components), and the dynamic population size $M(t)$.
* **`/All/data[i].dat`**: High-resolution snapshots of the entire population state at generation $i$. Each line represents an individual and contains: `NodeId`, spatial coordinates $x$ and $y$, `SpeciesId`, local degree $k_i$, ecological phenotype $P_i$, and fitness $w_i$.
* **`/diff/diff_boundary_[i].txt`**: Records the phenotypic differences between mating pairs located specifically at the environmental boundary ($L_1/2 \pm S$). This data is critical for analyzing selection against hybrids.
* **`/diff/diff_total_[i].txt`**: A global record of phenotypic differences between all mating pairs across the entire lattice.
* **`seed_final[run_id].dat`**: Stores the complete internal state of the `std::mt19937` pseudo-random number generator at the end of the run to ensure full auditability.

### **Compilation**

To compile the simulator, use the following command with the `-O1` optimization flag:

g++ -O1 main.cpp graph.cpp speciation_model.cpp auxiliar_functions.cpp -o speciation 

Run the executable by passing the required parameters:
./speciation <NumNodes> <G> <sigma> <directory> [run_id]

For additional information, please, contact Lara D. Hissa: l171513@dac.unicamp.br
