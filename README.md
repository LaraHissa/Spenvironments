# Speciation Dynamics Model

This repository contains a C++ simulator for an individual-based model (IBM) used to study speciation in sexually reproducing populations. The model runs on a discrete lattice with reflective boundaries and explores how individuals adapt and diversify across different environments.

## **Project Structure**

The project is organized into two  directories to separate data generation from analysis:

1. **`simulations/`**: Contains the C++ source code to run the evolutionary dynamics of speciation.
2. **`figures/`**: Contains the Python scripts used to process raw data and generate the plots for the paper.

---

## **1. Simulation of Speciation**

The simulator manages the spatial environment, genetic compatibility, and evolutionary dynamics.

* **`node.h`**: Defines individual $i$, storing its `ID`, position $(x, y)$, and the list of spatial neighbors.
* **`graph.h`**: Manages space and connectivity, including random initial positioning and neighborhood calculation via `region_S`.
* **`speciation_model.h`**: Dynamics biological logic, including fitness calculation (Gaussian), mating, mutations, and offspring dispersal.
* **`auxiliar_functions.h`**: Handles data export of phenotypes, fitness, and genomes into `.dat` files.
* **`main.cpp`**: The entry point where parameters (like $M_0, G, \sigma$) are set and the generational loop is executed.

---

## **2. Figures**

Once the simulation is complete, the scripts in the `figures/` folder are used to generate the paper's results:

* **`temporal.ipynb`**: Plots the evolution of species richness and population size over time from `evolution_global.dat`.
* **`histograms.ipynb`**: Generates 4x4 grids showing phenotypic distributions across different selection regimes.
* **`distance.ipynb`**: Figures of IBD and the phenotypic divergence and variance, focusing on environmental boundaries too.


---

## **How to Run**

Compile directly from `main.cpp` using the optimization flag:
```bash
g++ -O1 main.cpp -o speciation
```
Running a Simulation
```Bash 

./speciation <InitialPop> <G> <sigma> <output_directory> [run_id]
```

Data Outputs

The simulator generates the following files required by the Python scripts:

    evolution_global.dat: Records time, species richness, and population size M(t).

    /All/data[t].dat: Population snapshots including ID, position, species, phenotype, and fitness.

    /diff/: Files recording phenotypic differences between mating pairs at boundaries and globally.

Contact: Lara D. Hissa (l171513@dac.unicamp.br)


