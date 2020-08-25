# Numeric Simulation Laboratory

Exercise repository for Numeric Simulation Laboratory course, at Deparment of Physics, University of Milan.

## :blue_book: Description

The repo is structured as follow:

* Exercise number folder
  * `cpp/`
     * `main.cpp` the main simulation program that include required classes and external functions
     * `Makefile` ready to compile compile `main.cpp` with its classes and external functions
     * `func.(h|cpp)` appo functions
  * `out/` results folder, containing all the output data of cpp simulation needed for reports
  * `python/`
     * `ipynb` report of the exercise with graphs and parameters used in the exercise simulation
     * `default.py` custom style and default import packages inside `jupyter-notebook`
* Classes build to solve easily and clearly certain exercices, permitting code re-use
   * **RC | Random Class**, the main class, containing a lot of methods (the most important **blockingMethod**) for different types of exercises and minor classes
   * **MD | Molecolar Dynamic Class**, class for deterministic molecolar dynamic simulation, used for **Ex. 04** and **Ex. 07**
   * **I1D | Ising 1 Dimension Class**, class for Ising 1D Model simulation, used for **Ex. 06**
   * **MCMD | Monte Carlo Molecolar Dynamic Class**, class for dynamic molecolar dynamic simulation using Monte Carlo methods, partially derived from **MD** class, used for **Ex. 07**
   * **QMC_1D | Quantum Monte Carlo 1D**, not a real class, is the code provided for PIGS/PIMS, used for **Ex. 08**
   * **TSP | Traveling Salesman Problem Class**, class to solve **TSP** optimization problem, this class contain both **Genetic algo** and **Simulated Annealing algo** for the resolution of the problem. Used for **Ex. 09** and **Ex. 10**
* `custom.css` global custom template for jupyter-notebook
* `theme.mplstyle` global custom template for matplotlib
   
   **Exceptions:**
   
   * Some **class** folders and `cpp/` folders, can contains some input file for parameters, for example **RC** reads `Primes` and `seed.in`
   * In the **Ex. 04** I preferred to save the configuration position **r** and **r-old** inside `cpp/` instead `out/` because are more input parameters that output results needed for python report
   * In the **Ex. 06** there is a `theory.py` inside `python/` that is only a trivial file that the jupyter import, containing the theoretical curves used for checking the simulation results
   * In the **Ex. 07** I preferred to save the configuration position **r** inside `cpp/` instead `out/` because are more input parameters that output results needed for python report
   * In the **Ex. 07** there is a `correlation.cpp` source code inside `cpp/` that I used to compute auto-correlation function starting from raw `.out` data, there is no compilation rule inside `Makefile` for it, must be compiled manually
   * In the **Ex. 07** there is a `data-blocking-sm.py` inside `python/` that I used to compute multiple data-blocking starting from raw `.out` data
   * In the **Ex. 08** there is a `theory.py` inside `python/` that is only a trivial file that the jupyter import, containing the theoretical curves used for checking the simulation results
   * In the **Ex. 08** there is a `get-intg.py` inside `python/` that I used to normalize the eigenfunction of the exercise, and print out very long mathematical expressions as cpp code and/or python code, ready to be copyed and pasted
   * In the **Ex. 10** there is the standard `cpp/` and the `cpp-mpi/` dir containing the code used for the parallelization (MPI). Because I used a different implementation of **TSP** class (files `tsp-mpi.(h|cpp)`),\
I divided the `main.cpp` into two pieces that must be compiled using separated `Makefile`
   * In the **Ex. 11** `cpp/` and `out/` are empty, cause all exercise (simulation+report) is inside jupyter-notebook
   * In the **Ex. 12** `cpp/` and `out/` are empty, cause all exercise (simulation+report) is inside jupyter-notebook
   
## :star: Features

* As much as possible I used **stdlib** objects and methods (not always optimized)
* As much as possible I used a global `blockingMethod` pub function inside **RC** class, re-using code for different exercises.\
The main idea is to define a **lambda** function in the `main.cpp` of the exercise, that can be passed as `function` object (thanks to `functional` std library) to the `blockingMethod`, able to compute what needed at each simulation step.
* Custom style for jupyter-notebooks and matplotlib plots
  ```python
  from default import *
  ```
  > **Warning:** the GITHUB preview of jupyter-notebooks NOT RENDER the custom styles, that can only be visible if loading the notebooks using python kernel and run the 1st input cell
