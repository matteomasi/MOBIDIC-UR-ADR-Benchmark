# Software and data for paper: Coupled hydrologic, hydraulic, and surface water quality models for pollution management in urban-rural areas

This repository contains the software and data used in the paper "Coupled hydrologic, hydraulic, and surface water quality models for pollution management in urban-rural areas" by Matteo Masi, Daniele Masseroni, Fabio Castelli.

The repository file structure is organized as follows:
- `adr/`: contains the advection-diffusion-reaction model code.
- `benchmark_01_network_react/`: contains the benchmark of the advection-diffusion-reaction model vs US EPA WASP8 software.
- `benchmark_02_hfcw_1D/`: contains the benchmark of the transport modul of MOBIDIC–UR–CW model vs data from Boog et al. (2019).
- `benchmark_03_cwm1/`: contains the benchmark of reactive module of MOBIDIC–UR–CW model vs PHREEQC software.

To run benchmark simulations, the function `main.m` should be called from MATLAB, located in each benchmark folder. From the MATLAB command line, type:

```matlab
main
```

