# Dynamic Flux Balance Analysis 

This repository contains a lightweight and modular implementation of **Dynamic Flux Balance Analysis** (Dynamic FBA). It explicitly couples: 
- Biomass growth 
- Time-varying extracellular concentrations 
- Optional biomass inhibition 
- Flexible modeling via COBRApy's metabolic modeling backend 
- Built-in logging to CSV and plotting utilities 

## Installation 

Fork/clone and install. 
```bash
git clone https://github.com/YOUR_USERNAME/dynamicFBA.git
cd dynamic-fba
```

## Python Requirements 

We recommend Python 3.10+ and pip for package management. 

```bash
pip install -r requirements.txt
```

## Core Logic 

The core logic of the simulation is in the `DynamicFBASimulator.py` file. In a nutshell, this class does the following. 

1. Initializes DynamicFBASimulator object with information from a config file (e.g. `configs/ecoli_config1.json`). The config file specifies Michaelis-Menten constants, the name of the COBRA model (e.g. `textbook`) to be used, and other simulation parameters. 
2. Stores biomass and extracellular concentrations for exchange reactions. 
3. Uses Michaelis-Menten formulas to update flux bounds for uptake reactions based on extracellular concentrations. 
4. Solves a FBA instance (via Linear Programming) to obtain fluxes and biomass growth rate. 
5. Optionally computes inhibitory term for biomass growth due to external product (e.g. Acetate and Lactate). 
6. Updates biomass and extracellular concentrations based on Euler-Maruyama stepping. 

The simulation halts when the step limit is reached, OR when the FBA instance in Step (4) is infeasible. 

## Demos 

### Demo 1: E Coli Core in Batch Mode

This demo is specified in `configs/ecoli_config1.json`. Note that the media composition is taken from the following protocols: 
- [DSMZ Medium 382](https://www.dsmz.de/microorganisms/medium/pdf/DSMZ_Medium382.pdf)
- [Urniezius et al 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6833213/)

To run the demo, 
```bash
cd dynamicFBA/
python demo_ecoli1.py
```

Results will be plotted to `dynamicFBA/demo_ecoli1.csv`. 

