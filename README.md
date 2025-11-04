# Dynamic Flux Balance Analysis 

This repository contains a lightweight and modular implementation of **Dynamic Flux Balance Analysis** (Dynamic FBA). It explicitly handles: 
- Biomass growth 
- Time-varying extracellular concentrations 
- Optional biomass inhibition 
- Flexible modeling via COBRApy's metabolic modeling backend 
- Built-in logging to CSV and plotting utilities 

## Installation 

Fork/clone and install. 
```bash
git clone https://github.com/YOUR_USERNAME/dynamicFBA.git
cd dynamicFBA
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
3. Uses Michaelis-Menten formulas to update flux bounds for uptake reactions based on extracellular concentrations. For metabolite $m$ with concentration $C(m)$ and parameters $V_{max}(m), K_m(m)$, the flux lower bound $\ell_M$ becomes: 

$$\ell_M = -V_{max}(m) \frac{K_m(m)}{K_m(m) + C(m)}$$

4. Solves a FBA instance (via Linear Programming) to obtain fluxes and biomass growth rate. Note that flux lower bounds come from step 3. 
5. Optionally computes inhibitory term for biomass growth due to external product (e.g. Acetate and Lactate). If $M$ is the set of metabolites that are inhibitory, where $m \in M$ has concentration $C(m)$ and parameter $K_n(m)$, then the biomass inhibition term $\nu$ becomes: 

$$\nu = \prod\limits_{m \in M} \frac{K_n(m)}{K_{n}(m) + C(m)}$$

6. Updates biomass and extracellular concentrations based on Forward Euler stepping. 

The simulation halts when the step limit is reached, or when the FBA instance in Step (4) is infeasible. 

## Demo: E Coli in Batch Mode with DO Control  

This demo is specified in `configs/ecoli_config1.json`. Note that the media composition is taken from the following protocols: 
- [DSMZ Medium 382](https://www.dsmz.de/microorganisms/medium/pdf/DSMZ_Medium382.pdf)
- [Urniezius et al 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6833213/)

To run the demo, 
```bash
cd dynamicFBA/
python demo_ecoli1.py
```

In this demo, we control DO (dissolved oxygen). Typically setpoint control is done via PID control; for purposes of simulation, we assume that the lag time is negligible and therefore simply clamp setpoint values throughout the simulation. 

After running the demo script, results are saved to the `demo_results/demo_ecoli1/` directory. Your demo should produce plots in a file called `demo_results/demo_ecoli1/results_grid.pdf` that look like this. 

![image](figs/results_grid_ecoli1.jpg)
