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

