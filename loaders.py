import json
from DynamicFBASimulator import DynamicFBASimulator

def load_simulator_from_json(json_path: str, model_name: str = 'textbook'):
    with open(json_path, 'r') as f:
        cfg = json.load(f)
    return DynamicFBASimulator(
        model_name=cfg.get("model_name", model_name),
        dt=cfg.get("dt", 0.01),
        volume=cfg.get("volume", 1.0),
        initial_biomass=cfg.get("biomass_init", 0.1),
        ext_conc=cfg.get("ext_conc", {}),
        vmax_params=cfg.get("vmax_params", {}),
        km_params=cfg.get("km_params", {}),
        kn_params=cfg.get("kn_params", {})
    )
