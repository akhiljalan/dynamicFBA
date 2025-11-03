from loaders import load_simulator_from_json
from io_utils import save_results_to_csv, plot_results_from_csv_separate

if __name__ == "__main__":
    sim = load_simulator_from_json("configs/ecoli_config1.json")
    n_steps = 500
    sim.run(n_steps=n_steps, verbose=True)

    # TODO generate more descriptive filename 
    save_results_to_csv(sim, "dynamic_fba_results.csv")

    plot_results_from_csv_separate("dynamic_fba_results.csv")
