from loaders import load_simulator_from_json
from io_utils import save_results_to_csv, plot_results_from_csv_separate

if __name__ == "__main__":
    sim = load_simulator_from_json("configs/ecoli_config1.json")
    dt_seconds = sim.dt 

    # Set n_hours here for the simulation. 
    n_hours = 2.0 
    n_steps = int(3600.0 * n_hours / dt_seconds)
    sim.run(n_steps=n_steps, verbose=True)

    # TODO generate more descriptive filename 
    save_results_to_csv(sim, "demo_ecoli1.csv")

    plot_results_from_csv_separate("demo_ecoli1.csv")
