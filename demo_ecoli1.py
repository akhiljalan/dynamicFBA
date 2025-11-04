from loaders import load_simulator_from_json
import io_utils

if __name__ == "__main__":
    sim = load_simulator_from_json("configs/ecoli_config1.json")
    dt_seconds = sim.dt 

    # Set n_hours here for the simulation. 
    n_hours = 2.0 
    n_steps = int(3600.0 * n_hours / dt_seconds)
    sim.run(n_steps=n_steps, verbose=True)

    demo_folder_path = 'demo_results/demo_ecoli1'
    io_utils.save_results_to_csv(sim, f"{demo_folder_path}/results_timeseries.csv")
    io_utils.plot_results_to_pdf_grid(f"{demo_folder_path}/results_timeseries.csv", 
        f"{demo_folder_path}/results_grid.pdf")
