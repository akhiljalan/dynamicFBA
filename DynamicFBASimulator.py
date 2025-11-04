from cobra.io.web import load_model
from typing import Dict, List, Optional, Tuple 
import time 

class DynamicFBASimulator: 
    '''
    Simulator for dynamic FBA (Flux Balance Analysis). 
    Assumes a single compartment with uniform extracellular concentrations.
    Uses COBRApy metabolic model for stoichiometric information. 

    Units: 
    -- Fluxes (Cobra): mmol / gDW / hr 
    -- Concentrations: mmol / L 
    -- Biomass: gDW 
    -- Volume: L 
    -- dt: seconds 
    '''
    def __init__(
        self, 
        model_name: str = 'textbook',  # Default to textbook E coli 
        biomass_reaction_id: str = "Biomass_Ecoli_core", # Default to biomass reaction in E coli 
        dt: float = 0.01, # Seconds 
        volume: float = 1.0, # L 
        initial_biomass: float = 2.0, # gDw
        vmax_params: Optional[Dict[str, float]] = None, # Reaction ID -> V_max Value (Monod) 
        km_params: Optional[Dict[str, float]] = None, # Reaction ID -> Km Value (Monod)
        kn_params: Optional[Dict[str, float]] = None, # Reaction ID -> Kn Value (inhibitory)
        ext_conc: Optional[Dict[str, float]] = None, # Metabolite ID -> External concentration (mM)
        setpoints: Optional[Dict[str, float]] = None, # Metabolite ID -> External concentration (mM)
        essential_exchanges: Optional[List[str]] = None, 
        clip_negative: bool = True, # If True, clip negative concentrations to zero 
    ): 
        self.model = load_model(model_name)
        self.dt = dt 
        self.dt_hr = self.dt / 3600.0 
        self.volume = volume
        self.initial_biomass = initial_biomass
        self.biomass = self.initial_biomass
        self.vmax_params = vmax_params or {}
        self.km_params = km_params or {}
        self.kn_params = kn_params or {}
        self.ext_conc = ext_conc or {}
        self.setpoints = setpoints or {}
        self._init_missing_metabolites()
        self._check_setpoint_keys()
        self.exchange_reactions_map = self._get_exchange_map()

        self.essential_exchanges = essential_exchanges
        if not self.essential_exchanges: 
            self.essential_exchanges = ["EX_nh4_e","EX_pi_e","EX_h_e","EX_h2o_e",
                    "EX_k_e","EX_na1_e","EX_cl_e","EX_mg2_e","EX_ca2_e","EX_fe2_e"]
        self._apply_essential_bounds()
        self.BIOMASS_RXN = biomass_reaction_id
        self.results = self._init_timeseries()
        self.solution_feasible = True
        self.clip_negative = clip_negative

    def _init_missing_metabolites(self) -> None: 
        '''
        For any exchange reactions whose metabolites 
        are not tracked in ext_conc, set that initial 
        value to zero. 
        '''
        for rxn in self.model.exchanges: 
            for metabolite in rxn.metabolites: 
                if metabolite.id not in self.ext_conc: 
                    self.ext_conc[metabolite.id] = 0.0 

    def _check_setpoint_keys(self) -> None: 
        '''
        For each key in self.setpoints, ensure it is 
        present in self.ext_conc. 
        '''
        for k in self.setpoints: 
            if k not in self.ext_conc: 
                raise ValueError(f"{k} set to {self.setpoints[k]} but not found in external concentrations.")

    def _get_exchange_map(self) -> Dict[str, str]: 
        '''
        Create dictionary whose keys are 
        metabolite IDs, and values are reaction IDs,
        for all exchange reactions. 
        '''
        ex_reactions_map = {}
        for reaction in self.model.exchanges: 
            metabolites = list(reaction.metabolites.keys())
            if len(metabolites) == 1: 
                if metabolites[0].compartment == 'e': 
                    ex_reactions_map[metabolites[0].id] = reaction.id 
        return ex_reactions_map

    def _init_timeseries(self) -> Dict[str, List[float]]: 
        '''
        Initialize results dictionary 
        for tracking time series of simulation results. 
        '''
        return { 
            'time_s': [], 
            'biomass_gDW': [], 
            **{m: [] for m in self.ext_conc.keys()}, # track metabolites per time step 
        }

    def _apply_essential_bounds(self) -> None: 
        '''
        For all essential reactions, set 
        conservative bounds of [-1000.0, 1000.0]. 
        '''
        for rxn_id in self.essential_exchanges:
            try: 
                r = self.model.reactions.get_by_id(rxn_id)
                r.lower_bound = -1000.0
                r.upper_bound = 1000.0
            except KeyError: 
                pass 

    def _set_dynamic_bounds(self) -> float: 
        '''
        1. Update exchange bounds based on external concentrations. 
        The formula is based on Michaelis-Menten kinetics. 
        2. Compute a multiplicative biomass inhibition factor based on 
        concentration of inhibitory products, which are specified in 
        self.kn_params. 
        '''
        # 1. Compute update limits and udate reaction lower bounds accordingly 
        for metabolite_id, reaction_id in self.exchange_reactions_map.items(): 
            if metabolite_id not in self.ext_conc: 
                continue 
            reaction_current = self.model.reactions.get_by_id(reaction_id)
            conc = max(self.ext_conc[metabolite_id], 0.0)
            V_max = self.vmax_params.get(reaction_id, 10.0) # Default Vmax of 10
            Km = self.km_params.get(reaction_id, 0.01) # Default Km of 0.01 

            # Michaelis-Menten uptake lower bound 
            uptake_limit = 0.0 
            if conc > 0: 
                uptake_limit = -1.0 * V_max * conc / (Km + conc)

            reaction_current.lower_bound = uptake_limit 

        # 2. Compute inhibitory factors for biomass creation 
        biomass_inhibition_product_term = 1.0 
        for reaction_id, kn_value in self.kn_params.items(): 
            rxn_current = self.model.reactions.get_by_id(reaction_id)
            for metabolite in rxn_current.metabolites: 
                metabolite_id = metabolite.id 
                if metabolite_id.endswith('_e') and metabolite_id in self.ext_conc: 
                    conc = max(self.ext_conc[metabolite.id], 0.0)
                    if conc > 0: 
                        inhibitory_factor = kn_value / (kn_value + conc)
                        biomass_inhibition_product_term *= inhibitory_factor
        return biomass_inhibition_product_term

    def _update_concentrations(self, fba_soln) -> None: 
        '''
        Euler update for all extracellular concentrations. 

        fba_soln: Solution for FBA, keyed by metabolites. 
        '''
        for metabolite_id, reaction_id in self.exchange_reactions_map.items(): 
            # Only update concentrations for metabolites that are NOT
            # maintained at a setpoint value, AND are present in self.ext_conc. 
            if metabolite_id not in self.ext_conc or metabolite_id in self.setpoints: 
                continue 
            flux = fba_soln.fluxes[reaction_id] 
            # Note that Cobra fluxes are in 1/hr, not 1/s units. 
            delta = flux * self.biomass * self.dt_hr / self.volume 
            self.ext_conc[metabolite_id] += delta 
            # Clip to zero if needed 
            if self.clip_negative and self.ext_conc[metabolite_id] < 0: 
                self.ext_conc[metabolite_id] = 0.0 
            
    def _update_biomass(self, 
        mu: float, 
        inhibition_term: float) -> None: 
        '''
        Euler update for biomass. 
            dX / dt = mu * inhibition_term * X 

        mu: Instantaneous growth rate, computed from FBA solution. 
        inhibition_term: Inhibitory effect of external metabolites.  
        '''
        if not (0.0 <= inhibition_term and inhibition_term <= 1.0): 
            raise ValueError(f'Biomass inhibition term {inhibition_term} should be in [0, 1].')
        # NOTE: Biomass is in gDW, no need to normalize by volume. 
        self.biomass += mu * inhibition_term * self.biomass * self.dt_hr 

    def step(self, t: int) -> Tuple[float, float]:
        '''
        Compute single timestep of dynamic FBA. 
        ''' 
        # get inhibitory effects, apply dynamic exchange bounds 
        biomass_inhibition_product_term = self._set_dynamic_bounds()

        # solve for the model 
        sol = self.model.optimize()
        if sol.status.lower() != 'optimal': 
            print(f'Solver status at time {t} is {sol.status}. Expected "optimal".')
            self.solution_feasible = False

        mu = float(sol.fluxes[self.BIOMASS_RXN]) # In 1/hr units 

        self._update_concentrations(sol)
        self._update_biomass(mu, biomass_inhibition_product_term)

        # Record results 
        self.results['time_s'].append(self.dt * t)
        self.results['biomass_gDW'].append(self.biomass)
        for m in self.ext_conc: 
            self.results[m].append(self.ext_conc[m])

        return mu, biomass_inhibition_product_term
    
    def run(self, n_steps: int, verbose: bool = False) -> None:
        '''
        Run dynamic FBA simulation for n_steps steps. 
        ''' 
        # If verbose, then print every 10% of the steps. 
        print_every = int(n_steps / 10)
        start_time = time.time()

        if verbose: 
            print(f'Running dynamic FBA for {n_steps} steps with dt={self.dt}.')
        for timestep in range(n_steps): 
            if self.solution_feasible: 
                mu, biomass_inhibition_product_term = self.step(timestep)
                if verbose and timestep % print_every == 0: 
                    num_seconds_simulated = timestep * self.dt
                    num_hours_simulated = num_seconds_simulated / 3600.0
                    print(f'{num_hours_simulated:.4f} hours simulated.')
                    print(f'Biomass: {self.biomass:.4f} gDW.')
            else: 
                # Halt simulation when infeasible 
                break 
        if verbose: 
            cur_time = time.time()
            print(f'Final Biomass: {self.biomass:.4f} gDW.')
            print(f'Dynamic FBA simulation completed in {cur_time - start_time:.2f} seconds.')
