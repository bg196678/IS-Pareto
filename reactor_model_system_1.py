import pandas as pd
import numpy as np
from pyomo.environ import *
from pyomo.dae import *

kinetics_dataframe = pd.read_excel("data/system_1/kinetics.xlsx")
gsolv_dataframe = pd.read_excel("data/system_1/gsolv.xlsx")

R = 8.3145  # J/mol/K

species = [
    'Substrate',
    'Nucleophilic',
    'ITS1',
    'ITS2',
    'ITS3',
    'ITS4',
    'Product1',
    'Product2',
    'Product3',
    'LeavingGroup'
]

molar_masses = {
    "Substrate": 0.159,
    "Nucleophilic": 0.071,
    "ITS1": 0.159 + 0.071,
    "ITS2": 0.159 + 0.071,
    "ITS3": 0.159 + 0.071 + 0.071,
    "ITS4": 0.159 + 0.071 + 0.071,
    "Product1": 0.21008046,
    "Product2": 0.21008046,
    "Product3": 0.26114773,
    "LeavingGroup": 0.020,
}

reactions = [
    'R1_fwd',
    'R1_rev',
    'R2_fwd',
    'R2_rev',
    'R3_fwd',
    'R3_rev',
    'R4_fwd',
    'R4_rev',
    'R5_fwd',
    'R5_rev',
    'R6_fwd',
    'R6_rev',
    'R7_fwd',
    'R7_rev',
    'R8_fwd',
    'R8_rev'
]

stoich = {
    'R1_fwd': {'Substrate': -1, 'Nucleophilic': -1, 'ITS1': 1},
    'R1_rev': {'ITS1': -1, 'Substrate': 1, 'Nucleophilic': 1},
    'R2_fwd': {'Substrate': -1, 'Nucleophilic': -1, 'ITS2': 1},
    'R2_rev': {'ITS2': -1, 'Substrate': 1, 'Nucleophilic': 1},
    'R3_fwd': {'Product2': -1, 'Nucleophilic': -1, 'ITS3': 1},
    'R3_rev': {'ITS3': -1, 'Product2': 1, 'Nucleophilic': 1},
    'R4_fwd': {'Product1': -1, 'Nucleophilic': -1, 'ITS4': 1},
    'R4_rev': {'ITS4': -1, 'Product1': 1, 'Nucleophilic': 1},
    'R5_fwd': {'ITS1': -1, 'Product1': 1, 'LeavingGroup': 1},
    'R5_rev': {'Product1': -1, 'LeavingGroup': -1, 'ITS1': 1},
    'R6_fwd': {'ITS2': -1, 'Product2': 1, 'LeavingGroup': 1},
    'R6_rev': {'Product2': -1, 'LeavingGroup': -1, 'ITS2': 1},
    'R7_fwd': {'ITS3': -1, 'Product3': 1, 'LeavingGroup': 1},
    'R7_rev': {'Product3': -1, 'LeavingGroup': -1, 'ITS3': 1},
    'R8_fwd': {'ITS4': -1, 'Product3': 1, 'LeavingGroup': 1},
    'R8_rev': {'Product3': -1, 'LeavingGroup': -1, 'ITS1': 1},
}

reaction_TS_map = {
    'R1_fwd': 'TS1_fwd',
    'R1_rev': 'TS1_rev',
    'R2_fwd': 'TS2_fwd',
    'R2_rev': 'TS2_rev',
    'R3_fwd': 'TS3_fwd',
    'R3_rev': 'TS3_rev',
    'R4_fwd': 'TS4_fwd',
    'R4_rev': 'TS4_rev',
    'R5_fwd': 'TS12_fwd',
    'R5_rev': 'TS12_rev',
    'R6_fwd': 'TS22_fwd',
    'R6_rev': 'TS22_rev',
    'R7_fwd': 'TS32_fwd',
    'R7_rev': 'TS32_rev',
    'R8_fwd': 'TS42_fwd',
    'R8_rev': 'TS42_rev',
}

reaction_reactants_map = {
    'R1_fwd': ['Substrate', 'Nucleophilic'],
    'R1_rev': ['ITS1'],
    'R2_fwd': ['Substrate', 'Nucleophilic'],
    'R2_rev': ['ITS2'],
    'R3_fwd': ['Product2', 'Nucleophilic'],
    'R3_rev': ['ITS3'],
    'R4_fwd': ['Product1', 'Nucleophilic'],
    'R4_rev': ['ITS4'],
    'R5_fwd': ['ITS1'],
    'R5_rev': ['Product1', 'LeavingGroup'],
    'R6_fwd': ['ITS2'],
    'R6_rev': ['Product2', 'LeavingGroup'],
    'R7_fwd': ['ITS3'],
    'R7_rev': ['Product3', 'LeavingGroup'],
    'R8_fwd': ['ITS4'],
    'R8_rev': ['Product3', 'LeavingGroup'],
}

molar_masses = {
    "Substrate": 0.159,
    "Nucleophilic": 0.071,
    "ITS1": 0.159 + 0.071,
    "ITS2": 0.159 + 0.071,
    "ITS3": 0.159 + 0.071 + 0.071,
    "ITS4": 0.159 + 0.071 + 0.071,
    "Product1": 0.21008046,
    "Product2": 0.21008046,
    "Product3": 0.26114773,
    "LeavingGroup": 0.020,
}

def interpolate_value(
        df: pd.DataFrame,
        x_col: str,
        y_col: str,
        x_val: float,
) -> float:
    """Interpolates the grid values in the excel sheets"""
    df_sorted = df.sort_values(by=x_col)
    return float(np.interp(x_val, df_sorted[x_col], df_sorted[y_col]))

def get_k_gas(T: float, reaction: str) -> float:
    """Returns the gas phase rate constant for a reaction at temperature T"""
    return interpolate_value(
        kinetics_dataframe, "Temperature", reaction_TS_map[reaction], T
    )

def get_gsolv(T: float, species_name: str) -> float:
    return interpolate_value(
        gsolv_dataframe, "Temperature (K)", species_name, T
    )

def get_correction_factor(T, reaction):
    ts = reaction_TS_map[reaction]
    deltaG = (
                 get_gsolv(T, ts) - sum(
                 get_gsolv(T, r)
                 for r in reaction_reactants_map[reaction]
             )
             ) * 4184
    return np.exp(-deltaG / (R * T))

def simulate_reactor(T_C, conc_substrate, ratio, t_res_min):
    T = T_C + 273.15
    t_res = t_res_min * 60
    conc_nucleophilic = conc_substrate * ratio

    k_values = {
        r: get_k_gas(T, r) * get_correction_factor(T, r) for r in reactions
    }

    model = ConcreteModel()
    model.t = ContinuousSet(bounds=(0, t_res))
    model.C = Var(species, model.t, domain=NonNegativeReals)
    model.dCdt = DerivativeVar(model.C, wrt=model.t)

    def rate_rule(j, t):
        k = k_values[j]
        rate = k
        for sp in reaction_reactants_map[j]:
            rate *= model.C[sp, t]
        return rate

    def mass_balance_rule(m, i, t):
        return m.dCdt[i, t] == sum(
            stoich[j].get(i, 0) * rate_rule(j, t) for j in reactions
        )

    model.mass_bal = Constraint(
        species, model.t, rule=mass_balance_rule
    )

    def init_conds(m):
        m.C["Substrate", 0].fix(conc_substrate)
        m.C["Nucleophilic", 0].fix(conc_nucleophilic)
        for sp in species:
            if sp not in ["Substrate", "Nucleophilic"]:
                m.C[sp, 0].fix(0.0)

    model.init = BuildAction(rule=init_conds)

    TransformationFactory(
        "dae.finite_difference"
    ).apply_to(model, nfe=200, scheme="BACKWARD")
    solver = SolverFactory("ipopt")
    _ = solver.solve(model, tee=False)

    c_raw = {sp: value(model.C[sp, t_res]) for sp in species}
    mass_in = (
            conc_substrate * molar_masses["Substrate"] +
           conc_nucleophilic * molar_masses["Nucleophilic"]
    )
    mass_sim = sum(c_raw[sp] * molar_masses[sp] for sp in species)

    norm_factor = min(1.0, mass_in / mass_sim)

    c_end_normed = {sp: c_raw[sp] * norm_factor for sp in species}
    c_end = {sp: c_raw[sp] for sp in species}

    V = 1.0
    mass_product_normed = c_end_normed["Product1"] * molar_masses["Product1"]
    mass_waste_normed_total = sum(
        [c_end[sp] * molar_masses[sp] for sp in species if sp != "Product1"]
    )

    STY = 3600 * mass_product_normed / (V * t_res)
    E_factor = mass_waste_normed_total / mass_product_normed

    return STY, E_factor
