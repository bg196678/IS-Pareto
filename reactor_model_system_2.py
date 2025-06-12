import pandas as pd
import numpy as np
from pyomo.environ import *
from pyomo.dae import *

kinetics_dataframe = pd.read_excel("data/system_2/kinetics.xlsx")
gsolv_dataframe = pd.read_excel("data/system_2/gsolv.xlsx")
gsolv_dataframe = gsolv_dataframe.rename(columns={
    'HBr': 'LeavingGroup',
    "Reactant1": 'Reactant2',
    "Reactant2": 'Reactant1',
})

R = 8.3145  # J/mol/K

species = [
    'Reactant1',
    'Reactant2',
    'TS1',
    'TS2',
    'Product1',
    'Product2',
    'LeavingGroup'
]

reactions = [
    'R1',
    'R2',
    'R3',
    'R4'
]

stoich = {
    'R1': {'Reactant1': -1, 'Reactant2': -1, 'TS1': 1},
    'R2': {'Reactant2': -1, 'Product1': -1, 'TS2': 1},
    'R3': {'TS1': -1, 'Product1': 1, 'LeavingGroup': 1},
    'R4': {'TS2': -1, 'Product2': 1, 'LeavingGroup': 1},
}

reaction_TS_map = {
    'R1': 'TS1',
    'R2': 'TS2',
    'R3': 'TS1',
    'R4': 'TS2',
}

reaction_reactants_map = {
    'R1': ['Reactant1', 'Reactant2'],
    'R2': ['Product1', 'Reactant2'],
    'R3': ['TS1'],
    'R4': ['TS2'],
}

molar_masses = {
    'Reactant1': 0.159,
    'Reactant2': 0.071,
    'TS1': 0.230,
    'TS2': 0.230,
    'Product1': 0.21008046,
    'Product2': 0.21008046,
    'LeavingGroup': 0.020,
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

def simulate_reactor(T_C, conc_reactant1, ratio, t_res_min):
    T = T_C + 273.15
    t_res = t_res_min * 60
    conc_reactant2 = conc_reactant1 * ratio

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
        m.C["Reactant1", 0].fix(conc_reactant1)
        m.C["Reactant2", 0].fix(conc_reactant2)
        for sp in species:
            if sp not in ["Reactant1", "Reactant2"]:
                m.C[sp, 0].fix(0.0)

    model.init = BuildAction(rule=init_conds)

    TransformationFactory(
        "dae.finite_difference"
    ).apply_to(model, nfe=200, scheme="BACKWARD")
    solver = SolverFactory("ipopt")
    _ = solver.solve(model, tee=False)

    c_raw = {sp: value(model.C[sp, t_res]) for sp in species}
    mass_in = (
            conc_reactant1 * molar_masses["Reactant1"] +
           conc_reactant2 * molar_masses["Reactant2"]
    )
    mass_sim = sum(c_raw[sp] * molar_masses[sp] for sp in species)

    norm_factor = mass_in / mass_sim

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
