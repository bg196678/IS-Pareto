import pandas as pd
from summit.domain import ContinuousVariable, Domain
from summit.strategies import TSEMO, LHS
from summit.utils.dataset import DataSet
from reactor_model import simulate_reactor

###################
# Domain Definition
###################

domain = Domain()
domain += ContinuousVariable(
    "temperature",
    bounds=[60, 140],
    description="Temperature (°C)"
)
domain += ContinuousVariable(
    "conc",
    bounds=[100, 500],
    description="Concentration (mol/m3)"
)
domain += ContinuousVariable(
    "ratio",
    bounds=[1.0, 5.0],
    description="Molar Ratio"
)
domain += ContinuousVariable(
    "res_time",
    bounds=[0.5, 2.0],
    description="Residence Time (min)"
)
domain += ContinuousVariable(
    "STY",
    is_objective=True,
    bounds=[0, 1e8],
    maximize=True,
    description="Space-Time Yield"
)
domain += ContinuousVariable(
    "E_factor",
    is_objective=True,
    bounds=[0, 50],
    maximize=False,
    description="E-Factor"
)

columns = [v.name for v in domain.variables]

results_dataframe = pd.DataFrame(columns=columns)

# --- Evaluierungsfunktion ---
def reactor(suggestion: DataSet) -> pd.DataFrame:
    data = suggestion.to_dict()["data"][0]
    T = round(data[columns.index("temperature")], 4)
    conc = round(data[columns.index("conc")], 4)
    ratio = round(data[columns.index("ratio")], 4)
    res_time = round(data[columns.index("res_time")], 4)

    STY, E = simulate_reactor(T, conc, ratio, res_time)

    return pd.DataFrame(
        [{
                "temperature": T,
                "conc": conc,
                "ratio": ratio,
                "res_time": res_time,
                "STY": STY,
                "E_factor": E,
            }],
        columns=columns
    )

#######################
# Initial Sampling
#######################

num_initial_samples = 3
lhs_sampler = LHS(domain)
initial_suggestions = lhs_sampler.suggest_experiments(num_initial_samples)

for i in range(num_initial_samples):
    suggestion = initial_suggestions.iloc[[i]]
    evaluation = reactor(suggestion)
    print(f"Outcome of initial experiment {i+1}:\n{evaluation}\n")
    results_dataframe = pd.concat(
        [results_dataframe, evaluation], ignore_index=True
    )


###################
# Optimization Strategy
###################

num_samples = 100
tsemo_strategy = TSEMO(
    domain=domain
)

for i in range(num_samples):

    dataset = DataSet.from_df(results_dataframe)
    suggestion = tsemo_strategy.suggest_experiments(
        1, prev_res=dataset
    )
    evaluation = reactor(suggestion)

    print(f"Outcome of experiment:\n{evaluation}\n")

    results_dataframe = pd.concat(
        [results_dataframe, evaluation], ignore_index=True
    )
    results_dataframe.to_csv("data/tsemo_results.csv", index=False)
