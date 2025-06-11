import pandas as pd
from summit.domain import ContinuousVariable, Domain
from summit.strategies import TSEMO, LHS
from summit.utils.dataset import DataSet
from reactor_model import simulate_reactor

domain = Domain()
domain += ContinuousVariable("temperature", bounds=[60, 140], description="Temperature (°C)")
domain += ContinuousVariable("conc", bounds=[100, 500], description="Concentration (mol/m3)")
domain += ContinuousVariable("ratio", bounds=[1.0, 5.0], description="Molar Ratio")
domain += ContinuousVariable("res_time", bounds=[0.5, 2.0], description="Residence Time (min)")
domain += ContinuousVariable("STY", is_objective=True, bounds=[0, 1e8], maximize=True, description="Space-Time Yield")
domain += ContinuousVariable("E_factor", is_objective=True, bounds=[0, 50], maximize=False, description="E-Factor")

columns = [v.name for v in domain.variables]
strategy = TSEMO(
    domain=domain
)

df = pd.DataFrame(columns=columns)

# --- Evaluierungsfunktion ---
def evaluate_row(suggestion):
    data = suggestion.to_dict()["data"][0]
    T = round(data[columns.index("temperature")], 4)
    conc = round(data[columns.index("conc")], 4)
    ratio = round(data[columns.index("ratio")], 4)
    res_time = round(data[columns.index("res_time")], 4)

    STY, E = simulate_reactor(T, conc, ratio, res_time)

    return pd.DataFrame([{
        "temperature": T,
        "conc": conc,
        "ratio": ratio,
        "res_time": res_time,
        "STY": STY,
        "E_factor": E,
    }], columns=columns)

# Sample 3 initial points
print("Generating initial samples...")
lhs_sampler = LHS(domain)
initial_suggestions = lhs_sampler.suggest_experiments(3)
for i in range(len(initial_suggestions)):
    suggestion = initial_suggestions.iloc[[i]]
    evaluation = evaluate_row(suggestion)
    print(f"Outcome of initial experiment {i+1}:\n{evaluation}\n")
    df = pd.concat([df, evaluation], ignore_index=True)


# --- Optimierungsschleife ---
N_ITER = 100
for i in range(N_ITER):
    print(f"Iteration {i+1}")
    df_clean = df.dropna().copy()
    dataset = DataSet.from_df(df_clean)


    suggestion = strategy.suggest_experiments(1, prev_res=dataset)


    evaluation = evaluate_row(suggestion)
    print(f"Outcome of experiment:\n{evaluation}\n")

    df = pd.concat([df, evaluation], ignore_index=True)
    df.to_csv("data/tsemo_results.csv", index=False)

print("Optimierung abgeschlossen.")
