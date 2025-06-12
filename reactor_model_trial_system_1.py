import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from reactor_model_system_1 import simulate_reactor


experimental_data = pd.read_excel("data/system_1/experimental.xlsx")

# Simulated Reactor Data
E_reactor = []
STY_reactor = []

# Experimental Data
E_experimental = []
STY_experimental = []

# Conditions
temps = []
ratios = []
concs = []

for index, row in experimental_data.iterrows():

    tres_min = row["tres/min"]
    ratio = row["2:1"]
    concentration = row["Conc 1/M"] * 1000
    temperature = row["Temp/°C"]
    E_ex = row["E-factor"]
    STY_ex = row["STY/kg m-3 h-1"]

    print(f"Time Minutes: {tres_min}")
    print(f"Ratio: {ratio}")
    print(f"Concentration: {concentration}")
    print(f"Temperature: {temperature}")

    result = simulate_reactor(
        T_C=temperature,
        conc_substrate=concentration,
        ratio=ratio,
        t_res_min=tres_min,
    )

    STY, E = result

    print(f"E experimental: {E_ex}")
    print(f"E Reactor: {E}")
    print(f"E abs deviation = {np.abs(E - E_ex):.4f}")
    print(f"STY experimental: {STY_ex}")
    print(f"STY Reactor: {STY}")
    print(f"STY abs deviation = {np.abs(STY - STY_ex):.4f}")
    print("="*40)

    E_reactor.append(E)
    STY_reactor.append(STY)
    E_experimental.append(E_ex)
    STY_experimental.append(STY_ex)
    temps.append(temperature)
    ratios.append(ratio)
    concs.append(concentration)

E_reactor = np.array(E_reactor)
STY_reactor = np.array(STY_reactor)
E_experimental = np.array(E_experimental)
STY_experimental = np.array(STY_experimental)

E_abs_errors = np.abs(E_reactor - E_experimental)
E_rel_errors = np.abs(E_reactor - E_experimental) / E_experimental * 100
STY_abs_errors = np.abs(STY_reactor - STY_experimental)
STY_rel_errors = np.abs(
    STY_reactor - STY_experimental
) / STY_experimental * 100

E_correlation = np.corrcoef(E_reactor, E_experimental)[0, 1]
STY_correlation = np.corrcoef(STY_reactor, STY_experimental)[0, 1]
E_r2 = E_correlation**2
STY_r2 = STY_correlation**2

print("\n" + "="*50)
print("ERROR STATISTICS")
print("="*50)

print("\nE-FACTOR ERRORS:")
print(f"Mean Absolute Error: {np.mean(E_abs_errors):.4f}")
print(f"Std Dev of Absolute Error: {np.std(E_abs_errors):.4f}")
print(f"RMSE: {np.sqrt(np.mean((E_reactor - E_experimental)**2)):.4f}")
print(f"Mean Relative Error: {np.mean(E_rel_errors):.2f}%")
print(f"Max Absolute Error: {np.max(E_abs_errors):.4f}")
print(f"Min Absolute Error: {np.min(E_abs_errors):.4f}")
print(f"R²: {E_r2:.4f}")

print("\nSTY ERRORS:")
print(f"Mean Absolute Error: {np.mean(STY_abs_errors):.4f}")
print(f"Std Dev of Absolute Error: {np.std(STY_abs_errors):.4f}")
print(f"RMSE: {np.sqrt(np.mean((STY_reactor - STY_experimental)**2)):.4f}")
print(f"Mean Relative Error: {np.mean(STY_rel_errors):.2f}%")
print(f"Max Absolute Error: {np.max(STY_abs_errors):.4f}")
print(f"Min Absolute Error: {np.min(STY_abs_errors):.4f}")
print(f"R²: {STY_r2:.4f}")



fig, axes = plt.subplots(2, 3, figsize=(15, 10))

axes[0, 0].scatter(E_experimental, E_reactor, alpha=0.6)
axes[0, 0].plot([E_experimental.min(), E_experimental.max()],
                [E_experimental.min(), E_experimental.max()], 'r--')
axes[0, 0].set_xlabel('Experimental E-factor')
axes[0, 0].set_ylabel('Model E-factor')
axes[0, 0].set_title('E-factor: Model vs Experimental')

axes[0, 1].hist(E_rel_errors, bins=20, edgecolor='black')
axes[0, 1].set_xlabel('Relative Error (%)')
axes[0, 1].set_ylabel('Frequency')
axes[0, 1].set_title('E-factor Relative Error Distribution')

axes[0, 2].boxplot([E_abs_errors], labels=['E-factor'])
axes[0, 2].set_ylabel('Absolute Error')
axes[0, 2].set_title('E-factor Absolute Error Box Plot')

axes[1, 0].scatter(STY_experimental, STY_reactor, alpha=0.6)
axes[1, 0].plot([STY_experimental.min(), STY_experimental.max()],
                [STY_experimental.min(), STY_experimental.max()], 'r--')
axes[1, 0].set_xlabel('Experimental STY')
axes[1, 0].set_ylabel('Model STY')
axes[1, 0].set_title('STY: Model vs Experimental')

axes[1, 1].hist(STY_rel_errors, bins=20, edgecolor='black')
axes[1, 1].set_xlabel('Relative Error (%)')
axes[1, 1].set_ylabel('Frequency')
axes[1, 1].set_title('STY Relative Error Distribution')

axes[1, 2].boxplot([STY_abs_errors], labels=['STY'])
axes[1, 2].set_ylabel('Absolute Error')
axes[1, 2].set_title('STY Absolute Error Box Plot')
plt.savefig("plots/system_1/reactor_model_analysis.png")
plt.tight_layout()
plt.show()


df_sim = pd.DataFrame({
    "Temp": temps,
    "Ratio_2_1": ratios,
    "Conc": concs,
    "STY": STY_reactor,
    "E_factor": E_reactor
})

df_exp = pd.DataFrame({
    "Temp": temps,
    "Ratio_2_1": ratios,
    "Conc": concs,
    "STY": STY_experimental,
    "E_factor": E_experimental
})

params = ["Temp", "Ratio_2_1", "Conc"]

corr_exp = df_exp[params + ["STY", "E_factor"]].corr()
corr_sim = df_sim[params + ["STY", "E_factor"]].corr()

sty_corr = pd.DataFrame({
    "Experiment": corr_exp["STY"].drop("STY"),
    "Simulation": corr_sim["STY"].drop("STY")
})

e_corr = pd.DataFrame({
    "Experiment": corr_exp["E_factor"].drop("E_factor"),
    "Simulation": corr_sim["E_factor"].drop("E_factor")
})

plt.figure(figsize=(8, 5))
sty_corr.plot(kind="bar", ax=plt.gca())
plt.title("Korrelation der Parameter mit STY")
plt.ylabel("Korrelationskoeffizient")
plt.axhline(0, color="black", linewidth=0.8)
plt.xticks(rotation=0)
plt.tight_layout()
plt.savefig("plots/system_1/sty_correlation.png")
plt.show()

plt.figure(figsize=(8, 5))
e_corr.plot(kind="bar", ax=plt.gca())
plt.title("Korrelation der Parameter mit E-Factor")
plt.ylabel("Korrelationskoeffizient")
plt.axhline(0, color="black", linewidth=0.8)
plt.xticks(rotation=0)
plt.tight_layout()
plt.savefig("plots/system_1/e_factor_correlation.png")
plt.show()
