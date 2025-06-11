import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('data/tsemo_results.csv')

is_pareto = np.ones(df.shape[0], dtype=bool)
for i, row in df.iterrows():
    for j, other_row in df.iterrows():
        if i == j:
            continue
        # Check if other_row dominates row
        if (
                other_row['STY'] >= row['STY'] and
                other_row['E_factor'] <= row['E_factor']) and \
           (other_row['STY'] > row['STY'] or
            other_row['E_factor'] < row['E_factor']
           ):
            is_pareto[i] = False
            break

pareto_df = df[is_pareto].sort_values(by='E_factor')


plt.figure(figsize=(10, 6))
plt.scatter(
    df['STY'],
    df['E_factor'],
    c='blue', label='All Points', alpha=0.5
)
plt.scatter(
    pareto_df['STY'],
    pareto_df['E_factor'],
    c='red', label='Pareto Front', s=100, edgecolors='black'
)
plt.plot(
    pareto_df['STY'],
    pareto_df['E_factor']
    , c='red', linestyle='-', marker=''
)
plt.xlabel('STY')
plt.ylabel('E_factor')
plt.title('Pareto Front: E_factor vs. STY')
plt.legend()
plt.grid(True)
plt.show()