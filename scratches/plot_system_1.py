import imageio
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



def pareto_front(df: pd.DataFrame) -> pd.DataFrame:
    """Calculate the Pareto front from a DataFrame with 'STY' and 'E_factor'
    columns.
    """
    is_pareto = np.ones(df.shape[0], dtype=bool)
    for i, row in df.iterrows():
        for j, other_row in df.iterrows():
            if i == j:
                continue
            if (
                    other_row['STY'] >= row['STY'] and
                    other_row['E_factor'] <= row['E_factor']) and \
               (other_row['STY'] > row['STY'] or
                other_row['E_factor'] < row['E_factor']
               ):
                is_pareto[i] = False
                break

    return df[is_pareto].sort_values(by='E_factor')

def plot_pareto_front(
        df: pd.DataFrame,
        all_points_df: pd.DataFrame,
        iteration: int,
        num_lhs_points: int
) -> str:
    """Plot the Pareto front from a DataFrame with 'STY' and 'E_factor'
    columns.
    """
    pareto_df = pareto_front(df)
    initial_lhs_points = df.iloc[:num_lhs_points] if len(
        df
    ) >= num_lhs_points else df

    plt.figure(figsize=(8, 6))
    plt.scatter(
        initial_lhs_points['STY'], initial_lhs_points['E_factor'],
        marker='s', color='black', label=f'Initial LHS Points',
        s=64, zorder=4
    )
    plt.scatter(
        df['STY'], df['E_factor'],
        marker='x', color='blue', label=f'TSEMO Points',
        alpha=0.7, zorder=3
    )
    plt.scatter(
        pareto_df['STY'], pareto_df['E_factor'],
        facecolors='orange', edgecolors='red', s=100, linewidths=1.5,
        label='Pareto Front', zorder=5, marker='o'
    )
    plt.plot(
        pareto_df['STY'], pareto_df['E_factor'],
        color='red', linewidth=2, linestyle='--', zorder=2,
        label='Pareto Curve'
    )
    plt.xlim(0, 14500)
    plt.ylim(0, 3)
    plt.xlabel('STY [kg/m³/h]', fontsize=18)
    plt.ylabel('E-Factor', fontsize=18)
    plt.title(f'In-silicio closed loop multi-objective optimization - Iteration {iteration}', fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.legend(fontsize=15)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    file_path = (
        f"plots/system_1/pareto_fronts/"
        f"pareto_front_iteration_{iteration}.png"
    )
    plt.savefig(
        file_path,
        dpi=300
    )
    plt.close()
    return file_path

if __name__ == "__main__":

    results_dataframe = pd.read_csv('data/system_1/tsemo_results.csv')
    num_starting_points = 5
    num_lhs_points = 20

    frame_paths = []

    for i in range(num_starting_points, len(results_dataframe) + 1):
        current_df_subset = results_dataframe.iloc[:i]
        frame_path = plot_pareto_front(
            current_df_subset, results_dataframe, i, num_lhs_points
        )
        frame_paths.append(frame_path)

    with imageio.get_writer(
            "plots/system_1/pareto_front_animation.gif",
            mode='I', duration=0.4, loop=0,
    ) as writer:
        for frame_path in frame_paths:
            image = imageio.imread(frame_path)
            writer.append_data(image)



