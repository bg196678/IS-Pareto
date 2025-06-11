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
        iteration: int
) -> str:
    """Plot the Pareto front from a DataFrame with 'STY' and 'E_factor'
    columns.
    """
    pareto_df = pareto_front(df)

    plt.figure(figsize=(10, 6))
    plt.scatter(
        all_points_df['STY'], all_points_df['E_factor'],
        c='gray', label='All Data Points', alpha=0.3
    )
    plt.scatter(
        df['STY'], df['E_factor'],
        c='blue', label=f'Points Considered ({len(df)})', alpha=0.6
    )
    plt.scatter(
        pareto_df['STY'], pareto_df['E_factor'],
        c='red', label='Pareto Front', s=100, edgecolors='black', zorder=5
    )
    plt.plot(
        pareto_df['STY'], pareto_df['E_factor'],
        c='red', linestyle='-', marker='', zorder=4
    )
    plt.xlim(0, 14500)
    plt.ylim(0, 3)
    plt.xlabel('STY')
    plt.ylabel('E-Factor')
    plt.title('Pareto Front: E-Factor vs. STY')
    plt.legend()
    plt.grid(True)
    file_path = f"plots/pareto_fronts/pareto_front_iteration_{iteration}.png"
    plt.savefig(
        file_path,
        dpi=300
    )
    plt.close()
    return file_path

if __name__ == "__main__":

    results_dataframe = pd.read_csv('data/tsemo_results.csv')
    num_starting_points = 5

    frame_paths = []

    for i in range(num_starting_points, len(results_dataframe) + 1):
        current_df_subset = results_dataframe.iloc[:i]
        frame_path = plot_pareto_front(current_df_subset, results_dataframe, i,)
        frame_paths.append(frame_path)

    with imageio.get_writer(
            "plots/pareto_front_animation.gif", mode='I', duration=0.2, loop=0
    ) as writer:
        for frame_path in frame_paths:
            image = imageio.imread(frame_path)
            writer.append_data(image)



