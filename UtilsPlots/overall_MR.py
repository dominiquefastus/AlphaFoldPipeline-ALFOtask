import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from adjustText import adjust_text
import argparse

def read_data(file_path):
    """
    Reads the CSV file and returns the data.
    """
    return pd.read_csv(file_path)

def prepare_data(data):
    """
    Transposes and cleans the data for analysis.
    """
    data_transposed = data.transpose()
    data_transposed.columns = data_transposed.iloc[0]
    data_transposed = data_transposed.drop(data_transposed.index[0])
    return data_transposed.dropna(axis=1, how='all')

def create_scatter_plot(data):
    """
    Create and save the scatter plot.
    """
    scatter_data = data[['Type', 'LLG', 'TFZ']]
    scatter_data = scatter_data.dropna()
    scatter_data['LLG'] = pd.to_numeric(scatter_data['LLG'], errors='coerce')
    scatter_data['TFZ'] = pd.to_numeric(scatter_data['TFZ'], errors='coerce')

    plt.figure(figsize=(11, 8), dpi=300)
    ax = plt.subplot(111)
    types = scatter_data['Type'].unique()
    for t in types:
        subset = scatter_data[scatter_data['Type'] == t]
        plt.scatter(subset['TFZ'], subset['LLG'], s=150, label=t, alpha=0.6)
        for i in subset.index:
            plt.text(subset['TFZ'][i], subset['LLG'][i], i, fontsize=9, ha='right', va='bottom')

    plt.axhline(y=40, color='r', linestyle='-', label='LLG 40 (minimum for correct solution)')
    plt.axhline(y=60, color='orange', linestyle='-', label='LLG 60 (difficult problems)')
    plt.axhline(y=120, color='green', linestyle='-', label='LLG 120 (ideal minimum)')
    
    plt.axvline(x=5, color='r', linestyle='-', label='TFZ 5 (minimum for correct/solved solution)')
    plt.axvline(x=6, color='orange', linestyle='-', label='TFZ 6 possibly correct/solved')
    plt.axvline(x=8, color='green', linestyle='-', label='TFZ 8 (deffinately correct/solved)')
    
    # plt.title('Molecular Replacement Quality: TFZ Score vs LLG (with Protein Labels)')
    plt.xlabel('TFZ Score')
    plt.ylabel('LLG Score')
    
    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,
                    box.width, box.height * 0.9])

    # Put a legend below current axis
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),
            fancybox=True, shadow=True, ncol=3)

    plt.grid(True)
    plt.savefig('tfz_vs_llg.png')

def create_gradient_plot(data):
    """
    Create and save the gradient plot.
    """
    def draw_custom_gradient_bar_with_corrected_scales(ax, metric_values, y, cmap, label, worst_val, best_val):
        gradient = np.linspace(0, 1, 256)
        gradient = np.vstack((gradient, gradient))
        ax.imshow(gradient, aspect='auto', cmap=cmap, extent=(0.05, 0.95, y - 0.15, y + 0.15))

        sorted_values = sorted(metric_values.items(), key=lambda x: x[1])
        label_positions = np.linspace(y + 0.32, y + 0.64, len(metric_values))

        for (protein, value), label_y in zip(sorted_values, label_positions):
            value = 0.95 * value
            ax.plot(value, y, 'v', color='black', markersize=5)
            ax.text(value, label_y - 0.04, f'{protein}', va='center', ha='center', fontsize=5, color='black')
            ax.plot([value, value], [y, label_y - 0.15], color='black', linestyle='-', linewidth=0.5)

        ax.text(0.06, y - 0.25, f'Worst ({worst_val})', va='center', ha='center', fontsize=9, color='black')
        ax.text(0.94, y - 0.25, f'Best ({best_val})', va='center', ha='center', fontsize=9, color='black')
        ax.text(0.5, y + 0.55, label, va='center', ha='center', fontsize=8, color='black', fontweight='bold')

    custom_cmap_new = LinearSegmentedColormap.from_list("custom_cmap_new", ["tomato", "mediumseagreen"])

    metrics = {
        'Clashscore': data['Clashscore'].dropna().astype(float),
        'Ramachandran outliers': data['Ramachandran outliers'].dropna().str.rstrip('%').astype(float),
        'R-free': data['R-free'].dropna().astype(float),
        'Rotamer outliers': data['Rotamer outliers'].dropna().str.rstrip('%').astype(float)
    }

    normalized_metrics = {}
    normalized_metrics['R-free'] = ((metrics['R-free'] - 0.6) / (0 - 0.6)).to_dict()
    normalized_metrics['Clashscore'] = ((metrics['Clashscore'] - 20) / (0 - 20)).to_dict()
    normalized_metrics['Ramachandran outliers'] = ((metrics['Ramachandran outliers'] - 1) / (0 - 1)).to_dict()
    normalized_metrics['Rotamer outliers'] = ((metrics['Rotamer outliers'] - 5) / (0 - 5)).to_dict()

    fig, ax = plt.subplots(figsize=(9, 5), dpi=300)
    draw_custom_gradient_bar_with_corrected_scales(ax, normalized_metrics['R-free'], 4, custom_cmap_new, 'R-free', "0.6", "0")
    draw_custom_gradient_bar_with_corrected_scales(ax, normalized_metrics['Clashscore'], 1, custom_cmap_new, 'Clashscore', "20", "0")
    draw_custom_gradient_bar_with_corrected_scales(ax, normalized_metrics['Ramachandran outliers'], 3, custom_cmap_new, 'Ramachandran Outliers', "1%", "0%")
    draw_custom_gradient_bar_with_corrected_scales(ax, normalized_metrics['Rotamer outliers'], 2, custom_cmap_new, 'Rotamer Outliers', "5%", "0%")

    # ax.set_title('Protein Quality Metrics with Custom Gradient Scale')
    ax.set_xticks([0.05, 0.95])
    ax.set_xticklabels(['Worst', 'Best'])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 5)
    ax.set_yticks([])
    ax.set_xlabel('Metric Scale')
    ax.set_ylabel('Metrics')

    for spine in ax.spines.values():
        spine.set_visible(True)
    ax.tick_params(axis=u'both', which=u'both', length=0)
    plt.savefig('custom_gradient_plot.png')


def main():
    parser = argparse.ArgumentParser(description="Generate protein metrics plots from CSV data.")
    parser.add_argument("csv_file", help="Path to the CSV file containing protein metrics.")
    args = parser.parse_args()

    data = read_data(args.csv_file)
    prepared_data = prepare_data(data)

    create_scatter_plot(prepared_data)
    create_gradient_plot(prepared_data)
    print("Plots have been generated and saved.")

if __name__ == "__main__":
    main()
