from matplotlib import hatch
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
import numpy as np
import sys

def create_each_gene_color(df): #TODO 出現する遺伝子は20種類であることが前処理で判明している
    all_unique_gene = []
    for cls in df["class"]:
        file_path = f'classes_comp/array_cnt/{cls}_rna.csv'
        data = pd.read_csv(file_path)
        
        # unique_gene = list(set(trna for sequence in data['product_array'] for trna in sequence.split(', ')))
        # unique_gene = list(set([rna[11:] if rna[:11] == "complement-" else rna for rna in unique_gene]))
        # 重複する遺伝子を削除、complement-を外す
        unique_gene = list(
            {rna[11:] if rna[:11] == "complement-" else rna for sequence in data['product_array'] for rna in sequence.split(', ')}
        )

        all_unique_gene = list(set(all_unique_gene + unique_gene))
    
    # Use a predefined colormap (e.g., 'tab20') to generate distinct colors
    cmap = plt.get_cmap("tab20")
    all_unique_gene.sort()
    colors = {trna: cmap(i / len(all_unique_gene)) for i, trna in enumerate(all_unique_gene)}
        
    return colors
    
    
# Define a function to create a colored sequence plot
def plot_colored_sequences(data, colors, cls):
    plt.figure(figsize=(15, max(len(data) * 0.4, 8.0)))
    
    for idx, row in data.iterrows():
        sequence = row['product_array'].split(', ')
        for position, gene in enumerate(sequence):
            hatch_pattern = None
            # geneが相補鎖の情報の場合
            if gene[:11] == "complement-": 
                gene = gene[11:]
                hatch_pattern = "xxx"
            # Plot a colored rectangle for each tRNA
            plt.barh(idx, 1, left=position, 
                     color=colors[gene], height=0.8, 
                     edgecolor="black", linewidth=1,
                     hatch=hatch_pattern
                     )
       
    
    # Set y-axis to display sequence index or other identifiers
    plt.yticks(range(len(data)), data.index, fontsize=20)
    plt.xticks(fontsize=20)
    plt.xlabel('tRNA Position', fontsize=25)
    plt.ylabel('Sequence Index', fontsize=25)
    # plt.title('Colored tRNA Sequences')
    plt.tight_layout()
    
    
    # y軸を反転
    plt.gca().invert_yaxis()

    legend_patches = [mpatches.Patch(color=color, label=trna) for trna, color in sorted(colors.items())]
    # 相補鎖の説明用のハッチ模様を追加
    hatch_patch = mpatches.Patch(facecolor='white', edgecolor='black', hatch='xxx', label='gene on the other strand')
    legend_patches.append(hatch_patch)
    
    plt.legend(handles=legend_patches, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, fontsize=18)
    
    plt.savefig(f"result-test/graph/{cls}.png",bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    dist = "/home/nakanishi/M2/generate-img/data"
    df = pd.read_csv(f"{dist}/csv/class_sum.csv", encoding="shift-jis")
    
    unique_gene_color = create_each_gene_color(df)

    for cls in df["class"]:
        file_path = f'classes_comp/array_to_acc/{cls}_rna.csv'
        data = pd.read_csv(file_path)
        
        plot_colored_sequences(data, unique_gene_color, cls)