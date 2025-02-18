import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import re, os, sys
from pathlib import Path


def create_each_gene_color(df):
    all_unique_gene = []
    for cls in df["class"]:
        file_path = f'../classes_comp/array_cnt/{cls}_rna.csv'
        data = pd.read_csv(file_path)
        
        unique_gene = list(
            {rna[11:] if rna[:11] == "complement-" else rna for sequence in data['product_array'] for rna in sequence.split(', ')}
        )

        all_unique_gene = list(set(all_unique_gene + unique_gene))
    
    # Use a predefined colormap (e.g., 'tab20') to generate distinct colors
    cmap = plt.get_cmap("tab20")
    all_unique_gene.sort()
    colors = {trna: cmap(i / len(all_unique_gene)) for i, trna in enumerate(all_unique_gene)}
        
    return colors


def plot_colored_sequences(data, colors, cls, target):
    print(target)
    plt.figure(figsize=(15, 8))
    rna_coo = {}

    for idx, i in enumerate(target):
        sequence = data['product_array'][i].split(', ')
        for position, gene in enumerate(sequence):
            hatch_pattern = None
            # geneが相補鎖の情報の場合
            if gene[:11] == "complement-": 
                gene = gene[11:]
                hatch_pattern = "xxx"
            plt.barh(idx, 1, left=position, 
                     color=colors[gene], height=0.1, 
                     edgecolor="black", linewidth=1,
                     hatch=hatch_pattern
                    )
            rna_coo.setdefault(gene,[])

            rna_coo[gene].append({
                'y':idx,
                'x': position+0.5,
                'gene_type': hatch_pattern
            })

    
    # 各tRNAのペアを結ぶ
    for rna in rna_coo:
        y_1_rna_coo = [i for i in rna_coo[rna] if i['y'] == 1]
        y_0_rna_coo = [i for i in rna_coo[rna] if i['y'] == 0]
        for i in y_0_rna_coo:
            for j in y_1_rna_coo:
                x = [i['x'], j['x']]
                y = [i['y'], j['y']]
                # plt.plot(x,y,color=colors[rna])
                if i["gene_type"] == j['gene_type']:
                    plt.plot(x,y,color='blue')
                else:
                    plt.plot(x,y,color='red')
                plt.yticks(y,[f"arr_index_{target[0]}",f"arr_index_{target[1]}"])
    plt.yticks(fontsize=20)
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
    
    # plt.savefig(f"{cls}.png",bbox_inches='tight')
    new_dst = Path('compare_graph') / cls
    new_dst.mkdir(parents=True, exist_ok=True)

    plt.savefig(f"compare_graph/{cls}/{cls}_{target[0]}_{target[1]}.png", bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    dist = "/home/nakanishi/M2/generate-img/data"
    df = pd.read_csv(f"{dist}/csv/class_sum.csv", encoding="shift-jis")

    unique_gene_color = create_each_gene_color(df)

    new_dst = Path('compare_graph')
    new_dst.mkdir(parents=True, exist_ok=True)

    # 単数比較
    # cls = "Insecta"
    # target = [0,16]

    # file_path = f'../classes/array_cnt/{cls}_rna.csv'
    # data = pd.read_csv(file_path)
    
    # plot_colored_sequences(data, unique_gene_color, cls, target)

    # 複数比較
    flag = False
    for cls in df["class"]:
        # if cls == "Anthozoa":
        #     flag = True
        cls = "Aves"
        print(cls)
        # if flag:
        file_path = f'../classes_comp/array_to_acc/{cls}_rna.csv'
        data = pd.read_csv(file_path)
        for target in range(1,len(data)):
            plot_colored_sequences(data, unique_gene_color, cls, [0,target])
        sys.exit()
