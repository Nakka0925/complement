import csv
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from pandas.errors import EmptyDataError


# セルの色付け関数（マイナス値を赤色にする）
def color_negative(val):
    if val == "null":
        f'color: black'
        return
    
    val = float(val)  # 文字列を数値に変換
    color = 'red' if val < 0 else 'black'
    
    return f'color: {color}'


def gen_table(df, target, threshold):
    G = nx.Graph()

    # xとyを整数として抽出し、新しい変数に保存
    pairs = df['pair'].str.extract(r'\((\d+), (\d+)\)').astype(int)
    x = pairs[0]
    y = pairs[1]
    n = max(y)


    # 初期値 "inf" のデータテーブルを作成
    data_table = np.full((n+1, n+1), "null", dtype=object)

    #　全てのノードを追加 
    for i in range(n): G.add_node(i)

    for i, j, tau in zip(x, y, df['tau']):
        data_table[i][j] = f"{tau:.2f}" 
        if tau >= threshold:
            G.add_edge(i, j, weight=tau) # グラフにノードを追加

    indices = sorted(list(set(y)) + [0])
    out_df = pd.DataFrame(data_table, index=indices, columns=indices)

    # スタイルを適用して表示
    styled_df = out_df.style.applymap(color_negative)
    styled_df.to_html(f'table/styled_output_{target}.html')


    # グラフの描画
    plt.figure(figsize=(15, 8))
    pos = nx.spring_layout(G, weight='weight', seed=0, k=0.3)  # ノードの配置を決定
    # nx.draw(G, pos, with_labels=True, node_size=200, font_size=5, font_weight='bold', alpha=0.6)
    nx.draw_networkx_nodes(G, pos,alpha=0.6, node_size=300)
    nx.draw_networkx_labels(G, pos, font_size=10, font_weight="bold")
    nx.draw_networkx_edges(G, pos, alpha=0.4)
    # nx.draw_networkx_edge_labels(G, pos, font_size=7, edge_labels={(i, j): f"{G[i][j]['weight']:.2f}" for i, j in G.edges()})
    # plt.title(f"Graph of Nodes with Strong Relationships{threshold}")
    plt.axis("off")
    plt.savefig(f'grapth/{target}_{threshold}.png')
    plt.close()


if __name__ == "__main__":
    dist = "/home/nakanishi/M2/generate-img/data"
    df = pd.read_csv(f"{dist}/csv/class_sum.csv", encoding="shift-jis")
    
    threshold = 0.5

    for target in df["class"]:
        try:
            data = pd.read_csv(f'../result/score/{target}_kendall_score.csv')
            gen_table(data, target, threshold)
        except EmptyDataError as e:
            print(f"{target}はデータが1つしかないため計算できませんでした")
            print(f"エラー原因: {e}")
