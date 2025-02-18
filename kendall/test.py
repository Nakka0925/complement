import csv, sys
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from pandas.errors import EmptyDataError
from collections import deque

complement_sum = 0
incorrect_value = 0
non_data = 0
ha = 0

def find_complementary(grupe, df, df2, threshold):
    max_grupe = max(grupe, key=lambda x: x['cnt_sum'])
    max_node = 0
    tmp = 0
    result = []
    global complement_sum
    global incorrect_value
    global non_data
    global ha
    # max_node を見つけるループ
    # for node in max_grupe['grupe_list']:
    #     if tmp < df2['count'][node]:
    #         max_node = node
    #         tmp = df2['count'][node]
    
    
    # max_node と各ノードの tau 値を計算
    for node in range(len(df2)):
        if node == max_node:
            continue
        
        # if node not in max_grupe['grupe_list']:
        # max_node と node のペアを検索して tau 値を取得する
        pair = f"({min(max_node, node)}, {max(max_node, node)})"
        
        p_value = df.loc[(df['pair'] == pair), 'p_value'].values[0]
        tau_value = df.loc[(df['pair'] == pair), 'tau'].values[0]
        
        if pd.isna(p_value): # データがない場合(共通遺伝子がないor共通遺伝子が1個)
            non_data = non_data + df2['count'][node]
            continue
        
        if 0.05 <= p_value or pd.isna(p_value): # p値から有意差なしの場合
            ha = ha + 1
            incorrect_value = incorrect_value + df2['count'][node]
            continue
        
        # df から該当するペアの tau 値を取得
        count = df2['count'][node]
        # print(pair,tau_value)
        if tau_value <= -threshold:
            result.append({
                "max_node-node": max_node,
                "node": node, 
                "tau": tau_value,
                "count": count
            })
            complement_sum = complement_sum + count
        # if -0.2 < tau_value :
        #     result.append({
        #         "max_node-node": max_node,
        #         "node": node, 
        #         "tau": tau_value,
        #         "count": count
        #     })
            
    result_df = pd.DataFrame(result)
    result_df.to_csv(f'../complement/{target}_complement.csv', index=False)
    

def bfs(df, df2, target, threshold):
    # xとyを整数として抽出し、新しい変数に保存
    pairs = df['pair'].str.extract(r'\((\d+), (\d+)\)').astype(int)
    x = pairs[0]
    y = pairs[1]
    n = max(y) + 1 # 全ノード数
    m = len(pairs) # エッジの数
    dist = [-1] * n # -1はまだ探索されていないことを表す
    grupe = []
    grupe_index = 0

    
    g = [[] for _ in range(n)] # 隣接リスト: threshold以上の場合辺を結ぶ
    for i, j, tau in zip(x,y,df['tau']):
        if threshold <= tau:
            g[i].append(j)
            g[j].append(i)

    
    for k in range(n):
        grupe_list = [] 
        cnt = 0
        grupe_cnt_sum = 0
        if (dist[k] != -1):
            continue
        dist[k] = 0
        queue = deque([k])
        while queue:
            cnt = cnt + 1
            v = queue.popleft()
            grupe_list.append(v)
            grupe_cnt_sum = grupe_cnt_sum + df2['count'][v]

            for i in g[v]:
                if dist[i] == -1:
                    dist[i] = dist[v] + 1
                    queue.append(i)
        grupe.append({
            'index': grupe_index,
            'num': cnt,
            'cnt_sum': grupe_cnt_sum,
            'grupe_list': grupe_list
        })
        # print(len(grupe_list))
        grupe_index = grupe_index + 1
    find_complementary(grupe, df, df2, threshold)
    

if __name__ == "__main__":
    dist = "/home/nakanishi/M2/generate-img/data"
    df = pd.read_csv(f"{dist}/csv/class_sum.csv", encoding="shift-jis")
    
    threshold = 0.2

    for target in df["class"]:
        try:
            data = pd.read_csv(f'../result/score/{target}_kendall_score.csv')
            data2 = pd.read_csv(f'../classes/array_cnt/{target}_rna.csv')
            bfs(data, data2, target, threshold)
        except EmptyDataError as e:
            print(f"{target}はデータが1つしかないため計算できませんでした")
            print(f"エラー原因: {e}")
    print("相補鎖:",complement_sum)
    print("0.5<=p_value", incorrect_value)
    print("non_data:", non_data)
    print("相補鎖の確率:", complement_sum / (15887 - incorrect_value - non_data) * 100, "%")