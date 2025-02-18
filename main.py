import pandas as pd
import sys,re
import itertools
from Bio import SeqIO
import pandas as pd
import shutil


sequence_analysis_each_class = {}
array_type_each_species = {}
pattern = "voucher"

dist = "/home/nakanishi/M2/generate-img/data"
df = pd.read_csv(f"{dist}/dataset0/data_all.csv", encoding="shift-jis")
accs = df["accession"]

acc_to_class = {row["accession"][:-2] : row["class"] for _, row in df.iterrows()}


for acc in accs:
    # GenBankファイルのパスを指定 
    # NC_013934_0 -> NC_013934
    file_path = f"{dist}/gbk/{acc[:-2]}.gbk"

    # GenBankファイルを読み込む
    records = SeqIO.parse(file_path, "genbank")
    # 相補鎖のポジションを格納するリスト
    sequence_analysis_each_class.setdefault(acc_to_class[acc[:-2]], {}) # keyがない場合[]で要素を初期化

    array_type_each_species.setdefault(acc_to_class[acc[:-2]], {})

    for record in records:
        product_array = []
        # match_result = re.search(pattern,record.description)
        sequence_length = len(record.seq)
        for feature in record.features:
            if (feature.type == "tRNA"):
                product_type = feature.qualifiers["product"][0]
                if product_type == "tRNA-Asx":
                    # print("Asx")
                    continue
                elif product_type == "tRNA-OTHER":
                    # print("OTHER")
                    continue
                elif product_type == "tRNA-Glx":
                    # print("Glx")
                    continue
                elif product_type == "tRNA-Sec":
                    # print("Sec")
                    continue
                elif product_type == "tRNA-Xle":
                    # print("Xle")
                    continue
                elif product_type == "tRNA-iMet":
                    # print("iMet")
                    continue
                product_array.append(product_type)
        
        
        if len(product_array) == 0:
            print(acc_to_class[acc[:-2]])
            print(acc)
            continue

        sequence_analysis_each_class[acc_to_class[acc[:-2]]].setdefault(tuple(product_array), 0) #keyに動的な構造体はつかえない
        sequence_analysis_each_class[acc_to_class[acc[:-2]]][tuple(product_array)] += 1

        array_type_each_species[acc_to_class[acc[:-2]]].setdefault(tuple(product_array), [])
        array_type_each_species[acc_to_class[acc[:-2]]][tuple(product_array)].append(acc)

# 結果の表示
for key, value in sequence_analysis_each_class.items():
    data = {"product_array": list(value.keys()), "count": list(value.values())}
    positions_df = pd.DataFrame(data)
    
    # タプルを文字列に変換して保存する
    positions_df["product_array"] = positions_df["product_array"].apply(lambda x: ", ".join(x))
    positions_df = positions_df.sort_values(by="count", ascending=False).reset_index(drop=True)  # インデックスをリセット
    positions_df.to_csv(f"classes/array_cnt/{key}_rna.csv")

for key, value in array_type_each_species.items():
    data = {"product_array": list(value.keys()), "accessions": list(value.values())}
    positions_df = pd.DataFrame(data)
    
    # タプルを文字列に変換して保存する
    positions_df["product_array"] = positions_df["product_array"].apply(lambda x: ", ".join(x))
    positions_df = positions_df.sort_values(by="accessions", key=lambda x: x.str.len(), ascending=False).reset_index(drop=True)
    positions_df.to_csv(f"classes/array_to_acc/{key}_rna.csv")
