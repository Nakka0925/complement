import pandas as pd
import sys,re
import itertools
from Bio import SeqIO
import pandas as pd
import shutil


def complement_cacl(data, cls):
    dist = "/home/nakanishi/M2/generate-img/data"
    file_path = f'classes/array_to_acc/{cls}_rna.csv'
    df = pd.read_csv(file_path)

    for target in range(len(df)):
        accs = eval(df['accessions'][target])
        rna_list = df["product_array"][target].split(", ")
        strand_cnt = [0 for _ in range(len(rna_list))]
        for acc in accs:
            current = 0
            # GenBankファイルのパスを指定 
            # NC_013934_0 -> NC_013934
            file_path = f"{dist}/gbk/{acc[:-2]}.gbk"
            # GenBankファイルを読み込む
            records = SeqIO.parse(file_path, "genbank")

            for record in records:
                for feature in record.features:
                    if (feature.type == "tRNA"):
                        product_type = feature.qualifiers["product"][0]
                        if product_type == "tRNA-Asx":
                            print("Asx")
                            continue
                        elif product_type == "tRNA-OTHER":
                            print("OTHER")
                            continue
                        elif product_type == "tRNA-Glx":
                            print("Glx")
                            continue
                        elif product_type == "tRNA-Sec":
                            print("Sec")
                            continue
                        elif product_type == "tRNA-Xle":
                            print("Xle")
                            continue
                        elif product_type == "tRNA-iMet":
                            print("iMet")
                            continue
                        strand = feature.location.strand
                        if strand == 1:
                            strand_cnt[current] = strand_cnt[current] + 1
                        if strand == -1:
                            strand_cnt[current] = strand_cnt[current] - 1
                        current = current + 1

        for i in range(len(strand_cnt)):
            if strand_cnt[i] < 0:
                rna_list[i] = "complement-" + rna_list[i]

        # 更新されたリストをDataFrameに保存
        df.at[target, "product_array"] = ", ".join(rna_list)
    
    df.to_csv(f"classes_comp/array_to_acc/{cls}_rna.csv", index=False)


if __name__ == '__main__':
    dist = "/home/nakanishi/M2/generate-img/data"
    df = pd.read_csv(f"{dist}/csv/class_sum.csv", encoding="shift-jis")

    for cls in df["class"]:
        # file_path = f'result/score/{cls}_lcs_average.csv'
        file_path = f'classes_comp/array_cnt/{cls}_rna.csv'
        data = pd.read_csv(file_path)        
        
        complement_cacl(data, cls)