import pandas as pd
import sys,re
import itertools
from Bio import SeqIO
import pandas as pd
import shutil

dist = "/home/nakanishi/M2/generate-img/data"
cls = "Aves"
targets = [0,9]
file_path = f'../classes/array_to_acc/{cls}_rna.csv'
df = pd.read_csv(file_path)

for target in targets:
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
                    strand = feature.location.strand
                    if strand == 1:
                        strand_cnt[current] = strand_cnt[current] + 1
                    if strand == -1:
                        strand_cnt[current] = strand_cnt[current] - 1
                    current = current + 1

    for i in range(len(strand_cnt)):
        if 0 <= strand_cnt[i]:
            rna_list[i] = "strand-" + rna_list[i]

    print(rna_list)