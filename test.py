from Bio import SeqIO


# GenBankファイルのパスを指定 
file_path = "NC_012346.gbk"
# GenBankファイルを読み込む
records = SeqIO.parse(file_path, "genbank")

trna_list = []

for record in records:
    for feature in record.features:
        if (feature.type == "tRNA"):
            start = int(feature.location.start)
            end = int(feature.location.end)
            product_type = feature.qualifiers["product"][0]
            trna_list.append((product_type,end-start))

print(trna_list)