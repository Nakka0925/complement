import pandas as pd
import itertools
from scipy.stats import kendalltau
import matplotlib.pyplot as plt
import sys


under_5 = []
under_10 = []
under_15 = []
under_20 = []


def convert_to_numeric_order(species_order, reference_list):
    return [reference_list[trna] for trna in species_order if trna in reference_list]


def preproccesing(data):
    
    unique_gene = {}
    preproccesing_trna_list = []
    for row in data:
        genes = row.split(', ')
        gene_counts = {}
        preproccesing_trna = []
        for idx, gene in enumerate(genes):
            if gene not in gene_counts:
                gene_counts[gene] = 1
            else:
                gene_counts[gene] += 1
            unique_gene_name = f"{gene}{gene_counts[gene]}"
            preproccesing_trna.append(unique_gene_name)
            
            # if 1 < gene_counts[gene]: # 同じ遺伝子が複数の場合はその遺伝子を除外
            #     if f"{gene}1" in unique_gene:
            #         del unique_gene[f"{gene}1"]
            #     continue
            
            # unique_gene に全体の遺伝子の出現回数を記録
            if unique_gene_name not in unique_gene:
                unique_gene[unique_gene_name] = 1
            else:
                unique_gene[unique_gene_name] += 1
        preproccesing_trna_list.append(preproccesing_trna)
    filtered_genes = [gene for gene, count in unique_gene.items() if count == len(data)]
    
    genes_to_idx = {}
    
    for idx, gene in enumerate(filtered_genes):
        genes_to_idx[gene] = idx
        
    result_trna_0 =  convert_to_numeric_order(preproccesing_trna_list[0],genes_to_idx)
    result_trna_1 =  convert_to_numeric_order(preproccesing_trna_list[1],genes_to_idx)

    return result_trna_0,result_trna_1


def kendall_test(data, cls):
    product_array_list = data["product_array"]
    result = []
    each_seq_tau = {}
    
    for i in range(len(product_array_list)):
        for j in range(i+1, len(product_array_list)):
            product_array_i, product_array_j = preproccesing([product_array_list[i], product_array_list[j]])
            
            tau, p_value = kendalltau(product_array_i, product_array_j, method='exact')
            result.append({
                'pair': (i,j),
                'tau': tau,
                'p_value': p_value
            })
            if i == 0:
                product_len = len(product_array_i)
                if product_len <= 5:
                    under_5.append(p_value)
                elif product_len <= 10:
                    under_10.append(p_value)
                elif product_len <= 15:
                    under_15.append(p_value)
                elif 15 < product_len:
                    under_20.append(p_value)
            
            each_seq_tau.setdefault(i,[])
            each_seq_tau.setdefault(j,[])

            each_seq_tau[i].append(tau)
            each_seq_tau[j].append(tau)
    
    tau_ave = []
    for seq_id in each_seq_tau:
        ave = sum(each_seq_tau[seq_id]) / len(each_seq_tau[seq_id])
        # if ave < -0.2:
        tau_ave.append({'id': seq_id, 
                        'ave_tau': ave
        })
        # if ave < -0.3:
            
        
    tau_ave_df = pd.DataFrame(tau_ave)
    tau_ave_df.to_csv(f'../ave_kendall_test/{cls}_kendall_score.csv', index=False)
    result_df = pd.DataFrame(result)
    result_df.to_csv(f'../result/score/{cls}_kendall_score.csv', index=False)
    
            
if __name__ == '__main__':
    dist = "/home/nakanishi/M2/generate-img/data"
    df = pd.read_csv(f"{dist}/csv/class_sum.csv", encoding="shift-jis")

    for cls in df["class"]:
        file_path = f'../classes/array_cnt/{cls}_rna.csv'
        data = pd.read_csv(file_path)  
        kendall_test(data, cls)
    
    print("5以下", sum(under_5) / len(under_5))
    print("10以下", sum(under_10) / len(under_10))
    print("15以下", sum(under_15) / len(under_15))
    print("15以上", sum(under_20) / len(under_20))
    under_20.sort()
    print(under_20)
    print(len(under_15))
