#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from scipy.stats import ttest_ind

from sklearn.decomposition import PCA
# %%
geneset = pd.read_csv('GeneSetSignature.gmt', sep="\t", header=None, index_col=0)
geneset = geneset.iloc[:, 1:]
#%%
edema = "NGPT4,TJP2,OCLN,BNIP3,HMGB3,PIK3R3,AKT3,RICTOR,FOXO4,FOXO1,GSK3B,KRAS,MAPK8,MAPK9,IKBKB,CXCL10,CXCL11,ELK3,SDSL,TSPO,CASP8,CDKN1A,ARPC3,PFN1,TMSB10,TMSB4X,S100A11,S100A4,COTL1,APOA1,APOC1,APOBEC3G,TNFSF12,TNFRSF14,CASP4,CSTA,CTBS,CTSD,CD48"
enhancing ="SLC2A3,NDRG1,SNAI2,ACTN4,LIMK1,TESK1,DUSP3,PRKACA,MARK2,E2F8,NR4A3,TEAD3,TCF7,EIF3E"
nonenhancing = "FOXO1,MDM2,HIF3A,NGDN,CALM2,NCK2,MYOZ2,CGRRF1,TIMM9,EMX2,CLDN5,JAK3,FOXP3,DUSP2,DUSP18,PRICKLE1,WWTR1 (TAZ),SMAD7,RUNX3,NLRP3,ANGPT1,MMP9,VCL,EPS8,ARHGAP24,ALDH1L2,GADL1,VDR,TNS3,ITGA6,FN1,COL12A1,PTGER4,TGFBI,GLIPR1,GSN,MMP10"
#%%
#calculate percentage overlap of genes in edema with rows in geneset and which row overlaps with which element in edema
def calculate_overlap(geneset, compartment_genes):
    overlap_results = {}
    compartment_genes_set = set(compartment_genes.split(','))
    for index, row in geneset.iterrows():
        row_genes = set(row.dropna().values)
        overlap = compartment_genes_set.intersection(row_genes)
        overlap_percentage = len(overlap) / len(compartment_genes_set) * 100
        overlap_results[index] = {
            'overlap_percentage': overlap_percentage,
            'overlapping_genes': overlap
        }
    return overlap_results
#%%
# Calculate overlap for edema genes
overlap_results_edema = calculate_overlap(geneset, edema)
# Convert results to DataFrame for better visualization
overlap_df_edema = pd.DataFrame.from_dict(overlap_results, orient='index')
# Sort by overlap percentage
overlap_df_edema = overlap_df.sort_values(by='overlap_percentage', ascending=False)

overlap_results_enhancing = calculate_overlap(geneset, enhancing)
overlap_df_enhancing = pd.DataFrame.from_dict(overlap_results_enhancing, orient='index')
overlap_df_enhancing = overlap_df_enhancing.sort_values(by='overlap_percentage', ascending=False)

overlap_results_nonenhancing = calculate_overlap(geneset, nonenhancing)
overlap_df_nonenhancing = pd.DataFrame.from_dict(overlap_results_nonenhancing, orient='index')
overlap_df_nonenhancing = overlap_df_nonenhancing.sort_values(by='overlap_percentage', ascending=False)

#write all three dataframes to tsv files
overlap_df_edema.to_csv('overlap_edema.tsv', sep='\t')
overlap_df_enhancing.to_csv('overlap_enhancing.tsv', sep='\t')
overlap_df_nonenhancing.to_csv('overlap_nonenhancing.tsv', sep='\t')
# %%
df1_scores = pd.read_csv('AUCell_scores_GBM_pathways.txt', sep="\t", index_col=0)
df2_scores = pd.read_csv('AUCell_scores_LGG_pathways.txt', sep="\t", index_col=0)
df3_scores = pd.read_csv('AUCell_scores_BrainMet_pathways.txt', sep="\t", index_col=0)
#%%
concatenated_scores = pd.concat([df1_scores, df2_scores, df3_scores], axis=1)
concatenated_scores = concatenated_scores.T
# perform z-score normalization for each pathway across all samples
concatenated_scores = concatenated_scores.apply(lambda x: (x - x.mean()) / x.std(), axis=0)
#add labels to the concatenated scores 
labels = ['GBM'] * df1_scores.T.shape[0] + ['LGG'] * df2_scores.T.shape[0] + ['BrainMetastases'] * df3_scores.T.shape[0]
concatenated_scores['cancer_type'] = labels
# %%
# Prepare results directory
os.makedirs("pathway_barplots", exist_ok=True)
#%%
for pathway in concatenated_scores.columns[:-1]:  # exclude 'cancer_type'
    df_plot = concatenated_scores[[pathway, 'cancer_type']].copy()
    df_plot = df_plot.rename(columns={'cancer_type': 'group'})
    group_GBM= df_plot[df_plot['group'] == 'GBM']
    group_LGG = df_plot[df_plot['group'] == 'LGG']
    group_BrainMet = df_plot[df_plot['group'] == 'BrainMetastases']
    tGBM_LGG,pGBM_LGG = ttest_ind(group_GBM[pathway], group_LGG[pathway])
    # do ttest between GBM and BrainMet
    tGBM_BrainMet, pGBM_BrainMet = ttest_ind(group_GBM[pathway], group_BrainMet[pathway])
    # do ttest between LGG and BrainMet
    tLGG_BrainMet, pLGG_BrainMet = ttest_ind(group_LGG[pathway], group_BrainMet[pathway])
    print(f"{pathway} - GBM vs LGG: t={tGBM_LGG:.2f}, p={pGBM_LGG:.4f}")
    print(f"{pathway} - GBM vs BrainMet: t={tGBM_BrainMet:.2f}, p={pGBM_BrainMet:.4f}")
    print(f"{pathway} - LGG vs BrainMet: t={tLGG_BrainMet:.2f}, p={pLGG_BrainMet:.4f}")

    #save the ttest results to a file
    with open('ttest_results.txt', 'a') as f:
        f.write(f"{pathway} - GBM vs LGG: t={tGBM_LGG:.2f}, p={pGBM_LGG:.4f}\n")
        f.write(f"{pathway} - GBM vs BrainMet: t={tGBM_BrainMet:.2f}, p={pGBM_BrainMet:.4f}\n")
        f.write(f"{pathway} - LGG vs BrainMet: t={tLGG_BrainMet:.2f}, p={pLGG_BrainMet:.4f}\n")

    plt.figure(figsize=(5, 4))
    sns.barplot(x='group', y=pathway, data=df_plot, capsize=0.05, errorbar='sd', color='grey')
    plt.axhline(0, color='black', linestyle='--', linewidth=1)
    plt.title(pathway)
    plt.ylabel('Z-scored activity')
    plt.xlabel('Cancer type')
    plt.tight_layout()
    plt.savefig(f'pathway_barplots/{pathway}_barplot.png')
    plt.close()