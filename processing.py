#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from sklearn.decomposition import PCA
#%%
df1 = pd.read_csv('GBM.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt', sep= "\t", index_col=0)
df2 = pd.read_csv('LGG.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt', sep= "\t", index_col=0)
df3 = pd.read_csv('BrainMet.txt', sep= "\t", index_col=0, skiprows=5)
df1 = df1.iloc[30:]
df2 = df2.iloc[30:]
# %%
df1.index = df1.index.str.split('|').str[0]
df2.index = df2.index.str.split('|').str[0]
df1.index.name = 'gene_name'
df2.index.name = 'gene_name'
# %%
# Keep only the common genes in the three datasets
df1 = df1[df1.index.isin(df3.index)]
df2 = df2[df2.index.isin(df1.index)]
df3 = df3[df3.index.isin(df1.index)]
#%%
df1.to_csv('GBM_filtered.tsv', sep="\t")
df2.to_csv('LGG_filtered.tsv', sep="\t")
df3.to_csv('BrainMet_filtered.tsv', sep="\t")
# %%
#Load the AUCell results
df1_scores = pd.read_csv('AUCell_scores_GBM.txt', sep="\t", index_col=0)
df2_scores = pd.read_csv('AUCell_scores_LGG.txt', sep="\t", index_col=0)
df3_scores = pd.read_csv('AUCell_scores_BrainMet.txt', sep="\t", index_col=0)
# %%
df3a = pd.read_csv('cancer_type.txt', sep="\t", index_col=0).T
df3_scores = df3_scores.T
df3_scores = pd.concat([df3_scores, df3a], axis=1)
df2_scores = df2_scores.T
df2_scores['cancer_type'] = 'LGG'
df1_scores = df1_scores.T
df1_scores['cancer_type'] = 'GBM'
#%%
df = pd.concat([df1_scores.T, df2_scores.T, df3_scores.T], axis=1)
# create a label for the three datasets
labels = ['GBM'] * df1_scores.shape[0] + ['LGG'] * df2_scores.shape[0] + ['BrainMetastases'] * df3_scores.shape[0]
# %%
# drop the cancer type row
df_pca = df.drop(['cancer_type'])
# %%
# perform PCA
pca = PCA(n_components=2)
pca_result = pca.fit_transform(df_pca.T)
# create a dataframe with the PCA results and labels and cancer type
pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])
pca_df['label'] = labels
pca_df['cancer_type'] = df.T['cancer_type'].values
# %%
plt.figure(figsize=(10, 8))
plt.rcParams['font.size'] = 16
plt.rcParams['font.family'] = 'Arial'
sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='label', palette=['red', 'orange', 'blue'], s=20)
plt.title('PCA of primary and metastatic brain tumors')
plt.xlabel('PC1 (54.5% variance)')
plt.ylabel('PC2 (11.6% variance)')
plt.legend(title='Datasets')
plt.savefig('PCA_Primary_BrainMet.png', dpi=600, bbox_inches='tight')
plt.show()
# %%
# plot the PCA results
plt.figure(figsize=(10, 8))
sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='cancer_type', palette= ['red', 'orange','green','purple', 'dimgrey','saddlebrown','magenta','olive','cyan','pink','mediumslateblue'], s=30)
plt.title('PCA of primary and metastatic brain tumors')
plt.xlabel('PC1 (54.5% variance)')
plt.ylabel('PC2 (11.6% variance)')
plt.rcParams['font.size'] = 12.5
plt.legend(title='Dataset')
plt.savefig('PCA_GBMLGGBrainMet_cancertype.png', dpi=600, bbox_inches='tight')
plt.show()
# %%
loadings = pd.DataFrame(pca.components_.T, index=df_pca.index, columns=['PC1', 'PC2'])
#%%
# Sort features by loading values for PC1
plt.rcParams['font.size'] = 14
plt.rcParams['font.family'] = 'Arial'
loadings_sorted_PC1 = loadings.reindex(loadings['PC1'].sort_values(ascending=False).index)
plt.figure(figsize=(12.5, 10))
sns.barplot(x=loadings_sorted_PC1.index, y=loadings_sorted_PC1['PC1'], palette="viridis")
plt.title(f"PCA Loadings for PC1")
plt.xticks(rotation=90)
plt.xlabel("Gene Sets")
plt.ylabel("Loading")
plt.tight_layout()
plt.savefig("PCA_Loadings_PC1.png", dpi=600, bbox_inches='tight')
plt.show()
# %%
# Sort features by loading values for PC2
loadings_sorted_PC2 = loadings.reindex(loadings['PC2'].sort_values(ascending=False).index)
plt.figure(figsize=(12.5, 10))
sns.barplot(x=loadings_sorted_PC2.index, y=loadings_sorted_PC2['PC2'], palette="viridis")
plt.title(f"PCA Loadings for PC2")
plt.xticks(rotation=90)
plt.xlabel("Gene Sets")
plt.ylabel("Loading")
plt.tight_layout()
plt.savefig("PCA_Loadings_PC2.png", dpi=600, bbox_inches='tight')
plt.show()
#%%
# write the PC1 and PC2 loadings to a file
loadings.to_csv('PCA_loadings.tsv', sep="\t")
#%%
#write the PCA results to a file
pca_df.to_csv('PCA_results.tsv', sep="\t")