.. _`mining`:
========================================
PCG Mining Pipeline
========================================

.. note:: 
Find genes that could be PCGs by calculating and evaluating correlation coefficient values between known PCGs (from published studies) and target genes.

A. Use Spearman's rank correlation coefficient
===============================================

.. code-block:: python3
import scipy.stats
import scanpy as sc
import pandas as pd

.. code-block:: python3
adata = sc.read_h5ad(file_name)
df = adata.to_df()
target_genes = list(df.columns)  # or assigned by users
with open('known_pcgs.txt', 'r') as f:
    known_pcgs = f.read().splitlines()

Spearman's rank correlation coefficient
==============================================

.. code-block:: python3
df_data = {'gene':[], 'pcg':[], 'scc':[]}
for gene in target_genes:
    for pcg in known_pcgs:
        if pcg in df.columns and pcg != gene:
            df_data['gene'].append(gene)
            df_data['pcg'].append(pcg)
            scc = scipy.stats.spearmanr(df[gene], df[pcg]).correlation
            df_data['scc'].append(scc)

.. code-block:: python3
results = pd.DataFrame(df_data)
results = results.sort_values(by='scc', ascending=False)
results.to_csv('spearman.csv', index=False)
