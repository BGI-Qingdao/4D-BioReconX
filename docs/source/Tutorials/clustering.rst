.. _`clustering`:
========================================
Spatial pattern Clustering
========================================

Perform clustering via `hdbscan` and `LogicRegression`

Usage
---------------------------------

.. code-block:: python3

  from pcg_pattern import *
  from sklearn.linear_model import LogisticRegression
  import hdbscan 

Load data in anndata format
++++++++++++++++++++++++++++++++++++

.. code-block:: python3

  adata = sc.read_h5ad(file_name)
  df = adata.to_df()  # ensure rows are genes/obs and columns are bins/features, if not, transpose matrix first

First round of clustering using hdbscan
++++++++++++++++++++++++++++++++++++

.. code-block:: python3

  clusterer = hdbscan.HDBSCAN()
  clusterer.fit(df.to_numpy())
  labels = clusterer.labels_
  labels = labels.astype(int)

Second round of "clustering"
++++++++++++++++++++++++++++++++++++

for those genes that were not assigned into a cluster in the first round (labeled as -1), 
use LogicRegression to find which cluster's genes it is most similar to

.. code-block:: python3

  svc = LogisticRegression()
  svc.fit(df.to_numpy(),labels[labels!=-1])
  recall_non_labels = svc.predict(df[labels==-1].to_numpy())

Save results
++++++++++++++++++++++++++++++++++++

Visualization by heatmap
++++++++++++++++++++++++++++++++++++

.. code-block:: python3

  plot_heatmap(df.to_numpy(), labels)
