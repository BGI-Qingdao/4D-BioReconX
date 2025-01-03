.. _`pca`:
========================================
PCG PCA
========================================

PCG analysis in PCG PC space

usage
---------------------------------

.. code-block:: python3

  python3 pca.py -h

============================================= =============================================================================
argument                                      description
============================================= =============================================================================
-h, --help                                     show this help message and exit
--name NAME, -n NAME                           sample name
--symbol SYMBOL, -s SYMBOL                     symbol
--trainingData TRAININGDATA, -t TRAININGDATA   training dataset, cols are genes, rows are bins                               
--pcg PCG, -p PCG                              known pcg list file
--genes GENES, -g GENES                        target gene list file               
--classes CLASSES, -c CLASSES                  pcg class file
--data DATA, -i DATA                           data expression file to project in PC space, cols are genes, rows are bins
============================================= =============================================================================

example
---------------------------------

.. code-block:: python3

  python3 pca.py -i total.csv -t WT.csv -p pcg.txt -c pcg_class.txt -n test_WT

.. code-block:: python3

  python3 pca.py -i total.csv -t WT.csv -g filtered_gene.txt -p pcg.txt -c pcg_class.txt -n test_WT -s symbols.txt
