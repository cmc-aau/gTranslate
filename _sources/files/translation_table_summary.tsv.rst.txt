.. _files/translation_table_summary.tsv:


translation_table_summary.tsv
========================================

The output of gTranslate is provided in a tab-separated file (typically gtranslate.tsv). This file contains the final taxonomic predictions, the biological features calculated from the sequences, and the specific outputs of the machine learning ensemble.

* user_genome: Unique identifier of the query genome.
* best_tln_table: The final predicted Genetic Translation Table (GTT). This is the "consensus" choice recommended by the tool (e.g., 11 for standard, 4 for UGA=Trp, or 25 for UGA=Gly).
* coding_density_4: The calculated coding density of the genome assuming the translation table is Table 4.
* coding_density_11: The calculated coding density assuming the translation table is Table 11 (Standard code).
* gc_percent: The percentage of Guanine and Cytosine in the genome.
* n50: N50 of the genome assembly.
* genome_size: The total length of the genome assembly in base pairs.
* contig_count: The total number of contigs (fragments) in the user's genome file.
* confidence: A score ranging from 0.0 to 1.0 representing the level of agreement across the internal machine learning models. A score of 1.0 indicates that all classifiers in the ensemble (Adaboost, MLP, etc.) predicted the same table.
* adaboost_pred: The specific GTT prediction made by the AdaBoost classifier.
* decisiontree_pred: The specific GTT prediction made by the Decision Tree classifier.
* kneighbors_pred: The specific GTT prediction made by the K-Nearest Neighbors (KNN) classifier.
* xgb_pred: The specific GTT prediction made by the XGBoost (Extreme Gradient Boosting) classifier.
* mlp_pred: The specific GTT prediction made by the Multi-Layer Perceptron (Neural Network) classifier.
* warnings: Flags any unusual features or inconsistencies, such as extreme coding density differences or low confidence scores, that may suggest the prediction should be manually reviewed. "N/A" indicates no issues were detected.


Produced by
-----------

* :ref:`commands/detect_table`

Example
-------

.. code-block:: text

    user_genome	best_tln_table	coding_density_4	coding_density_11	gc_percent	n50	genome_size	contig_count	confidence	adaboost_pred	decisiontree_predkneighbors_pred	xgb_pred	mlp_pred	warnings
    GCA_000145705.TT4	4	90.29257	64.14511	25.88353	839615	839615	1	1.0	4	4	4	4	4	N/A
    GCA_027046965.1_TT25	25	94.21967	56.0663	32.96163	50860	854813	32	1.0	25	25	25	25	25	N/A
    GCA_910575315.1_TT11	11	87.75128	87.91231	44.09593	73133	4990987	224	1.0	11	11	11	11	11	N/A