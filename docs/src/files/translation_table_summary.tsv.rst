.. _files/translation_table_summary.tsv:

translation_table_summary.tsv
=============================

The main output of gTranslate: a tab-separated file containing the predicted translation table of each genome,
summary statistics of the genome, and the prediction of each classifier in the ensemble.

It is written to ``predict/[prefix].translation_table_summary.tsv`` (default prefix: ``gtranslate``), with a
symbolic link of the same name in the output directory. Genomes skipped because Prodigal could not call genes are
not included (see :ref:`failed_genomes.tsv <files/failed_genomes.tsv>`).

The columns in this file are as follows:

* ``user_genome``: unique identifier of the genome.
* ``best_tln_table``: the predicted translation table: 11 (standard bacterial code), 4 (UGA = tryptophan) or
  25 (UGA = glycine).
* ``coding_density_4``: gene coding density (%) when predicting genes with Prodigal using translation table 4.
* ``coding_density_11``: gene coding density (%) when predicting genes with Prodigal using translation table 11.
* ``gc_percent``: percentage of guanine (G) and cytosine (C) nucleotides in the genome.
* ``n50``: N50 of the genome assembly.
* ``genome_size``: total length of the genome assembly in base pairs.
* ``contig_count``: number of contigs in the genome assembly.
* ``confidence``: support for the prediction, see :ref:`below <files/translation_table_summary.tsv#confidence>`.
* ``adaboost_pred``: translation table predicted by the AdaBoost classifier.
* ``decisiontree_pred``: translation table predicted by the decision tree classifier.
* ``kneighbors_pred``: translation table predicted by the k-nearest neighbours classifier.
* ``xgb_pred``: translation table predicted by the XGBoost classifier.
* ``mlp_pred``: translation table predicted by the multi-layer perceptron classifier.
* ``warnings``: semicolon-separated warnings for the genome, or ``N/A`` if there are none, see
  :ref:`below <files/translation_table_summary.tsv#warnings>`.


.. _files/translation_table_summary.tsv#confidence:

Prediction and confidence
-------------------------

The translation table is predicted by majority vote of the five classifiers, and ``confidence`` is the
fraction of classifiers that voted for the predicted table. When table 4 or 25 is predicted, 0.05 is added for
each vote for the other UGA reassignment table, as these votes still support a reassignment of the UGA stop codon.
For example, four votes for table 4 and one vote for table 25 gives a confidence of 0.85.

If no table receives more than two votes, gTranslate predicts a reassignment of the UGA stop codon: table 4
or 25, whichever received more votes, or table 4 if both received two votes. For example, votes of 11, 11, 4, 4,
and 25 predict table 4 with a confidence of 0.45.


.. _files/translation_table_summary.tsv#warnings:

Warnings
--------

The ``warnings`` column can contain the following:

* ``Low confidence: maximum model agreement was only 2/5.``: no table received more than two votes.
* ``Tie broken: ...``: how the prediction was made when no table received more than two votes (see above).
* ``Used Prodigal 'meta' mode fallback for TT4`` (or ``TT11``, ``TT25``): Prodigal failed to train on the
  genome and genes were called using its metagenome mode. This can occur for small or fragmented genomes.
  Genomes under 100 kb always use metagenome mode and do not receive this warning.

Predictions with warnings should be reviewed manually.


Produced by
-----------

* :ref:`commands/detect_table`

Example
-------

.. code-block:: text

    user_genome	best_tln_table	coding_density_4	coding_density_11	gc_percent	n50	genome_size	contig_count	confidence	adaboost_pred	decisiontree_pred	kneighbors_pred	xgb_pred	mlp_pred	warnings
    GCA_000145705.TT4	4	90.29257	64.14511	25.88353	839615	839615	1	1.0	4	4	4	4	4	N/A
    GCA_027046965.1_TT25	25	94.21967	56.0663	32.96163	50860	854813	32	1.0	25	25	25	25	25	N/A
    GCA_910575315.1_TT11	11	87.75128	87.91231	44.09593	73133	4990987	224	1.0	11	11	11	11	11	N/A
