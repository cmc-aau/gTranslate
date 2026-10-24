.. _about/validation:

Validation
==========

The results below are from the gTranslate manuscript
(`Chaumeil et al., 2026 <https://doi.org/10.7717/peerj.21825>`_), where the full results are available. They apply to
the models distributed with gTranslate 0.0.4 and later. Models were trained on 356,020 genomes from GTDB R09-RS220 (see :ref:`about/models`), and the expected translation
table of each genome was established from its taxonomic classification. gTranslate was compared to the coding
density heuristic used by CheckM, `Codetta <https://doi.org/10.7554/eLife.71402>`_ (Shulgina & Eddy, 2021), and
the translation table information provided by NCBI.


Cross-validation
----------------

Using 5-fold cross-validation across the 356,020 GTDB R09-RS220 genomes, the best individual classifier
(k-nearest neighbours) made 15 misclassifications and the best ensembles 14. By comparison, the CheckM heuristic
made 101 misclassifications, and Codetta made 444 or 78 when requiring 8 or 34 alignment columns, respectively, to
support a UGA reassignment. NCBI translation table information disagreed with the expected translation table for 62
of the 207,495 genomes for which it was available.


Taxonomically novel genomes
---------------------------

To assess performance on genomes from taxa absent from the training data, all genomes from a single species, genus
or family were removed, the classifiers were trained on the remaining genomes, and the removed genomes were
classified. Performance decreased with taxonomic novelty, and ensembles of classifiers performed more consistently
than individual classifiers. The ensemble used by gTranslate (AdaBoost, decision tree, k-nearest neighbours,
multi-layer perceptron and XGBoost) was selected as it had the fewest misclassifications across the cross-validation
and leave-one-taxon-out tests.


Independent test sets
---------------------

.. list-table::
   :widths: 40 20 20 20
   :header-rows: 1

   * - Test set
     - Genomes
     - Correctly predicted
     - Accuracy
   * - New genomes in GTDB R10-RS226
     - 116,508
     - 116,502
     - 99.995%
   * - GlobDB release 226
     - 134,075
     - 134,030
     - 99.97%

These test sets include genomes from many taxa absent from the training data, and the GlobDB genomes come from
metagenomic studies not represented in the GTDB.

Most disagreements appear to be errors in the expected translation table rather than in the predictions. For
example, five *Ca.* Stammera capleta genomes predicted to use table 11 form a basal lineage of this species with
genomic evidence of using table 11, and 22 *Minisyncoccia* genomes predicted to use table 4 have genomic evidence of
a UGA reassignment to tryptophan. After accounting for such cases, gTranslate made 16 misclassifications across the
250,583 test genomes (99.994% accuracy), compared to 25 for the CheckM heuristic and 34 for Codetta. Twelve of the 16
misclassifications were between tables 4 and 25, giving an accuracy of 98.9% across the 1,077 genomes not using
table 11. Unlike the CheckM heuristic, gTranslate distinguishes between tables 4 and 25.


Computational performance
-------------------------

gTranslate requires approximately 20 seconds per genome per CPU, approximately 100 times less than Codetta.
