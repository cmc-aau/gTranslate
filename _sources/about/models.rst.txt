.. _about/models:

Models
======

gTranslate predicts translation tables using an ensemble of five classifiers. The models used by ``detect_table`` are
provided as a separate download (see :ref:`installing#gtranslate-models`).


Classifiers
-----------

The ensemble was selected by comparing seven classifiers and all combinations of three and five classifiers in
the gTranslate manuscript (see :ref:`about/validation`). It consists of:

.. list-table::
   :widths: 25 25 50
   :header-rows: 1

   * - Classifier
     - Model file
     - Hyperparameters
   * - AdaBoost
     - ``ada_multi_class.pkl.gz``
     - 50 estimators of decision trees with a maximum depth of 3, learning rate of 1
   * - Decision tree
     - ``dt_multi_class.pkl.gz``
     - Gini criterion, balanced class weights
   * - k-nearest neighbours
     - ``knn_multi_class.pkl.gz``
     - 2 neighbours, Manhattan distance, uniform weights
   * - Multi-layer perceptron
     - ``mlp_multi_class.pkl.gz``
     - Hidden layers of 20 and 10 neurons, tanh activation, Adam solver, alpha of 0.0001
   * - XGBoost
     - ``xgb_multi_class.pkl.gz``
     - 50 estimators, learning rate of 0.1, subsample of 0.5

Hyperparameter values were selected by grid search with 5-fold cross-validation. Features are scaled to between
0 and 1 (min-max scaling) for the k-nearest neighbours and multi-layer perceptron classifiers. The classifiers are
implemented with scikit-learn and XGBoost, and the file ``label_encoder.pkl.gz`` maps the classifier outputs to
translation tables.


Training data
-------------

The models distributed with gTranslate 0.0.4 and later are the models evaluated in the gTranslate manuscript (see
:ref:`about/validation`). They were trained on genomes from GTDB release R09-RS220, dereplicated
to at most 100 genomes per species to avoid biasing the models towards species with many sequenced genomes. This gave
356,020 genomes, of which 354,299 (99.52%) use table 11, 1,521 (0.43%) use table 4, and 200 (0.06%) use table 25.

As the GTDB is updated, the models can be retrained to include new genomes and taxa, including newly identified
reassignments of the UGA stop codon.


Ground truth
------------

The expected translation table of each training genome is assigned from its taxonomic classification by the
``ground_truth`` command:

* **Table 25:** GTDB orders Absconditabacterales and BD1-5 (class JAEDAM01).
* **Table 4:**

  * GTDB order Mycoplasmatales, species *Zinderia insecticola*, and *Eggerthellaceae* genera *Equivita* and
    *Gorillivita*.
  * GTDB *Minisyncoccia* family JAKLIH01 and four species of the genus GCA-2747955, in which a UGA reassignment was
    identified in the gTranslate manuscript.
  * NCBI species *Ca.* Hodgkinia cicadicola, *Ca.* Nasuia deltocephalinicola, *Ca.* Stammera capleta,
    *Ca.* Organicella extenuata and *Ca.* Pinguicoccus supinus, which are not in the GTDB as their highly reduced
    genomes cannot be robustly placed in the GTDB reference trees.

* **Unresolved:** GTDB species *Providencia_A siddallii* and *Providencia_A siddallii_A*, whose translation table is
  unclear (`Manzano-Marín et al., 2023 <https://doi.org/10.1093/gbe/evad164>`_). These genomes are excluded from
  training.
* **Table 11:** all other genomes.

The *Minisyncoccia* and *Providencia_A* rules were added after the analyses in the gTranslate manuscript, so the
ground truth produced by the ``ground_truth`` command differs slightly from the one used to train the distributed
models.
The ground truth for individual genomes can be set with the ``--manual_gt_file`` argument.


Training your own models
------------------------

The models can be retrained, for example to include genomes from a newer GTDB release:

1. :ref:`commands/ground_truth`: assign the expected translation table to each genome from its taxonomy.
2. :ref:`commands/build_features`: call genes and calculate the features for each genome.
3. :ref:`commands/fit_models`: train the five classifiers using the features and ground truth.

The resulting directory can be used with ``detect_table`` through the ``--custom_model_path`` argument or
the ``GTRANSLATE_MODEL_PATH`` environment variable.
