.. _faq:

FAQ
===


Which translation tables can gTranslate predict?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Translation tables 11 (the standard bacterial code), 4 (UGA reassigned to tryptophan) and 25 (UGA reassigned to
glycine). Other genetic codes, such as the reassignment of UAG to pyrrolysine, are not predicted. gTranslate is
designed for prokaryotic genomes and should not be used for eukaryotic or organelle genomes.


How accurate is gTranslate?
~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the gTranslate manuscript, the distributed models correctly predicted the translation table of more than 99.99% of
250,583 test genomes, including many from taxa absent from the training data. Most errors were between tables 4 and
25. See :ref:`about/validation`.


How should I interpret the confidence and warnings?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A confidence of 1.0 means all five classifiers predicted the same translation table. Lower values indicate
disagreement between the classifiers, and a warning is reported when no table received more than two votes.
Predictions with a low confidence or warnings should be reviewed, e.g. by checking the predictions for related
genomes or using :ref:`commands/generate_plot` to compare the genome to others. See
:ref:`files/translation_table_summary.tsv` for how the confidence is calculated and the possible warnings.


Can I use gTranslate on metagenome-assembled genomes (MAGs) or incomplete genomes?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Yes. The training data included genomes as low as 18% complete, and gTranslate was tested on genomes from metagenomic
studies that were at least 50% complete. Genomes under 100 kb are processed using Prodigal's metagenome mode.
Predictions for very small or highly fragmented genomes should be treated with more caution.


Why does the prediction disagree with the taxonomy of my genome, or with NCBI?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The translation table is not always the same for all members of a taxon. For example, the gTranslate manuscript
identified a lineage of *Ca.* Stammera capleta using table 11 although the species uses table 4, and lineages of
*Minisyncoccia* using table 4 within a phylum otherwise using tables 11 and 25. Disagreements can also arise from
genomes that are misclassified. The translation table information at NCBI is based on taxonomy and was incorrect for
some genomes in the manuscript's analyses.


Why is a genome missing from the results?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Genomes for which Prodigal could not call genes, e.g. empty genomes, are skipped. Skipped genomes and the reason
they were skipped are listed in :ref:`files/failed_genomes.tsv` and ``gtranslate.warnings.log``.


Can I use gzipped genomes?
~~~~~~~~~~~~~~~~~~~~~~~~~~

Yes. With ``--genome_dir``, set ``--extension`` to the full extension of the files, e.g. ``-x fna.gz``, so that the
genome identifiers do not include the ``.fna`` extension. With ``--batchfile``, gzipped and uncompressed files can be
mixed.


gTranslate reports missing model files, what should I do?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``GTRANSLATE_MODEL_PATH`` environment variable, or the ``--custom_model_path`` argument, must point to the
directory containing the six extracted model files. Run :ref:`commands/check_install` to check which files are
missing, and see :ref:`installing#gtranslate-models` to download the models.


Do I need to download new models when upgrading gTranslate?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Not currently: gTranslate 0.0.4 and later use the same models. Any release that requires new models will be listed
in the :doc:`announcements`.


How long does gTranslate take to run?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Approximately 20 to 30 seconds per genome per CPU, most of which is spent calling genes with Prodigal. Use ``--cpus``
to process genomes in parallel.


Can I train my own models?
~~~~~~~~~~~~~~~~~~~~~~~~~~

Yes, using the :ref:`commands/ground_truth`, :ref:`commands/build_features` and :ref:`commands/fit_models` commands.
See :ref:`about/models` for how the distributed models were trained.


How do I cite gTranslate?
~~~~~~~~~~~~~~~~~~~~~~~~~

Please cite the gTranslate manuscript and the third-party software it uses, see :ref:`references`.
