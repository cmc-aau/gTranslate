.. _commands/build_features:

build_features
==============

``build_features`` turns raw genomic sequences into numbers.
it calculates the specific signatures—like amino acid frequencies and density metrics the model trains on to distinguish translation tables.


Arguments
---------

.. argparse::
   :module: gtranslate.cli
   :func: get_main_parser
   :prog: gtranslate
   :path: build_features
   :nodefaultconst:

Files output
------------

* :ref:`gtranslate_training.feature_summary.tsv <files/feature_summary.tsv>`
* :ref:`[prefix].log <files/gtranslate.log>`
* :ref:`[prefix].warnings.log <files/gtranslate.warnings.log>`


Example
-------

Input
^^^^^

.. code-block:: bash

    gtranslate build_features --batchfile 1000genomes/5K_batchfile.tsv --out_dir features_test --cpus 90

Output
^^^^^^

.. code-block:: text

    [2026-04-09 23:09:40] INFO: gTranslate v0.0.2
    [2026-04-09 23:09:40] INFO: gtranslate build_features --batchfile 1000genomes/5K_batchfile.tsv --out_dir features_test --cpus 90
    [2026-04-09 23:09:40] INFO: Generating feature vectors for training models.
    [2026-04-09 23:10:45] TASK: Running Prodigal V2.6.3 to identify genes.
    [2026-04-09 23:36:35] INFO: Completed 5,000 genomes in 25.84 minutes (193.47 genomes/minute).
    [2026-04-09 23:36:37] INFO: Done.