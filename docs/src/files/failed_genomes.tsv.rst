.. _files/failed_genomes.tsv:

failed_genomes.tsv
==================

Tab-separated file listing the genomes which were skipped, with the reason they were skipped. It is written to
``predict/[prefix].failed_genomes.tsv`` (default prefix: ``gtranslate``) on every run, and is empty if no genomes
were skipped. Skipped genomes are not included in the
:ref:`translation_table_summary.tsv <files/translation_table_summary.tsv>` and
:ref:`feature_summary.tsv <files/feature_summary.tsv>` files.

The file has no header and two columns: the genome identifier and the reason the genome was skipped, e.g.:

* ``Cannot call Prodigal on an empty genome``
* ``No genes were called by Prodigal``
* ``No glycine codons in the genes called by Prodigal, so the genome could not be classified``
* the error message reported by Prodigal, if Prodigal failed


Produced by
-----------

* :ref:`commands/build_features`
* :ref:`commands/detect_table`

Example
-------

.. code-block:: text

    genome_A	No genes were called by Prodigal
    genome_B	Cannot call Prodigal on an empty genome
