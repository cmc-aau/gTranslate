
Change log
==========

0.0.7
-----

Bug Fixes:

* The genomes used by the ``test`` command are now included in the PyPI and Bioconda packages.
  Previously, ``gtranslate test`` failed on installed copies of gTranslate as these genomes were
  missing.


0.0.6
-----

Bug Fixes:

* Genes reported for genomes predicted to use translation table 25 are now called on the
  uncompressed genome. Previously, gzipped genomes predicted to use table 25 had their genes
  called on the compressed file, producing only a handful of spurious genes in the files kept
  with ``--keep_called_genes``. Predicted translation tables were not affected.

Minor Changes:

* ``check_install`` now verifies that ``GTRANSLATE_MODEL_PATH`` is set and that all model files
  are present.


0.0.5
-----

Bug Fixes:

* gTranslate no longer exits when Prodigal fails to call genes for a genome. The genome is
  now skipped, the number of skipped genomes is reported in ``gtranslate.log``, and each
  skipped genome is listed with the reason it failed in ``gtranslate.warnings.log`` and in
  the ``failed_genomes.tsv`` report.
* Errors reported by Prodigal are written to the log files instead of standard output.
* Fixed a division by zero error which occurred when no genes were called for a genome.
* Fixed the reporting of the Prodigal ``meta`` mode fallback, which duplicated earlier
  warnings in the translation table summary.


0.0.4
-----

Major Changes:

* version 0.4.0 of gTranslate does not include the models used in the previous version. Users will need to download the new models and set the environment variable ``GTRANSLATE_MODEL_PATH`` to point to the new models.


0.0.2
-----

in Dev

