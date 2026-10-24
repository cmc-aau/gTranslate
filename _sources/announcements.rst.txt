Announcements
=============

gTranslate 0.0.6
----------------

*October, 2026*

* Fixed the genes reported for gzipped genomes predicted to use translation table 25 when
  ``--keep_called_genes`` is used. Predicted translation tables were not affected.
* ``check_install`` now also verifies that the gTranslate models are present.

gTranslate 0.0.5
----------------

*September, 2026*

* gTranslate no longer exits when Prodigal fails to call genes for a genome. Such genomes
  are skipped and reported in ``gtranslate.log``, ``gtranslate.warnings.log``, and the
  ``failed_genomes.tsv`` report.

gTranslate 0.0.4
----------------

*May, 2026*

* The new version of gTranslate does not include the models used in the previous version. Users will need to download the new models and set the environment variable ``GTRANSLATE_MODEL_PATH`` to point to the new models.
