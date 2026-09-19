Announcements
=============

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
