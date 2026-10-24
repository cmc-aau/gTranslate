.. _files/gtranslate.warnings.log:

gtranslate.warnings.log
=======================

Detailed warnings encountered while running gTranslate, such as each genome that was skipped and the reason it
was skipped. If warnings are written to this file, the console output will indicate that it should be checked.
It is written to the output directory as ``gtranslate.warnings.log`` (regardless of ``--prefix``) by commands with
an ``--out_dir`` argument, and output from subsequent runs is appended to the file.


Produced by
-----------

* :ref:`commands/build_features`
* :ref:`commands/detect_table`
* :ref:`commands/fit_models`
* :ref:`commands/test`
