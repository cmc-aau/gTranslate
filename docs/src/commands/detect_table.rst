.. _commands/detect_table:

detect_table
============

The ``detect_table`` command is the main command of gtranslate.
It is used to predict the translation table (TT) used by prokaryotic genomes. By analyzing specific sequence features — such as coding density differences, Trp ratios, and Gly ratios — `gTranslate` can accurately distinguish between the genetic codes associated with reassignment of the UGA stop codon, i.e. translation tables 11 (standard prokaryotic code), 4 (UGA=Trp), and 25 (UGA=Gly).

Arguments
---------

.. argparse::
   :module: gtranslate.cli
   :func: get_main_parser
   :prog: gtranslate
   :path: detect_table
   :nodefaultconst:


Files output
------------

* :ref:`[prefix].translation_table_summary.tsv <files/translation_table_summary.tsv>`
* :ref:`[prefix].log <files/gtranslate.log>`
* :ref:`[prefix].warnings.log <files/gtranslate.warnings.log>`


Example
-------

Input
^^^^^

.. code-block:: bash

    gtranslate detect_table --genome_dir test_command/genomes --out_dir test_command/output_multi --cpus 3

Output
^^^^^^

.. code-block:: text

    [2026-03-05 20:36:09] INFO: gTranslate v0.0.2
    [2026-03-05 20:36:09] INFO: gtranslate detect_table --genome_dir test_command/genomes --out_dir test_command/output_multi --cpus 3
    [2026-03-05 20:36:09] INFO: Running detect_table
    [2026-03-05 20:36:10] TASK: Running Prodigal V2.6.3 to identify genes.
    [2026-03-05 20:36:27] INFO: Completed 3 genomes in 17.09 seconds (5.70 seconds/genome).
    [2026-03-05 20:36:27] TASK: Summarising translation tables.
    [2026-03-05 20:36:27] INFO: Completed 3 genomes in 0.00 seconds (31,615.36 genomes/second).
    [2026-03-05 20:36:27] INFO: Cleaning up intermediate files.
    [2026-03-05 20:36:27] INFO: Deleting generated gene files.
    [2026-03-05 20:36:27] INFO: gene files deleted.
    [2026-03-05 20:36:27] INFO: Done.
