.. _commands/detect_table:

detect_table
============

Detect the correct translation table for a given sequence file.

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

* output_dir
    * :ref:`[prefix].translation_table_summary.tsv <files/gtranslate.translation_table_summary.tsv>`
    * :ref:`[prefix].log <files/gtranslate.log>`
    * :ref:`[prefix].warnings.log <files/gtranslate.warnings.log>`
    * predict
        * :ref:`[prefix].failed_genomes.tsv <files/failed_genomes.tsv>`
        * :ref:`[prefix].translation_table_summary.tsv <files/gtranslate.translation_table_summary.tsv>`
        * intermediate_results
            * called_genes
                * <individual genome files id>
                    * :ref:`<genome_id>.protein.faa <files/genome_id.protein.faa>`
                    * :ref:`<genome_id>.protein.fna <files/genome_id.protein.fna>`
                    * :ref:`<genome_id>.protein.gff <files/genome_id.protein.gff>`
                    * :ref:`<genome_id>_translation_table.tsv <files/genome_id_translation_table.tsv>`

Example
-------


Input
^^^^^


.. code-block:: bash

     gtranslate detect_table --genome_dir gtranslate_tests/genomes/ --out_dir gtranslate_tests/test



Output
^^^^^^

.. code-block:: text

    [2025-03-13 15:39:36] INFO: gTranslate v0.0.1
    [2025-03-13 15:39:36] INFO: gtranslate detect_table --genome_dir gtranslate_tests/genomes/ --out_dir gtranslate_tests/test
    [2025-03-13 15:39:36] INFO: Running detect_table
    [2025-03-13 15:39:36] INFO: All dependencies are correctly installed.
    [2025-03-13 15:39:36] TASK: Running Prodigal V2.6.3 to identify genes.
    [2025-03-13 15:39:36] TASK: Running Jellyfish V2.3.1 to calculate kmers.
    [2025-03-13 15:40:15] INFO: Completed 4 genomes in 39.16 seconds (9.79 seconds/genome).
    [2025-03-13 15:40:15] TASK: Summarising translation tables.
    [2025-03-13 15:40:15] INFO: Completed 4 genomes in 0.00 seconds (45,221.61 genomes/second).
    [2025-03-13 15:40:15] INFO: Cleaning up intermediate files.
    [2025-03-13 15:40:15] INFO: Deleting generated gene files.
    [2025-03-13 15:40:15] INFO: gene files deleted.
    [2025-03-13 15:40:15] INFO: Done.