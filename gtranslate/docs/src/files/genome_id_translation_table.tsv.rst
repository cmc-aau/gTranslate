.. _files/genome_id_translation_table.tsv:

genome_id_translation_table.tsv
===============================

Individual file storing information about the translation table used for each genome.

.. list-table::
   :header-rows: 1

   * - Column Name
     - Description
    * - **best_translation_table**
        - The best translation table detected by gTranslate.
    * - **coding_density_4**
        - The coding density of the genome using translation table 4.
    * - **coding_density_11**
        - The coding density of the genome using translation table 11.
    * - **gc_percent**
        - The GC content of the genome.
    * - **n50**
        - The N50 value of the genome.
    * - **genome_size**
        - The size of the genome in base pairs.
    * - **contig_count**
        - The number of contigs in the genome.
    * - **probability_4_11**
        - The highest probability of the genome being translated by either translation table 4 or 11.
    * - **probability_4_25**
        - The highest probability of the genome being translated by either translation table 4 or 25 (if the first classifier predicts translation table 4).

Produced by
-----------

* :ref:`commands/detect_table`




Example
-------

.. code-block:: text

    [2025-03-17 16:07:53] INFO: gTranslate v0.0.1
    [2025-03-17 16:07:53] INFO: gtranslate detect_table --genome_dir test_command/genomes --out_dir test_command/output --cpus 1
    [2025-03-17 16:07:53] INFO: Running detect_table
    [2025-03-17 16:07:53] TASK: Running Prodigal V2.6.3 to identify genes.
    [2025-03-17 16:07:53] TASK: Running Jellyfish V2.3.1 to calculate kmers.
    [2025-03-17 16:08:22] INFO: Completed 3 genomes in 29.13 seconds (9.71 seconds/genome).
    [2025-03-17 16:08:22] TASK: Summarising translation tables.
    [2025-03-17 16:08:22] INFO: Completed 3 genomes in 0.00 seconds (36,792.14 genomes/second).
    [2025-03-17 16:08:22] INFO: Cleaning up intermediate files.
    [2025-03-17 16:08:22] INFO: Deleting generated gene files.
    [2025-03-17 16:08:22] INFO: gene files deleted.
    [2025-03-17 16:08:22] INFO: Done.