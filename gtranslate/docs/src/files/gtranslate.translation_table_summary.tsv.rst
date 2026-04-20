.. _files/gtranslate.translation_table_summary.tsv:


gtranslate.translation_table_summary.tsv
=========================================

.. list-table::
   :header-rows: 1

   * - Column Name
     - Description
   * - **user_genome**
     - Unique identifier of query genome taken from the FASTA file of the genome.
   * - **best_tln_table**
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

.. list-table:: Translation Table Summary File
   :header-rows: 1
   :widths: auto

   * - user_genome
     - best_tln_table
     - coding_density_4
     - coding_density_11
     - gc_percent
     - n50
     - genome_size
     - contig_count
     - probability_4_11
     - probability_4_25
   * - GCA_000145705.TT4
     - 4
     - 90.29257
     - 64.14511
     - 25.88353
     - 839615
     - 839615
     - 1
     - 0.70856
     - 0.77237
   * - GCA_027046965.1_TT25
     - 25
     - 94.21967
     - 56.0663
     - 32.96163
     - 50860
     - 854813
     - 32
     - 0.68575
     - 0.60195
   * - GCA_910575315.1_TT11
     - 11
     - 87.75128
     - 87.91231
     - 44.09593
     - 73133
     - 4990987
     - 224
     - 0.80575
     - N/A
