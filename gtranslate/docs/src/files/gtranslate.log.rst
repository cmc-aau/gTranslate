.. _files/gtranslate.log:

gtranslate.log
==============

The console output of GTDB-Tk saved to disk.

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
