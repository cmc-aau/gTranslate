.. _commands/ground_truth:

ground_truth
============

Before training, you need to know the "right" answers. This command looks at the taxonomy of your genomes to assign them their confirmed genetic codes. This creates the labeled dataset that the model uses to learn.


Arguments
---------

.. argparse::
   :module: gtranslate.cli
   :func: get_main_parser
   :prog: gtranslate
   :path: ground_truth
   :nodefaultconst:

Files output
------------

* :ref:`[ground_truth_results].tsv <files/ground_truth_results.tsv>`


Example
-------

Input
^^^^^

.. code-block:: bash

    gtranslate ground_truth --taxonomy_file taxonomy_file_r226.tsv.gz --output_file ground_truth_results.tsv

Output
^^^^^^

.. code-block:: text

    [2026-05-01 01:50:19] INFO: gTranslate v0.0.2
    [2026-05-01 01:50:19] INFO: gtranslate ground_truth --taxonomy_file taxonomy_file_r226.tsv.gz --output_file ground_truth_results.tsv
    [2026-05-01 01:50:19] INFO: Selecting Ground Truth translation tables based on taxonomic classification.
    [2026-05-01 01:50:19] INFO: Determining ground truth for genomes:
    [2026-05-01 01:50:20] INFO:  - determined ground truth for 116,508 genomes
    [2026-05-01 01:50:20] INFO: Table 11: 115,941 (99.51%)
    [2026-05-01 01:50:20] INFO: Table 25: 121 (0.10%)
    [2026-05-01 01:50:20] INFO: Table 4: 445 (0.38%)
    [2026-05-01 01:50:20] INFO: Table UNRESOLVED: 1 (0.00%)
    [2026-05-01 01:50:20] INFO: Done.