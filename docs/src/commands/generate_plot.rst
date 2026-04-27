.. _commands/generate_plot:

generate_plot
=============

Machine learning can sometimes feel like a "black box". ``generate_plot`` lets you see exactly why the model made a certain call. It creates an interactive HTML dashboard where you can see your genomes mapped out. You can visually explore how they cluster based on their Trp/Gly ratios and coding density, making it much easier to verify outliers or double-check tricky predictions.


Arguments
---------

.. argparse::
   :module: gtranslate.cli
   :func: get_main_parser
   :prog: gtranslate
   :path: generate_plot
   :nodefaultconst:


Files output
------------





Example
-------

Input
^^^^^

.. code-block:: bash

    gtranslate generate_plot --feature_file gtranslate_featured.tsv --output_file gtranslate_plot.html

Output
^^^^^^

.. code-block:: text

    [2026-04-27 22:30:33] INFO: gTranslate v0.0.2
    [2026-04-27 22:30:33] INFO: gtranslate generate_plot --feature_file gtranslate_featured.tsv --output_file gtranslate_plot.html
    [2026-04-27 22:30:33] INFO: Running generate_plot
    [2026-04-27 22:30:49] INFO: Interactive HTML dashboard generated: gtranslate_plot.html
    [2026-04-27 22:30:49] INFO: Done.
