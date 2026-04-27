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