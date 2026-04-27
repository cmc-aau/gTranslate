.. _commands/build_features:

build_features
==============

``build_features`` turns raw genomic sequences into numbers.
it calculates the specific signatures—like amino acid frequencies and density metrics the model trains on to distinguish translation tables.


Arguments
---------

.. argparse::
   :module: gtranslate.cli
   :func: get_main_parser
   :prog: gtranslate
   :path: build_features
   :nodefaultconst: