.. _commands/fit_models:

fit_models
==========

``fit_models`` takes the features built from ``build_features`` and trains the model to recognize the boundaries between Table 11, 4, and 25.
Once it’s done, you’ll have fresh model files ready for the ``detect_table`` command.


Arguments
---------

.. argparse::
   :module: gtranslate.cli
   :func: get_main_parser
   :prog: gtranslate
   :path: fit_models
   :nodefaultconst: