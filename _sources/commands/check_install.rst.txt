.. _commands/check_install:

check_install
=============

The ``check_install`` command runs a series of checks to ensure that gTranslate is properly installed and configured on the user's system: all third-party software is on the system path and
all model files are present in the directory pointed to by ``GTRANSLATE_MODEL_PATH``.

Arguments
---------

.. argparse::
   :module: gtranslate.cli
   :func: get_main_parser
   :prog: gtranslate
   :path: check_install
   :nodefaultconst:

Example
-------

Input
^^^^^

.. code-block:: bash

    gtranslate check_install

Output
^^^^^^

.. code-block:: text

    [2026-05-04 17:07:18] INFO: gTranslate v{{ gtranslate_version }}
    [2026-05-04 17:07:18] INFO: gtranslate check_install
    [2026-05-04 17:07:18] INFO: Running install verification
    [2026-05-04 17:07:18] INFO: Checking that all third-party software are on the system path:
    [2026-05-04 17:07:18] INFO:          |-- prodigal         OK
    [2026-05-04 17:07:18] INFO: Checking that all model files are present:
    [2026-05-04 17:07:18] INFO: Check folder GTRANSLATE_MODEL_PATH (/path/to/gtranslate/models/): OK
    [2026-05-04 17:07:18] INFO: Check file ada_multi_class.pkl.gz (/path/to/gtranslate/models/ada_multi_class.pkl.gz): OK
    [2026-05-04 17:07:18] INFO: Check file dt_multi_class.pkl.gz (/path/to/gtranslate/models/dt_multi_class.pkl.gz): OK
    [2026-05-04 17:07:18] INFO: Check file knn_multi_class.pkl.gz (/path/to/gtranslate/models/knn_multi_class.pkl.gz): OK
    [2026-05-04 17:07:18] INFO: Check file xgb_multi_class.pkl.gz (/path/to/gtranslate/models/xgb_multi_class.pkl.gz): OK
    [2026-05-04 17:07:18] INFO: Check file mlp_multi_class.pkl.gz (/path/to/gtranslate/models/mlp_multi_class.pkl.gz): OK
    [2026-05-04 17:07:18] INFO: Check file label_encoder.pkl.gz (/path/to/gtranslate/models/label_encoder.pkl.gz): OK
    [2026-05-04 17:07:18] INFO: Done.