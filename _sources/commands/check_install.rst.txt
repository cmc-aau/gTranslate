.. _commands/check_install:

check_install
=============

The ``check_install`` command runs a series of checks to ensure that gTranslate is properly installed and configured on the user's system.

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

    [2026-05-04 17:07:18] INFO: gTranslate v0.0.2
    [2026-05-04 17:07:18] INFO: gtranslate check_install
    [2026-05-04 17:07:18] INFO: Running install verification
    [2026-05-04 17:07:18] INFO: Checking that all third-party software are on the system path:
    [2026-05-04 17:07:18] INFO:          |-- prodigal         OK
    [2026-05-04 17:07:18] INFO: Done.