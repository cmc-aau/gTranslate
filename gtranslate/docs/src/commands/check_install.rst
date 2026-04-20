.. _commands/check_install:

check_install
=============

The `check_install` command is used to verify the integrity of the gTranslate.

If any inconsistencies are identified then the program will exit with code 1.

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

    [2025-03-20 11:51:30] INFO: gTranslate v0.0.1
    [2025-03-20 11:51:30] INFO: gtranslate check_install
    [2025-03-20 11:51:30] INFO: Running install verification
    [2025-03-20 11:51:30] INFO: Checking that all third-party software are on the system path:
    [2025-03-20 11:51:30] INFO:          |-- jellyfish        OK
    [2025-03-20 11:51:30] INFO:          |-- prodigal         OK
    [2025-03-20 11:51:30] INFO: Done.

