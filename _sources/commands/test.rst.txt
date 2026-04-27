.. _commands/test:

test
====

The ``test`` command is used to run three genomes with 3 different translation tables with the ``detect_table`` command.

If your installation is unable to run the ``test`` command with an exit code of ``0``, then
there is an issue with your installation.

Arguments
---------

.. argparse::
   :module: gtranslate.cli
   :func: get_main_parser
   :prog: gtranslate
   :path: test
   :nodefaultconst:


Files output
------------



Example
-------

Input
^^^^^

.. code-block:: bash

    gtranslate test --out_dir /tmp/test --cpus 3

Output
^^^^^^

.. code-block:: text
    
    [2026-04-08 17:05:24] INFO: gTranslate v0.0.2
    [2026-04-08 17:05:24] INFO: gtranslate test --out_dir run_test_command
    [2026-04-08 17:05:24] INFO: Command: gtranslate detect_table --genome_dir run_test_command/genomes --out_dir run_test_command/output --cpus 1
    [2026-04-08 17:06:22] INFO: Test has successfully finished.
