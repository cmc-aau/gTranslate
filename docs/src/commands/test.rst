.. _commands/test:

test
====

The ``test`` command is used to run three genomes with 3 different translation tables with the ``detect_table`` command.

If your installation is unable to run the ``test`` command with an exit code of ``0``, then
there is an issue with your installation.

.. note::
   The ``test`` command runs ``detect_table`` and therefore requires the gTranslate models.
   Download the models and set the ``GTRANSLATE_MODEL_PATH`` environment variable before running it
   (see :ref:`installing#gtranslate-models`).

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

Files are written to ``--out_dir``, or to a temporary directory which is deleted afterwards if ``--out_dir`` is
not specified:

* ``genomes/``: the three test genomes.
* ``output/``: the output of ``detect_table`` run on the test genomes, see :ref:`commands/detect_table`.
* ``test_execution.log``: the console output of ``detect_table``, which should be checked if the test fails.
* :ref:`gtranslate.log <files/gtranslate.log>`


Example
-------

Input
^^^^^

.. code-block:: bash

    gtranslate test --out_dir /tmp/test --cpus 3

Output
^^^^^^

.. code-block:: text
    
    [2026-04-08 17:05:24] INFO: gTranslate v{{ gtranslate_version }}
    [2026-04-08 17:05:24] INFO: gtranslate test --out_dir run_test_command
    [2026-04-08 17:05:24] INFO: Command: gtranslate detect_table --genome_dir run_test_command/genomes --out_dir run_test_command/output --cpus 1
    [2026-04-08 17:06:22] INFO: Test has successfully finished.
