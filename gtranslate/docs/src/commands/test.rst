.. _commands/test:

test
====

The ``test`` command is used to run three genomes ( one with translation table 11, one with translation table 4 and
one with translation table 25 ) through the detect_table workflow.

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

* :ref:`[prefix].log <files/gtranslate.warnings.log>`
* :ref:`[prefix].warnings.log <files/gtranslate.warnings.log>`
* :ref:`output/ <commands/detect_table>`
* :ref:`test_execution.log <files/test_execution.log>`

Example
-------

Input
^^^^^

.. code-block:: bash

    gtranslate test --out_dir test_command --cpus 1

Output
^^^^^^

.. code-block:: text

    [2025-03-17 15:45:31] INFO: gTranslate v0.0.1
    [2025-03-17 15:45:31] INFO: gtranslate test --out_dir test_command
    [2025-03-17 15:45:31] INFO: Command: gtranslate detect_table --genome_dir test_command/genomes --out_dir test_command/output --cpus 1
    [2025-03-17 15:46:01] INFO: Test has successfully finished.
