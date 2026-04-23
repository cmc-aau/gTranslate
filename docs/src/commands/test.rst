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
    
    [2020-04-13 09:50:58] INFO: GTDB-Tk v1.1.0
    [2020-04-13 09:50:58] INFO: gtdbtk test --out_dir /tmp/test --cpus 3
    [2020-04-13 09:50:58] INFO: Using GTDB-Tk reference data version r89: /release89
    [2020-04-13 09:50:58] INFO: Command: gtdbtk classify_wf --genome_dir /tmp/test/genomes --out_dir /tmp/test/output --cpus 3
    [2020-04-13 09:52:35] INFO: Test has successfully finished.
