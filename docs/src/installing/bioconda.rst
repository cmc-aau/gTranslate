.. _installing/bioconda:

Bioconda
========

Step 1: Install conda (if not already done)
-------------------------------------------

We strongly recommend using `Mamba <https://mamba.readthedocs.io/en/latest/installation.html>`_ (much faster!) over `miniconda <https://docs.conda.io/en/latest/miniconda.html>`_/`conda <https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html>`_, but all will work.


Step 2: Create the gTranslate environment
-----------------------------------------

It is strongly recommended to create a new conda environment for each version of gTranslate released.


gTranslate requires third-party packages from the ``conda-forge`` and ``bioconda`` channels, make sure to
specify those channels in that order!

.. include:: ../includes/install_block.rst
