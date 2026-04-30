.. _installing/pip:

pip
===


Step 1: Install third-party dependencies
----------------------------------------

Ensure that the software described in :ref:`installing#third-party-software` are on the system path.


Step 2: Install gTranslate via pip
-------------------------------

Note: It is strongly recommended to create a `virtual environment <https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/>`_
for each version of gTranslate installed.


Once the third-party dependencies have been installed, install gTranslate via `pip https://test.pypi.org/project/gtdb-gtranslate/`_:

.. code-block:: bash

    python -m pip install -i https://test.pypi.org/simple/ gtdb-gtranslate

Alternatively, if you have a previously installed gTranslate, you can upgrade the latest version:

.. code-block:: bash

    python -m pip install -i https://test.pypi.org/simple/ --upgrade

