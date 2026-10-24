.. _installing/pip:

pip
===


Step 1: Install third-party dependencies
----------------------------------------

Ensure that the software described in :ref:`installing#third-party-software` are on the system path.


Step 2: Install gTranslate via pip
----------------------------------

Note: It is strongly recommended to create a `virtual environment <https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/>`_
for each version of gTranslate installed.


Once the third-party dependencies have been installed, install gTranslate via `pip <https://pypi.org/project/gtdb-gtranslate/>`_:

.. code-block:: bash

    python -m pip install gtdb-gtranslate

Alternatively, if you have a previously installed gTranslate, you can upgrade the latest version:

.. code-block:: bash

    python -m pip install --upgrade gtdb-gtranslate


Step 3: Download and alias the gTranslate models
------------------------------------------------

gTranslate locates the models using the ``GTRANSLATE_MODEL_PATH`` environment variable, which must be set to the
directory containing the models. Alternatively, the ``--custom_model_path`` argument of ``detect_table`` can be used.


.. code-block:: bash

    export GTRANSLATE_MODEL_PATH=/path/to/models/

You can permanently save this variable as described `here <https://unix.stackexchange.com/questions/26047/how-to-correctly-add-a-path-to-path>`_.
