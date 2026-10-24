.. _installing:

Installing gTranslate
=====================

gTranslate is available through multiple sources, you only need to choose one.
If you are unsure which one to choose, Bioconda is generally the easiest.


Sources
-------


.. toctree::
   :maxdepth: 1

   bioconda
   pip


Hardware requirements
---------------------

Processing time is dominated by gene calling with Prodigal
and is roughly 30 seconds per genome per CPU, e.g. 5,000 genomes took 26 minutes using 90 CPUs. Use ``--cpus`` to
process genomes in parallel.


Python libraries
----------------

gTranslate requires Python >= 3.12 and the following libraries, which are installed automatically:

.. list-table::
   :widths: 10 10 80
   :header-rows: 1

   * - Library
     - Version
     - Reference
   * - `NumPy <https://numpy.org/>`_
     - ~= 1.26
     - Harris, C.R., Millman, K.J., van der Walt, S.J. et al. Array programming with NumPy. Nature 585, 357–362 (2020). DOI: `10.1038/s41586-020-2649-2 <https://doi.org/10.1038/s41586-020-2649-2>`_
   * - `Pandas <https://pandas.pydata.org/>`_
     - ~= 2.2
     - McKinney W. 2010. Data Structures for Statistical Computing in Python. Proceedings of the 9th Python in Science Conference, 51-56.
   * - `scikit-learn <https://scikit-learn.org/stable/>`_
     - == 1.6.1
     - Pedregosa F, et al. 2011. Scikit-learn: Machine Learning in Python. Journal of Machine Learning Research, 12, 2825-2830.
   * - `xgboost <https://xgboost.readthedocs.io/en/stable/>`_
     - == 2.1.4
     - Chen T, et al. 2016. XGBoost: A Scalable Tree Boosting System. Proceedings of the 22nd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining, 785–794. DOI: `10.1145/2939672.2939785 <https://doi.org/10.1145/2939672.2939785>`_
   * - `joblib <https://joblib.readthedocs.io/en/latest/>`_
     - ~= 1.3
     - Joblib: https://joblib.readthedocs.io/en/latest/
   * - `tqdm <https://github.com/tqdm/tqdm>`_
     - ~= 4.67
     - DOI: `10.5281/zenodo.595120 <https://doi.org/10.5281/zenodo.595120>`_
   * - `plotly <https://plotly.com/>`_
     - ~= 5.15
     - Plotly Technologies Inc. 2015. Collaborative data science. Montréal, QC: Plotly Technologies Inc. https://plot.ly

scikit-learn and XGBoost are pinned to exact versions as the gTranslate models were saved with these versions
and may not load correctly with other versions.

The ``lightgbm``, ``mlxtend``, ``requests`` and ``scipy`` packages are also currently installed as dependencies,
but are not used by gTranslate.

Please cite the libraries above if you use gTranslate in your work.


.. _installing#third-party-software:

Third-party software
--------------------

gTranslate makes use of the following 3rd party dependencies and assumes they are on your system path:

.. tip::
   The :ref:`commands/check_install` command will verify that all of the programs are on the path.

.. list-table::
   :widths: 10 10 80
   :header-rows: 1

   * - Software
     - Version
     - Reference
   * - `Prodigal <http://compbio.ornl.gov/prodigal/>`_
     - >= 2.6.2
     - Hyatt D, et al. 2010. `Prodigal: prokaryotic gene recognition and translation initiation site identification <https://www.ncbi.nlm.nih.gov/pubmed/20211023>`_. *BMC Bioinformatics*, 11:119. doi: 10.1186/1471-2105-11-119.


Please cite these tools if you use gTranslate in your work.


.. _installing#gtranslate-models:

gTranslate Models
-----------------

Due to file size limits, the gTranslate models are provided as a separate download. gTranslate 0.0.4 and later
use the same models; models downloaded for earlier versions are not compatible.

Downloading the Models
^^^^^^^^^^^^^^^^^^^^^^

The models are available from three mirrors. For optimal speeds, choose the mirror geographically closest to you:

* **Europe (GTDB Main):** https://data.gtdb.aau.ecogenomic.org/tools/gtranslate/
* **Australia (UQ Nectar):** https://data.gtdb.ecogenomic.org/tools/gtranslate/
* **Australia (UQ ACE):** https://data.ace.uq.edu.au/public/gtdb/tools/gtranslate/

Download the archive from **one** of the mirrors and extract it:

.. code-block:: bash

    # Europe (GTDB Main)
    wget https://data.gtdb.aau.ecogenomic.org/tools/gtranslate/gtranslate_models.tar.gz

    # OR Australia (UQ Nectar)
    # wget https://data.gtdb.ecogenomic.org/tools/gtranslate/gtranslate_models.tar.gz

    # OR Australia (UQ ACE)
    # wget https://data.ace.uq.edu.au/public/gtdb/tools/gtranslate/gtranslate_models.tar.gz

    # Extract the downloaded archive
    tar xvzf gtranslate_models.tar.gz

Configuring the Model Path
^^^^^^^^^^^^^^^^^^^^^^^^^^

gTranslate requires the ``GTRANSLATE_MODEL_PATH`` environment variable, or the ``--custom_model_path`` argument
of ``detect_table``, to point to the directory containing the extracted model files:

* ``ada_multi_class.pkl.gz``
* ``dt_multi_class.pkl.gz``
* ``knn_multi_class.pkl.gz``
* ``mlp_multi_class.pkl.gz``
* ``xgb_multi_class.pkl.gz``
* ``label_encoder.pkl.gz``

Run :ref:`commands/check_install` to verify that the model files are found.

For detailed instructions on setting this variable, please refer to the documentation for your specific setup:

* :ref:`pip installation <installing/pip>`
* :ref:`Bioconda installation <installing/bioconda>`
