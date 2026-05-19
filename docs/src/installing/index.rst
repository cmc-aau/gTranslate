.. _installing:

Installing gTranslate
=====================

IN PROGRESS

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


Python libraries
----------------

gTranslate is designed for Python >=3.12 and requires the following libraries, which will be automatically installed:

.. list-table::
   :widths: 10 10 80
   :header-rows: 1

   * - Library
     - Version
     - Reference
   * - `NumPy <https://numpy.org/>`_
     - >= 1.26.0
     - Harris, C.R., Millman, K.J., van der Walt, S.J. et al. Array programming with NumPy. Nature 585, 357–362 (2020). DOI: `0.1038/s41586-020-2649-2 <https://doi.org/10.1038/s41586-020-2649-2>`_
   * - `tqdm <https://github.com/tqdm/tqdm>`_
     - >= 4.67.0
     - DOI: `10.5281/zenodo.595120 <https://doi.org/10.5281/zenodo.595120>`_
   * - `Pandas <https://pandas.pydata.org/>`_
     - >= 2.2.0
     - McKinney W. 2010. Data Structures for Statistical Computing in Python. Proceedings of the 9th Python in Science Conference, 51-56.
   * - `scikit-learn <https://scikit-learn.org/stable/>`_
     - >= 1.6.1
     - Pedregosa F, et al. 2011. Scikit-learn: Machine Learning in Python. Journal of Machine Learning Research, 12, 2825-2830.
   * - `joblib <https://joblib.readthedocs.io/en/latest/>`_
     - >= 1.3.2
     - Joblib: https://joblib.readthedocs.io/en/latest/
   * - `scipy <https://www.scipy.org/>`_
     - >= 1.12.0
     - Virtanen P, et al. 2020. SciPy 1.0: fundamental algorithms for scientific computing in Python. Nature Methods, 17, 261–272 (2020). DOI: `10.1038/s41592-019-0686-2 <https://doi.org/10.1038/s41592-019-0686-2>`_
   * -`mlxtend <http://rasbt.github.io/mlxtend/>`_
     - >= 0.22.0
     - Raschka S. 2018. MLxtend: Providing machine learning and data science utilities and extensions to Python's scientific computing stack. Journal of Open Source Software, 3(24), 638, https://doi.org/10.21105/joss.00638
   * - `plotly <https://plotly.com/>`_
     - >= 5.15.0
     - Plotly Technologies Inc. 2015. Collaborative data science. Montréal, QC: Plotly Technologies Inc. https://plot.ly
   * - `xgboost <https://xgboost.readthedocs.io/en/stable/>`_
     - >= 2.0.0
     - Chen T, et al. 2016. XGBoost: A Scalable Tree Boosting System. Proceedings of the 22nd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining, 785–794. DOI: `10.1145/2939672.2939785 <https://doi.org/10.1145/2939672.2939785>`_
   * - `lightgbm <https://lightgbm.readthedocs.io/en/latest/>`_
     - >= 3.3.5
     - Ke G, et al. 2017. LightGBM: A Highly Efficient Gradient Boosting Decision Tree. Advances in Neural Information Processing Systems, 30, 3146–3154.
   * - `requests <https://docs.python-requests.org/en/latest/>`_
     - >= 2.31.0
     - Reitz K. and Kenneth Reitz. 2023. Requests: HTTP for Humans. https://docs.python-requests.org/en/latest/

Please cite these libraries if you use gTranslate in your work.


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

gTranslate models
-----------------

Due to size limit, most of the models are available through a separate download.

.. code-block:: bash

    wget https://data.ace.uq.edu.au/public/to/finalise/gtranslate_models/gtranslate_models.tar.gz
    wget https://data.gtdb.ecogenomic.org/public/to/finalise/gtranslate_models/gtranslate_models.tar.gz ( mirror for Australia)
    tar xvzf gtranslate_models.tar.gz


**Alias the gTranslate models:**

gTranslate requires an environment variable named ``GTRANSLATE_MODEL_PATH`` to be set to the directory
containing the unarchived reference data. This is documented under:

- :ref:`pip installation <installing/pip>`
- :ref:`Bioconda installation <installing/bioconda>`


