.. _files/model_multi_class.pkl.gz:

model_multi_class.pkl.gz
========================

The five multi-class classifiers trained to predict the translation table from the features produced by gTranslate,
saved as compressed scikit-learn pipelines:

* ``ada_multi_class.pkl.gz``: AdaBoost
* ``dt_multi_class.pkl.gz``: decision tree
* ``knn_multi_class.pkl.gz``: k-nearest neighbours
* ``mlp_multi_class.pkl.gz``: multi-layer perceptron
* ``xgb_multi_class.pkl.gz``: XGBoost

Together with :ref:`label_encoder.pkl.gz <files/label_encoder.pkl.gz>`, these are the files ``detect_table`` loads
from ``--custom_model_path`` or ``GTRANSLATE_MODEL_PATH``.


Produced by
-----------

* :ref:`commands/fit_models`
