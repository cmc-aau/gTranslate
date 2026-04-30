.. _commands/fit_models:

fit_models
==========

``fit_models`` takes the features built from ``build_features`` and trains the model to recognize the boundaries between Table 11, 4, and 25.
Once it’s done, you’ll have fresh model files ready for the ``detect_table`` command.


Arguments
---------

.. argparse::
   :module: gtranslate.cli
   :func: get_main_parser
   :prog: gtranslate
   :path: fit_models
   :nodefaultconst:


Files output
------------

* :ref:`label_encoder.pkl.gz <files/label_encoder.pkl.gz>`
* :ref:`[model]_multi_class.pkl.gz <files/model_multi_class.pkl.gz>`
* :ref:`[prefix].log <files/gtranslate.log>`
* :ref:`[prefix].warnings.log <files/gtranslate.warnings.log>`

Example
-------

Input
^^^^^

.. code-block:: bash

    gtranslate fit_models --feature_file features_test/gtranslate_training.feature_summary.tsv --tt_file 1000genomes/5000_genomes_r226.tsv --out_dir models

Output
^^^^^^

.. code-block:: text


[2026-04-10 19:22:07] INFO: gTranslate v0.0.2
[2026-04-10 19:22:07] INFO: gtranslate fit_models --feature_file features_test/gtranslate_training.feature_summary.tsv --tt_file 1000genomes/5000_genomes_r226.tsv --out_dir models
[2026-04-10 19:22:07] INFO: Training models based on training data.
[2026-04-10 19:22:07] INFO: We use seed: 427
[2026-04-10 19:22:07] INFO: All feature columns present in input. Skipping feature engineering.
[2026-04-10 19:22:07] INFO: Merged dataset shape after target filtering: (5000, 12)
[2026-04-10 19:22:07] INFO: First few rows of the preprocessed dataframe:
[2026-04-10 19:22:07] INFO:    Coding_density_4  Coding_density_11  Density_Diff        GC  Trp_ratio  Trp_magnitude  Gly_ratio  UGG_density
0          84.29031           85.94244      -1.65213  56.58347   -1.84939       10.08464   -3.39509      0.21315
1          86.46786           88.43902      -1.97116  57.48985   -1.46175       10.24775   -3.11130      0.19213
2          86.37359           88.21455      -1.84096  57.46309   -1.46933       10.24881   -3.12142      0.19164
3          83.01639           91.96291      -8.94652  71.41138   -1.36587        9.99150   -3.37393      0.13424
4          84.64863           86.23985      -1.59122  56.69494   -1.90344       10.05767   -3.46512      0.20978
[2026-04-10 19:22:07] INFO: Class counts: [4988    3    9]
[2026-04-10 19:22:07] INFO: Class weights mapping: {'11': 0.33413525795241916, '25': 555.5555555555555, '4': 185.1851851851852}
[2026-04-10 19:22:07] INFO: Using the full training data
[2026-04-10 19:22:07] INFO: --- Training KNeighbors ---
[2026-04-10 19:22:07] INFO: Saved and compressed KNeighbors pipeline to models/knn_multi_class.pkl.gz
[2026-04-10 19:22:07] INFO: Balanced accuracy for KNeighbors: 0.8889

[2026-04-10 19:22:07] INFO: --- Training AdaBoost ---
[2026-04-10 19:22:07] INFO: Saved and compressed AdaBoost pipeline to models/ada_multi_class.pkl.gz
[2026-04-10 19:22:07] INFO: Balanced accuracy for AdaBoost: 1.0000

[2026-04-10 19:22:07] INFO: --- Training XGBoost ---
[2026-04-10 19:22:08] INFO: Saved and compressed XGBoost pipeline to models/xgb_multi_class.pkl.gz
[2026-04-10 19:22:08] INFO: Balanced accuracy for XGBoost: 0.8889

[2026-04-10 19:22:08] INFO: --- Training MLP ---
[2026-04-10 19:22:08] INFO: Saved and compressed MLP pipeline to models/mlp_multi_class.pkl.gz
[2026-04-10 19:22:08] INFO: Balanced accuracy for MLP: 0.3333

[2026-04-10 19:22:08] INFO: --- Training DecisionTree ---
[2026-04-10 19:22:08] INFO: Saved and compressed DecisionTree pipeline to models/dt_multi_class.pkl.gz
[2026-04-10 19:22:08] INFO: Balanced accuracy for DecisionTree: 1.0000

[2026-04-10 19:22:08] INFO: Done.