"""
Configuration file for Machine Learning Classifiers and their best parameters.
This structure allows you to easily add new models, custom metadata (like 'name'), 
and update parameters.
"""
from sklearn.tree import DecisionTreeClassifier

classifier_configs = [
    {
        "name": "KNeighbors",
        "short_name": "knn",
        "custom_label": "Baseline KNN Model",
        "best_params": {
            'model__metric': 'manhattan',
             'model__n_jobs': 1,
             'model__n_neighbors': 2,
             'model__weights': 'uniform'}
    },
    # {
    #     "name": "RandomForest",
    #     "short_name": "rf",
    #     "custom_label": "Production RF - V1",
    #     "best_params":
    #         {'model__class_weight': 'balanced',
    #          'model__min_samples_leaf': 2,
    #          'model__min_samples_split': 6,
    #          'model__n_estimators': 50,
    #          'model__n_jobs': 1}
    # },
    {
        "name": "AdaBoost",
        "short_name": "ada",
        "custom_label": "AdaBoost Ensemble",
        "best_params": {
            'model__estimator': DecisionTreeClassifier(max_depth=3),
            'model__learning_rate': 1,
            'model__n_estimators': 50,
        }
    },
    {

        "name": "XGBoost",
        "short_name": "xgb",
        "custom_label": "XGBoost Gradient Boosting",
        "best_params": {
            'model__eval_metric': 'mlogloss',
            'model__learning_rate': 0.1,
            'model__n_estimators': 50,
            'model__n_jobs': 1,
            'model__subsample': 0.5
        }
    },
    {
        "name": "MLP",
        "short_name": "mlp",
        "custom_label": "Multi-layer Perceptron",
        "best_params": {
            'model__activation': 'tanh',
            'model__alpha': 0.0001,
            'model__hidden_layer_sizes': (20, 10),
            'model__learning_rate': 'constant',
            'model__max_iter': 2000,
            'model__solver': 'adam'
        }
    },
    {
        "name": "DecisionTree",
        "short_name": "dt",
        "custom_label": "Simple Interpretable Tree",
        "best_params": {
            'model__class_weight': 'balanced',
            'model__criterion': 'gini',
            'model__min_samples_leaf': 1,
            'model__min_samples_split': 2,
        }
    }
]