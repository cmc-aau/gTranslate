#!/usr/bin/env python

###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__prog_name__ = 'generate_training_set.py'
__prog_desc__ = ('Generate the a new training set based on a gtdb taxonomy. Parse the metadata file from gtdb,'
                 'select the representative genomes and up to 99 additional genomes per species. For this set of genomes '
                 'The script will calculate the feature values - coding density, GC content, Trp and Gly ratios - and save the resulting table for training the classifiers.'
                 'The script will aslo create a file use to pick the ground truth Translation table.')

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2026'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.3'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'uqpchaum@uq.edu.au'
__status__ = 'Development'

import argparse
import os
import sys

import joblib
import numpy as np
import pandas as pd

# Sklearn Imports
from sklearn.ensemble import AdaBoostClassifier, RandomForestClassifier
from sklearn.metrics import balanced_accuracy_score
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.preprocessing import MinMaxScaler, LabelEncoder
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import Pipeline
from sklearn.utils import compute_class_weight
from lightgbm import LGBMClassifier

from mlxtend.feature_selection import SequentialFeatureSelector
from script_logger import CustomLogger

# Import classifier configurations
from classifier_parameters import classifier_configs

# prevents warnings from being displayed.
import warnings

warnings.simplefilter("ignore")


class TrainingSetGenerator(object):
    """
    A class to merge data, perform feature engineering, and train multiple
    classifiers with pipeline-based feature scaling and model persistence.
    """

    def __init__(self, feature_file, output_dir, threads_count, seed=None):

        self.feature_file = feature_file
        self.output_dir = output_dir
        self.thread_count = threads_count
        self.seed = seed

        # Initialize the logger
        logger_instance = CustomLogger(output_dir, __prog_name__)
        self.logger = logger_instance.get_logger()

        # make sure the output directory exists
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        # Map for instantiating raw models
        self.model_mapping = {
            "KNeighbors": KNeighborsClassifier(n_jobs=self.thread_count),
            "RandomForest": RandomForestClassifier(random_state=self.seed, n_jobs=self.thread_count),
            "AdaBoost": AdaBoostClassifier(random_state=self.seed),
            "LGBM": LGBMClassifier(random_state=self.seed, n_jobs=self.thread_count, verbose=-1),
            "DecisionTree": DecisionTreeClassifier(random_state=self.seed)
        }

    def load_and_preprocess_data(self):
        """Loads the 4 input TSVs, merges them, and engineers features."""
        self.logger.info("Loading and merging data files...")

        try:
            df_f = pd.read_csv(self.feature_file, sep='\t')
        except FileNotFoundError as e:
            self.logger.error(f"Error loading files: {e}")
            raise

        df_f['Genome'] = df_f['Genome'].astype(str).str.strip()


        # Filter Target
        valid_targets = ['11', '4', '25']
        df = df_f[df_f['Ground_truth'].astype(str).isin(valid_targets)].copy()

        self.logger.info(f"Merged dataset shape after target filtering: {df.shape}")

        # Feature Engineering
        df['Coding_density_4'] = df['Coding_density_4']
        df['Coding_density_11'] = df['Coding_density_11']
        df['Density_Diff'] = df['Coding_density_4'] - df['Coding_density_11']

        df['Trp_ratio'] = np.log((df['tt4_uga_count'] + 1) / (df['tt4_ugg_count'] + 1)).clip(-6.0, 5.0)
        df['Trp_magnitude'] = np.log(df['tt4_uga_count'] + df['tt4_ugg_count'] + 1)
        df['Gly_ratio'] = np.log((df['tt4_uga_count'] + 1) / (df['tt4_std_gly_count'] + 1)).clip(-10.0, 0.0)

        #df['feature_uga_density'] = df['tt4_uga_count'] / (df['Coding_density_4'] + 0.001)
        df['UGG_density'] = df['tt4_ugg_count'] / df['tt4_std_gly_count']
        #df['feature_uga_vs_gly_pool'] = df['tt4_uga_count'] / (df['tt4_std_gly_count'] + df['tt4_uga_count'] + 1)

        self.feature_cols = [
            'Coding_density_4', 'Coding_density_11', 'Density_Diff', 'GC',
            'Trp_ratio', 'Trp_magnitude', 'Gly_ratio',
            'UGG_density'
        ]

        # Handle NaNs from mathematical operations
        if df[self.feature_cols].isnull().any().any():
            self.logger.warning("NaN values generated during feature engineering. Filling with 0.")
            df[self.feature_cols] = df[self.feature_cols].fillna(0)

        return df

    def run(self, split_data=False):
        """Execute the full training pipeline across multiple classifiers."""
        self.logger.info(f'{__prog_name__} {" ".join(sys.argv[1:])}')

        seed_to_use = self.seed if self.seed is not None else np.random.randint(0, 1000)
        self.logger.info(f"We use seed: {seed_to_use}")

        # 1. Load & Preprocess
        df = self.load_and_preprocess_data()





        # 2. Prepare Data
        X = df[self.feature_cols]
        labels_raw = df['Ground_truth'].astype(str)

        # lets show the first few rows of the dataframe to check if everything is correct
        self.logger.info("First few rows of the preprocessed dataframe:")
        self.logger.info(X.head().to_string())

        # Encode multi-class labels safely
        le = LabelEncoder()
        y = le.fit_transform(labels_raw)

        # Save LabelEncoder to decode predictions later
        joblib.dump(le, os.path.join(self.output_dir, 'label_encoder.pkl'))

        # count class distribution
        class_counts = np.bincount(y)
        self.logger.info(f"Class counts: {class_counts}")

        class_weights = compute_class_weight(class_weight='balanced', classes=np.unique(y), y=y)
        weights_dict = dict(zip(le.classes_, class_weights))
        self.logger.info(f"Class weights mapping: {weights_dict}")

        if split_data:
            self.logger.info("Splitting the data into training and validation sets")
            X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, stratify=y,
                                                              random_state=seed_to_use)
            self.logger.info(f"Training data shape: {X_train.shape}")
        else:
            self.logger.info("Using the full training data")
            X_train, X_val, y_train, y_val = X, X, y, y

        # 3. Iterate Over Configurations
        for config in classifier_configs:
            model_name = config["name"]
            short_name = config["short_name"]
            best_params = config["best_params"]

            self.logger.info(f"--- Training {model_name} ---")

            base_model = self.model_mapping.get(model_name)
            if base_model is None:
                self.logger.warning(f"Model {model_name} not found in mapping. Skipping.")
                continue

            # Build Pipeline (MinMaxScaler for distance-based models)
            if model_name in ['KNeighbors']:
                pipe = Pipeline([
                    ('scaler', MinMaxScaler(feature_range=(0, 1))),
                    ('model', base_model)
                ])
            else:
                pipe = Pipeline([
                    ('model', base_model)
                ])

            # Apply best parameters dynamically
            pipe.set_params(**best_params)

            X_train_current = X_train.copy()
            X_val_current = X_val.copy()

            # Fit the pipeline
            pipe.fit(X_train_current, y_train)

            # Save the fitted pipeline
            model_out_path = os.path.join(self.output_dir, f'{short_name}_multi_class.pkl')
            try:
                joblib.dump(pipe, model_out_path)
                self.logger.info(f"Saved {model_name} pipeline to {model_out_path}")
            except Exception as e:
                self.logger.error(f"Error saving {model_name}: {e}")
                continue

            # Validation
            predictions = pipe.predict(X_val_current)
            balanced_accuracy = balanced_accuracy_score(y_val, predictions)
            self.logger.info(f"Balanced accuracy for {model_name}: {balanced_accuracy:.4f}\n")


if __name__ == "__main__":
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Input files matching your execution loop script
    parser.add_argument('--metadata_file_bac', required=True, help='GTDB metadata file for bacteria (tsv).')
    parser.add_argument('--metadata_file_arc', required=True, help='GTDB metadata file for archaea (tsv).')
    parser.add_argument('--genome_path_file', required=True, help='File containing paths to genome assemblies (tsv).')
    parser.add_argument('--output_dir', '-od', required=True, help='Directory to save the features and ground_truth tables.')
    args = parser.parse_args()

    try:
        classifier = TrainingSetGenerator(
            feature_file=args.features,
            output_dir=args.od,
            threads_count=args.threads,
            seed=args.seed
        )
        classifier.run(args.split_data)
    except Exception as e:
        print(f"ERROR: {str(e)}", file=sys.stderr)
        sys.exit(1)