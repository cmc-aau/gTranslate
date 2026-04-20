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

__prog_name__ = 'train_classifier_multi.py'
__prog_desc__ = ('Generates multiple multi-class models and scalers using engineered '
                 'codon/glycine features. Includes data preprocessing, multi-model '
                 'training pipelines, and evaluation.')

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
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import MinMaxScaler, LabelEncoder
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import Pipeline
from sklearn.utils import compute_class_weight
from xgboost import XGBClassifier


from script_logger import CustomLogger

# Import classifier configurations
from classifier_parameters import classifier_configs

# prevents warnings from being displayed.
import warnings

warnings.simplefilter("ignore")


class ScalerClassifierMulti(object):
    """
    A class to merge data, perform feature engineering, and train multiple
    classifiers with pipeline-based feature scaling and model persistence.
    """

    def __init__(self, feature_file, output_dir, threads_count, seed=None,ledger=None):

        self.feature_file = feature_file
        self.output_dir = output_dir
        self.thread_count = threads_count
        self.seed = seed
        self.ledger = ledger

        # Initialize the logger
        logger_instance = CustomLogger(output_dir, __prog_name__)
        self.logger = logger_instance.get_logger()

        # make sure the output directory exists
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)


        # Map for instantiating raw models
        self.model_mapping = {
            "KNeighbors": KNeighborsClassifier(n_jobs=self.thread_count),
            "MLP": MLPClassifier(random_state=self.seed),
            "AdaBoost": AdaBoostClassifier(random_state=self.seed),
            "XGBoost": XGBClassifier(random_state=self.seed, n_jobs=self.thread_count, verbose=-1),
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

        # load the ledger if provided and merge it to overwrite the original ground truth for UNRESOLVED genomes
        if self.ledger is not None:
            try:
                # 1. Read and clean the ledger
                df_ledger = pd.read_csv(self.ledger, sep='\t')
                df_ledger['Genome'] = df_ledger['Genome'].astype(str).str.strip()
                df_ledger = df_ledger.rename(columns={'Translation_table': 'Ledger_TT'})

                # 2. Clean the main dataframe's Genome column to guarantee a perfect match
                df_f['Genome'] = df_f['Genome'].astype(str).str.strip()

                # 3. Merge the ledger data into the main dataframe
                df_f = pd.merge(df_f, df_ledger[['Genome', 'Ledger_TT']], on='Genome', how='left')

                # 4. New Mask: Target EVERY row where the ledger provided a value (not NaN)
                update_mask = df_f['Ledger_TT'].notna()

                # 5. Overwrite the Ground_truth unconditionally for those specific genomes
                df_f.loc[update_mask, 'Ground_truth'] = df_f.loc[update_mask, 'Ledger_TT']

                # 6. Clean up the temporary column
                df_f = df_f.drop(columns=['Ledger_TT'])

                self.logger.info(
                    f"Ledger merged: {df_ledger.shape[0]} entries. Force-updated {update_mask.sum()} genomes.")
            except Exception as e:
                self.logger.error(f"Error loading or merging ledger: {e}")
                raise

        # Define the target feature columns
        self.feature_cols = [
            'Coding_density_4', 'Coding_density_11', 'Density_Diff', 'GC',
            'Trp_ratio', 'Trp_magnitude', 'Gly_ratio', 'UGG_density'
        ]

        # =========================
        # Conditional Feature Engineering
        # =========================

        # Check if all required feature columns already exist
        missing_cols = [col for col in self.feature_cols if col not in df_f.columns]

        if not missing_cols:
            self.logger.info("All feature columns present in input. Skipping feature engineering.")
        else:
            self.logger.info(f"Missing columns: {missing_cols}. Performing feature engineering...")

            # Ensure base columns exist before calculating
            required_base_cols = ['Coding_density_4', 'Coding_density_11', 'tt4_uga_count', 'tt4_ugg_count',
                                  'tt4_std_gly_count']
            if not all(col in df_f.columns for col in required_base_cols):
                self.logger.error("Base columns for feature engineering are missing. Check input data.")
                raise ValueError(f"Missing base columns. Required: {required_base_cols}")

            df_f['Density_Diff'] = df_f['Coding_density_4'] - df_f['Coding_density_11']
            df_f['Trp_ratio'] = np.log((df_f['tt4_uga_count'] + 1) / (df_f['tt4_ugg_count'] + 1)).clip(-6.0,
                                                                                                       5.0)
            df_f['Trp_magnitude'] = np.log(df_f['tt4_uga_count'] + df_f['tt4_ugg_count'] + 1)
            df_f['Gly_ratio'] = np.log((df_f['tt4_uga_count'] + 1) / (df_f['tt4_std_gly_count'] + 1)).clip(
                -10.0, 0.0)
            df_f['UGG_density'] = (df_f['tt4_ugg_count'] / df_f['tt4_std_gly_count'])

        # Handle NaNs before saving
        if df_f[self.feature_cols].isnull().any().any():
            self.logger.warning("NaN values generated or found. Filling with 0.")
            df_f[self.feature_cols] = df_f[self.feature_cols].fillna(0)

        # =========================
        # Save FULL dataset (includes unresolved)
        # =========================

        cols_to_save = ['Genome', 'Ground_truth'] + self.feature_cols
        existing_cols = [c for c in cols_to_save if c in df_f.columns]

        output_path = os.path.join(self.output_dir, 'preprocessed_with_unresolved.tsv')
        df_f[existing_cols].to_csv(output_path, sep='\t', index=False)
        self.logger.info(f"Full dataset (including unresolved genomes) saved to {output_path}")

        # =========================
        # Filter valid targets
        # =========================

        valid_targets = ['11', '4', '25']
        df = df_f[df_f['Ground_truth'].astype(str).isin(valid_targets)].copy()

        self.logger.info(f"Merged dataset shape after target filtering: {df.shape}")

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
            if model_name in ['KNeighbors','MLP']:
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
    parser.add_argument('--features', required=True, help='Path to input_table.tsv')

    parser.add_argument('-o', '--output_dir', dest="od", required=True, help='Path to the output directory.')
    parser.add_argument('-s', '--seed', dest="seed", type=int, default=None, help='Seed for reproducibility.')
    parser.add_argument('-t', '--threads', dest="threads", type=int, default=1, help='Number of threads to use.')
    parser.add_argument('-l', '--ledger', dest="ledger", default=None, help='For genomes flagged as UNRESOLVED, this ledger give a manual ground truth to overwrite the original ground truth. '
                                                                            'It is a TSV file with 2 columns: Genome and Ground_truth. The Ground_truth column should have the same format as the original ground truth (e.g., 4, 11, 25,UNRESOLVED).')
    parser.add_argument('--split_data', dest="split_data", action='store_true',
                        help='Enable data splitting into training and validation sets.')

    args = parser.parse_args()

    try:
        classifier = ScalerClassifierMulti(
            feature_file=args.features,
            output_dir=args.od,
            threads_count=args.threads,
            seed=args.seed,
            ledger=args.ledger
        )
        classifier.run(args.split_data)
    except Exception as e:
        print(f"ERROR: {str(e)}", file=sys.stderr)
        sys.exit(1)