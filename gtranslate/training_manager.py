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

import logging
import os
import gzip
import shutil
from collections import defaultdict

import joblib
import numpy as np
import pandas as pd

from gtranslate.classifiers.classifiers_parameters.classifier_parameters import classifier_configs
from gtranslate.config.common import CONFIG
from gtranslate.config.output import *
from gtranslate.external.prodigal import Prodigal
from gtranslate.files.featurefile import FeatureFile
from gtranslate.tools import tqdm_log, symlink_f, remove_intermediate_files

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


class TrainingManager(object):
    """Determine the ground truth for genomes based on their taxonomic classification."""

    def __init__(self, cpus: int = 1 , seed: int | None = None) -> None:
        """Initialize."""

        self.cpus = cpus
        self.seed = seed if seed is not None else np.random.randint(0, 1000)

        # Ground truth from GTDB classifications
        self.GTDB_TT25 = set(['o__Absconditabacterales', 'o__BD1-5'])
        self.GTDB_TT4 = set(['o__Mycoplasmatales', 's__Zinderia insecticola'])

        # Eggerthellacea genera using table 4; will need to be updated to names in Parks et al., 2026
        # once these appear in GTDB
        self.GTDB_TT4.update(set(['g__CAVGFB01', 'g__JAUNQF01']))

        # Minisyncoccia family identified in gTranslate manuscript that uses table 4. The majority, but not
        # all genomes in g__GCA-2747955 were also identified as using table 4. Currently, this is handled by
        # explicitly indicating the species in this genus identified as using table 4.
        self.GTDB_TT4.update(set(['f__JAKLIH01', 's__GCA-2747955 sp027024305', 's__GCA-2747955 sp027039745',
                                  's__GCA-2747955 sp947311625']))

        # Must include the Fastidiosibacteraceae XS4 species cluster once (if) this genome appears in GTDB:
        #  - https://www.ncbi.nlm.nih.gov/nuccore/AP038919.1
        #  - https://pmc.ncbi.nlm.nih.gov/articles/PMC12213064

        # Ground truth from NCBI classifications
        self.NCBI_TT4 = set(['s__Candidatus Hodgkinia cicadicola', 's__Candidatus Nasuia deltocephalincola', 's__Candidatus Stammera capleta'])
        self.NCBI_TT4.update(set(['s__Hodgkinia cicadicola', 's__Nasuia deltocephalincola', 's__Stammera capleta']))

        # species appears to have incorrect spelling at NCBI
        self.NCBI_TT4.update(set(['s__Candidatus Nasuia deltocephalinicola', 's__Nasuia deltocephalinicola']))

        # These species clusters have an unclear ground truth, see https://doi.org/10.1093/gbe/evad164
        self.GTDB_UNRESOLVED = set(['s__Providencia_A siddallii', 's__Providencia_A siddallii_A'])

        self.prefix = "gtranslate_training"

        self.logger = logging.getLogger('timestamp')

        self.model_mapping = {
            "KNeighbors": KNeighborsClassifier(n_jobs=self.cpus),
            "MLP": MLPClassifier(random_state=self.seed),
            "AdaBoost": AdaBoostClassifier(random_state=self.seed),
            "XGBoost": XGBClassifier(random_state=self.seed, n_jobs=self.cpus),
            "DecisionTree": DecisionTreeClassifier(random_state=self.seed)
        }

        self.feature_cols = [
            'Coding_density_4', 'Coding_density_11', 'Density_Diff', 'GC',
            'Trp_ratio', 'Trp_magnitude', 'Gly_ratio', 'UGG_density'
        ]

        self.genome_file_suffix = GENOME_FILE_SUFFIX
        self.protein_file_suffix = PROTEIN_FILE_SUFFIX
        self.nt_gene_file_suffix = NT_GENE_FILE_SUFFIX
        self.gff_file_suffix = GFF_FILE_SUFFIX
        self.checksum_suffix = CHECKSUM_SUFFIX

    def parse_manual_ground_truth_file(self, manual_gt_file: str) -> dict:
        """Parse manual ground truth file."""

        manual_ground_truth = {}
        open_file = gzip.open if manual_gt_file.endswith('.gz') else open
        with open_file(manual_gt_file, 'rt') as f:
            header = f.readline().strip().split('\t')
            gid_idx = header.index("Genome ID")
            gt_idx = header.index("Translation table")

            for line in f:
                tokens = line.strip().split('\t')
                manual_ground_truth[tokens[gid_idx]] = tokens[gt_idx]

        return manual_ground_truth

    def select_ground_truth(self, taxonomy_file: str, out_file: str, manual_gt_file: str | None = None) -> None:
        """Determine the ground truth for genomes based on their taxonomic classification."""

        # read files with manually specific ground truth
        manual_ground_truth = {}
        if manual_gt_file:
            self.logger.info('Parsing manual ground truth file:')
            manual_ground_truth = self.parse_manual_ground_truth_file(manual_gt_file)
            self.logger.info(f' - identified manual ground truth for {len(manual_ground_truth):,} genomes')

        # determine ground truth for genomes based on their taxonomic classification
        self.logger.info('Determining ground truth for genomes:')
        total_genomes = 0
        gt_table_count = defaultdict(int)
        num_by_manual_gt = 0

        open_file = gzip.open if taxonomy_file.endswith('.gz') else open
        with open_file(taxonomy_file, 'rt') as f:
            header = f.readline().strip().split('\t')

            if 'Genome ID' in header:
                gid_idx = header.index("Genome ID")
            else:
                self.logger.error("Taxonomy file must have a 'Genome ID' column.")

            fout = open(out_file, 'w')
            fout.write('Genome ID\tGround truth table')

            taxonomy_idx = None
            if 'Taxonomy' in header:
                taxonomy_idx = header.index("Taxonomy")
                fout.write('\tTaxonomy')

            gtdb_taxonomy_idx = None
            if 'GTDB taxonomy' in header:
                gtdb_taxonomy_idx = header.index('GTDB taxonomy')
                fout.write('\tGTDB taxonomy')

            ncbi_taxonomy_idx = None
            if 'NCBI taxonomy' in header:
                ncbi_taxonomy_idx = header.index('NCBI taxonomy')
                fout.write('\tNCBI taxonomy')

            if taxonomy_idx is None and gtdb_taxonomy_idx is None and ncbi_taxonomy_idx is None:
                self.logger.error("Taxonomy file must have a 'Taxonomy' column or a 'GTDB taxonomy' and 'NCBI taxonomy' columns.")

            fout.write('\n')

            for line in f:
                tokens = line.strip().split('\t')

                total_genomes += 1

                gid = tokens[gid_idx]

                if gid in manual_ground_truth:
                    ground_truth_tt = manual_ground_truth[gid]
                    num_by_manual_gt += 1
                else:
                    # determine ground truth translation table based on GTDB
                    # or NCBI taxonomic classification of genome
                    taxa = set()
                    if taxonomy_idx:
                        taxa.update(set(tokens[taxonomy_idx].split(';')))

                    gtdb_taxa = set()
                    if gtdb_taxonomy_idx:
                        gtdb_taxa.update(set(tokens[gtdb_taxonomy_idx].split(';')))

                    ncbi_taxa = set()
                    if ncbi_taxonomy_idx:
                        ncbi_taxa.update(set(tokens[ncbi_taxonomy_idx].split(';')))

                    if taxa.intersection(self.GTDB_TT25) or gtdb_taxa.intersection(self.GTDB_TT25):
                        ground_truth_tt = '25'
                    elif taxa.intersection(self.GTDB_TT4) or gtdb_taxa.intersection(self.GTDB_TT4):
                        ground_truth_tt = '4'
                    elif taxa.intersection(self.NCBI_TT4) or ncbi_taxa.intersection(self.NCBI_TT4):
                        ground_truth_tt = '4'
                    elif taxa.intersection(self.GTDB_UNRESOLVED) or gtdb_taxa.intersection(self.GTDB_UNRESOLVED):
                        ground_truth_tt = 'UNRESOLVED'
                    else:
                        ground_truth_tt = '11'

                gt_table_count[ground_truth_tt] += 1

                # write out ground truth results
                fout.write(f'{gid}\t{ground_truth_tt}')

                if taxonomy_idx:
                    fout.write(f'\t{tokens[taxonomy_idx]}')

                if gtdb_taxonomy_idx:
                    fout.write(f'\t{tokens[gtdb_taxonomy_idx]}')

                if ncbi_taxonomy_idx:
                    fout.write(f'\t{tokens[ncbi_taxonomy_idx]}')

                fout.write('\n')

        fout.close()

        # write out number of genomes assigned to each translation table
        self.logger.info(f' - determined ground truth for {total_genomes:,} genomes')
        if manual_gt_file:
            self.logger.info(f' - ground truth set manually for {num_by_manual_gt:,} genomes')

        for tran_table, genome_count in sorted(gt_table_count.items()):
            self.logger.info(f'Table {tran_table}: {genome_count:,} ({100*genome_count/total_genomes:.2f}%)')

    def build_features(self, genomes: str, out_dir: str,force:bool) -> None:
        """Generate feature vectors for training models."""

        reports = {}
        self.called_gene_dir = os.path.join(out_dir, DIR_PREDICT_GENES)
        self.failed_genomes = os.path.join(
            out_dir, PATH_FAILS.format(prefix=self.prefix))
        reports.setdefault('all', []).append(self.failed_genomes)

        prodigal = Prodigal(self.cpus,
                            self.failed_genomes,
                            self.called_gene_dir,
                            self.protein_file_suffix,
                            self.nt_gene_file_suffix,
                            self.gff_file_suffix,
                            force)
        self.logger.log(
            CONFIG.LOG_TASK, f'Running Prodigal {prodigal.version} to identify genes.')

        genome_dictionary = prodigal.run(genomes)

        gene_files = [(db_genome_id, genome_dictionary[db_genome_id]['aa_gene_path'])
                      for db_genome_id in genome_dictionary.keys()]
        # if gene_files is empty, we have no genomes to process
        # we still need to write all the genomes failing the identify step to the bac120_summary file
        if len(gene_files) == 0:
            self.logger.warning('All genomes failed the identify step.')
            # Exit gracefully
            return

        feature_file = FeatureFile(out_dir, self.prefix)

        for db_genome_id, info in tqdm_log(sorted(genome_dictionary.items()), unit='genome'):
            # Write the best translation table to disk for this genome.

            features = info.get("feature_vector", {})
            if features:
                feature_file.add_row(db_genome_id,0, features)

        feature_file.write(training_features=True)

        #we remove the outdir from the path to get a relative path
        relative_path_feature_file = os.path.relpath(feature_file.path, out_dir)
        # copy relative_path_feature_file to out_dir
        print(feature_file.path)
        print(os.path.join(out_dir, relative_path_feature_file))
        shutil.copyfile(feature_file.path, os.path.join(out_dir,os.path.basename(feature_file.path)))

        # remove predict folder
        self.logger.info('Cleaning up intermediate files.')
        remove_intermediate_files(out_dir)

    def fit_models(self, feature_file: str, tt_file: str, out_dir: str, split_data: bool) -> None:
        """Train models based on the feature vector tables."""

        seed_to_use = self.seed
        self.logger.info(f"We use seed: {seed_to_use}")

        df_gt_features = self._load_and_preprocess_data(feature_file, tt_file)

        # 2. Prepare Data
        X = df_gt_features[self.feature_cols]
        labels_raw = df_gt_features['Ground truth table'].astype(str)

        # lets show the first few rows of the dataframe to check if everything is correct
        self.logger.info("First few rows of the preprocessed dataframe:")
        self.logger.info(X.head().to_string())

        # Encode multi-class labels safely
        le = LabelEncoder()
        y = le.fit_transform(labels_raw)

        # Save LabelEncoder to decode predictions later
        joblib.dump(le, os.path.join(out_dir, 'label_encoder.pkl.gz'))

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
            if model_name in ['KNeighbors', 'MLP']:
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
            model_out_path = os.path.join(out_dir, f'{short_name}_multi_class.pkl.gz')
            try:
                # joblib detects the .gz extension and automatically compresses it
                joblib.dump(pipe, model_out_path)
                self.logger.info(f"Saved and compressed {model_name} pipeline to {model_out_path}")
                # gzip the model

            except Exception as e:
                self.logger.error(f"Error saving {model_name}: {e}")
                continue

            # Validation
            predictions = pipe.predict(X_val_current)
            balanced_accuracy = balanced_accuracy_score(y_val, predictions)
            self.logger.info(f"Balanced accuracy for {model_name}: {balanced_accuracy:.4f}\n")

    def _load_and_preprocess_data(self, feature_file: str, tt_file: str):
        """Load and preprocess the feature and ground truth data."""

        # Load feature data
        df_features = pd.read_csv(feature_file, sep='\t')
        df_gt = pd.read_csv(tt_file, sep='\t')

        # Merge on Genome ID
        df_merged = pd.merge(df_gt, df_features, on='Genome ID')


        # Check if all required feature columns already exist
        missing_cols = [col for col in self.feature_cols if col not in df_merged.columns]

        if not missing_cols:
            self.logger.info("All feature columns present in input. Skipping feature engineering.")
        else:
            self.logger.warning(f"Missing columns: {missing_cols}. Performing feature engineering...")

        # Handle NaNs before saving
        if df_merged[self.feature_cols].isnull().any().any():
            self.logger.warning("NaN values generated or found. Filling with 0.")
            df_merged[self.feature_cols] = df_merged[self.feature_cols].fillna(0)


        # =========================
        # Filter valid targets
        # =========================

        valid_targets = ['11', '4', '25']
        df = df_merged[df_merged['Ground truth table'].astype(str).isin(valid_targets)].copy()

        self.logger.info(f"Merged dataset shape after target filtering: {df.shape}")

        return df



