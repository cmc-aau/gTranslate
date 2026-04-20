import csv
import os
from typing import Dict, Any
from gtranslate.biolib_lite.common import make_sure_path_exists
from gtranslate.config.output import PATH_FEATURES_SUMMARY
from gtranslate.exceptions import GTranslateExit


class FeatureFile(object):
    """Records the features used by classifiers and the final predicted translation table."""

    def __init__(self, out_dir: str, prefix: str):
        # Name the output file
        self.path: str = os.path.join(out_dir, PATH_FEATURES_SUMMARY.format(prefix=prefix))
        self.rows: Dict[str, Dict[str, Any]] = dict()
        self.feature_columns = []
        self.none_value = 'N/A'

    def add_row(self, gid: str, predicted_table: int, features: Dict[str, Any]):
        """Add a genome's features and prediction to the file."""
        if gid in self.rows:
            raise GTranslateExit(f'Attempting to add duplicate row to feature file: {gid}')

        # Track unique feature columns dynamically
        for key in features.keys():
            if key not in self.feature_columns:
                self.feature_columns.append(key)

        # Build the row dictionary
        row_data = {
            'user_genome': gid,
            'predicted_tln_table': predicted_table
        }
        row_data.update(features)

        self.rows[gid] = row_data

    def write(self,training_features=False) -> None:
        """Writes the feature file using csv.DictWriter."""
        if not self.rows:
            return  # Nothing to write

        make_sure_path_exists(os.path.dirname(self.path))

        # Put genome and prediction first, followed by all dynamically found features
        if training_features:
            columns_names = ['Genome ID']+ self.feature_columns
            data_key_map = {'Genome ID': 'user_genome'}
        else:
            columns_names = ['user_genome', 'predicted_tln_table'] + self.feature_columns
            data_key_map = {}

        with open(self.path, 'w', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=columns_names, delimiter='\t')
            writer.writeheader()

            for gid, row_dict in sorted(self.rows.items()):
                clean_dict = {}
                for col in columns_names:
                    lookup_key = data_key_map.get(col, col)
                    val = row_dict.get(lookup_key, self.none_value)

                    # Clean up data types for export
                    if val is None:
                        val = self.none_value
                    elif isinstance(val, float):
                        val = round(val, 5)

                    clean_dict[col] = val

                writer.writerow(clean_dict)