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

import csv
import os
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional

from gtranslate.biolib_lite.common import make_sure_path_exists
from gtranslate.config.output import PATH_TLN_TABLE_SUMMARY
from gtranslate.exceptions import GTranslateExit


@dataclass(slots=True)
class TranslationSummaryFileRow:
    """A row contained within the ClassifySummaryFile object."""
    gid: str
    best_tln_table: Optional[int] = None
    coding_density_4: Optional[float] = None
    coding_density_11: Optional[float] = None
    gc_percent: Optional[float] = None
    n50: Optional[int] = None
    genome_size: Optional[int] = None
    contig_count: Optional[int] = None
    confidence: Optional[float] = None
    warnings: List[str] = field(default_factory=list)
    ensemble_preds: dict = field(default_factory=dict)

    def __post_init__(self):
        """Automatically enforce types upon creation."""
        # Cleanly cast floats-as-strings (e.g., "4.0") to true integers
        if self.best_tln_table is not None: self.best_tln_table = int(float(self.best_tln_table))
        if self.n50 is not None: self.n50 = int(float(self.n50))
        if self.genome_size is not None: self.genome_size = int(float(self.genome_size))
        if self.contig_count is not None: self.contig_count = int(float(self.contig_count))

        if self.coding_density_4 is not None: self.coding_density_4 = float(self.coding_density_4)
        if self.coding_density_11 is not None: self.coding_density_11 = float(self.coding_density_11)
        if self.gc_percent is not None: self.gc_percent = float(self.gc_percent)
        if self.confidence is not None: self.confidence = float(self.confidence)


class TranslationSummaryFile(object):
    """Records the translation table for one or more genomes."""
    fixed_columns_names = [
        'user_genome', 'best_tln_table', 'coding_density_4', 'coding_density_11',
        'gc_percent', 'n50', 'genome_size', 'contig_count','confidence','warnings']

    def __init__(self, out_dir: str, prefix: str):
        """Configure paths and initialise storage dictionary."""
        self.path: str = os.path.join(out_dir, PATH_TLN_TABLE_SUMMARY.format(prefix=prefix))
        self.rows: Dict[str, TranslationSummaryFileRow] = dict()
        self.none_value = 'N/A'  # Value to use for None values in the file.


    def add_row(self, row: TranslationSummaryFileRow):
        if row.gid in self.rows:
            raise GTranslateExit(f'Attempting to add duplicate row: {row.gid}')
        self.rows[row.gid] = row

    def get_row(self, gid: str) -> TranslationSummaryFileRow:
        if gid not in self.rows:
            raise GTranslateExit(f'Attempting to get non-existent row: {gid}')
        return self.rows[gid]

    def update_row(self, row: TranslationSummaryFileRow):
        if row.gid not in self.rows:
            raise GTranslateExit(f'Attempting to update non-existent row: {row.gid}')
        self.rows[row.gid] = row

    def has_row(self) -> bool:
        if self.rows.items():
            return True
        return False

    def write(self):
        """Writes the summary file using csv.DictWriter."""
        make_sure_path_exists(os.path.dirname(self.path))

        model_columns = []
        for row in self.rows.values():
            for model_name in row.ensemble_preds.keys():
                if model_name not in model_columns:
                    model_columns.append(model_name)

        # put the warnings column at the end
        columns_names = self.fixed_columns_names[:-1] + model_columns + [self.fixed_columns_names[-1]]


        with open(self.path, 'w', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=columns_names, delimiter='\t')
            writer.writeheader()

            for gid, row in sorted(self.rows.items()):
                # Convert dataclass to a dictionary
                row_dict = asdict(row)

                # Map 'gid' to 'user_genome' for the output file
                row_dict['user_genome'] = row_dict.pop('gid')

                # Format warnings list to string
                row_dict['warnings'] = ';'.join(row_dict['warnings']) if row_dict['warnings'] else self.none_value
                preds = row_dict.pop('ensemble_preds')
                for model in model_columns:
                    row_dict[model] = preds.get(model, self.none_value)

                # Replace None with 'N/A' and round floats
                for key, val in row_dict.items():
                    if val is None:
                        row_dict[key] = self.none_value
                    elif isinstance(val, float):
                        row_dict[key] = round(val, 5)

                writer.writerow(row_dict)

    def read(self):
        """Reads the summary file and dynamically repacks model columns into a dictionary."""
        if not os.path.isfile(self.path):
            raise GTranslateExit(f'Error, classify summary file not found: {self.path}')

        with open(self.path, newline='') as fh:
            reader = csv.DictReader(fh, delimiter='\t')

            # Validate that the core standard columns exist
            for col in self.fixed_columns_names:
                if col not in reader.fieldnames:
                    raise GTranslateExit(f'Missing standard column in file: {col}')

            # 1. Identify which columns are dynamic model columns (anything not standard)
            model_columns = [col for col in reader.fieldnames if col not in self.fixed_columns_names]

            for row_data in reader:
                clean_data = {k: (None if v == self.none_value or v == '' else v) for k, v in row_data.items()}

                warnings_str = clean_data.get('warnings')
                warnings_list = warnings_str.split(';') if warnings_str else []

                # 2. Rebuild the ensemble_preds dictionary from the trailing columns
                preds_dict = {}
                for model in model_columns:
                    val = clean_data.get(model)
                    if val is not None:
                        # Safely cast '4.0' or '4' strings back to clean integers
                        preds_dict[model] = int(float(val))

                row = TranslationSummaryFileRow(
                    gid=clean_data['user_genome'],
                    best_tln_table=clean_data['best_tln_table'],
                    coding_density_4=clean_data['coding_density_4'],
                    coding_density_11=clean_data['coding_density_11'],
                    gc_percent=clean_data['gc_percent'],
                    n50=clean_data['n50'],
                    genome_size=clean_data['genome_size'],
                    contig_count=clean_data['contig_count'],
                    confidence=clean_data['confidence'],
                    warnings=warnings_list,
                    ensemble_preds=preds_dict
                )
                self.add_row(row)