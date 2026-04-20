import json
import os
from typing import Optional, List

from gtranslate.config.output import TRANSLATION_TABLE_SUFFIX
from gtranslate.exceptions import GTranslateExit


class TlnTableFile(object):
    """
    A class to handle the translation table summary file for a genome.
    """

    def __init__(self, out_dir: str, gid: str,
                 best_tln_table: Optional[int] = None,
                 coding_density_4: Optional[float] = None,
                 coding_density_11: Optional[float] = None,
                 gc_percentage: Optional[float] = None,
                 n50_value: Optional[int] = None,
                 genome_length: Optional[int] = None,
                 contig_count: Optional[int] = None,
                 confidence: Optional[float] = None,
                 ensemble_preds: Optional[dict] = None,
                 feature_vector:Optional[dict] = None,
                 warnings: Optional[List[str]] = None):

        self.path = self.get_path(out_dir, gid)
        self._best_tln_table = best_tln_table
        self._coding_density_4 = coding_density_4
        self._coding_density_11 = coding_density_11
        self._gc_percent = gc_percentage
        self._n50 = n50_value
        self._genome_size = genome_length
        self._contig_count = contig_count
        self._confidence = confidence
        self._ensemble_preds = ensemble_preds or {}
        self._feature_vector = feature_vector or {}
        self._warnings = warnings or []

    def _validate_and_set(self, attribute, value, expected_type):
        # 1. Handle missing/empty data gracefully
        if value is None or str(value).strip().upper() in ['N/A', 'NONE', 'NAN', '']:
            setattr(self, f'_{attribute}', None)
            return

        try:
            # If expected_type is int, casting a string like '4.0' directly
            # to int crashes. Casting to float first safely strips the decimal.
            if expected_type is int:
                clean_value = int(float(value))
            else:
                clean_value = expected_type(value)

            setattr(self, f'_{attribute}', clean_value)

        except ValueError:
            raise GTranslateExit(f'Invalid {attribute} value: {value} for {self.path}')

    @property
    def best_tln_table(self):
        return self._best_tln_table

    @best_tln_table.setter
    def best_tln_table(self, v):
        self._validate_and_set('best_tln_table', v, int)

    @property
    def coding_density_4(self):
        return self._coding_density_4

    @coding_density_4.setter
    def coding_density_4(self, v):
        self._validate_and_set('coding_density_4', v, float)

    @property
    def coding_density_11(self):
        return self._coding_density_11

    @coding_density_11.setter
    def coding_density_11(self, v):
        self._validate_and_set('coding_density_11', v, float)

    @property
    def gc_percent(self):
        return self._gc_percent

    @gc_percent.setter
    def gc_percent(self, v):
        self._validate_and_set('gc_percent', v, float)

    @property
    def n50(self):
        return self._n50

    @n50.setter
    def n50(self, v):
        self._validate_and_set('n50', v, int)

    @property
    def genome_size(self):
        return self._genome_size

    @genome_size.setter
    def genome_size(self, v):
        self._validate_and_set('genome_size', v, int)

    @property
    def contig_count(self):
        return self._contig_count

    @contig_count.setter
    def contig_count(self, v):
        self._validate_and_set('contig_count', v, int)

    @property
    def confidence(self):
        return self._confidence

    @confidence.setter
    def confidence(self, v):
        self._validate_and_set('confidence', v, float)

    @property
    def ensemble_preds(self):
        return self._ensemble_preds 

    @ensemble_preds.setter
    def ensemble_preds(self, v):
        if v is None or str(v).strip().upper() in ['N/A', 'NONE', 'NAN', '']:
            self._ensemble_preds = {}
        elif isinstance(v, dict):
            self._ensemble_preds = v
        else:
            raise GTranslateExit(f'Invalid ensemble_preds value: {v} for {self.path}')

    @property
    def feature_vector(self):
        return self._feature_vector

    @feature_vector.setter
    def feature_vector(self, v):
        if v is None or str(v).strip().upper() in ['N/A', 'NONE', 'NAN', '']:
            self._feature_vector = {}
        elif isinstance(v, dict):
            self._feature_vector = v
        else:
            raise GTranslateExit(f'Invalid feature_vector value: {v} for {self.path}')

    @property
    def warnings(self):
        return self._warnings

    @warnings.setter
    def warnings(self, v):
        # Custom logic to handle strings (from reading the file) or lists
        if v is None or str(v).strip().upper() in ['N/A', 'NONE', 'NAN', '']:
            self._warnings = []
        elif isinstance(v, str):
            self._warnings = [w.strip() for w in v.split(';') if w.strip()]
        elif isinstance(v, list):
            self._warnings = v
        else:
            raise GTranslateExit(f'Invalid warnings value: {v} for {self.path}')

    @staticmethod
    def get_path(out_dir: str, gid: str):
        return os.path.join(out_dir, f'{gid}{TRANSLATION_TABLE_SUFFIX}')

    def read(self):
        try:
            with open(self.path, 'r') as fh:
                for line in fh.readlines():
                    # Handle empty values cleanly by using .split('\t', 1)
                    parts = line.strip('\n').split('\t', 1)
                    if len(parts) != 2:
                        continue

                    idx, val = parts

                    if idx == 'best_translation_table':
                        self.best_tln_table = val
                    elif idx == 'coding_density_4':
                        self.coding_density_4 = val
                    elif idx == 'coding_density_11':
                        self.coding_density_11 = val
                    elif idx == 'gc_percent':
                        self.gc_percent = val
                    elif idx == 'n50':
                        self.n50 = val
                    elif idx == 'genome_size':
                        self.genome_size = val
                    elif idx == 'contig_count':
                        self.contig_count = val
                    elif idx == 'confidence':
                        self.confidence = val
                    elif idx == 'warnings':
                        self.warnings = val
                    elif idx == 'ensemble_preds':
                        if val.strip().upper() in ['N/A', 'NONE', 'NAN', '']:
                            self.ensemble_preds = {}
                        else:
                            self.ensemble_preds = json.loads(val)
                    elif idx == 'feature_vector':
                        if val.strip().upper() in ['N/A', 'NONE', 'NAN', '']:
                            self.feature_vector = {}
                        else:
                            self.feature_vector = json.loads(val)

        except FileNotFoundError:
            raise GTranslateExit(f'Translation table summary file not found: {self.path}')
        except Exception as e:
            raise GTranslateExit(f'Error parsing file: {self.path} - {e}')

    def write(self):
        try:
            with open(self.path, 'w') as fh:
                fh.write(f'best_translation_table\t{self.best_tln_table}\n')
                fh.write(f'coding_density_4\t{self.coding_density_4}\n')
                fh.write(f'coding_density_11\t{self.coding_density_11}\n')
                fh.write(f'gc_percent\t{self.gc_percent}\n')
                fh.write(f'n50\t{self.n50}\n')
                fh.write(f'genome_size\t{self.genome_size}\n')
                fh.write(f'contig_count\t{self.contig_count}\n')
                fh.write(f'confidence\t{self.confidence}\n')

                # Format warnings list to a semicolon-separated string, or 'N/A' if empty
                warnings_str = ';'.join(self.warnings) if self.warnings else 'N/A'
                fh.write(f'warnings\t{warnings_str}\n')

                ensemble_preds_str = json.dumps(self.ensemble_preds) if self.ensemble_preds else 'N/A'
                fh.write(f'ensemble_preds\t{ensemble_preds_str}\n')

                feature_vector_str = json.dumps(self.feature_vector) if self.feature_vector else 'N/A'
                fh.write(f'feature_vector\t{feature_vector_str}\n')


        except Exception as e:
            raise GTranslateExit(f'Error writing file: {self.path} - {e}')