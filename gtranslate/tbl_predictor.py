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

from gtranslate.config.common import CONFIG
from gtranslate.config.output import *
from gtranslate.external.prodigal import Prodigal
from gtranslate.files.prodigal.tln_table_summary import TranslationSummaryFile, TranslationSummaryFileRow
from gtranslate.files.featurefile import FeatureFile
from gtranslate.tools import tqdm_log, symlink_f


class TablePredictor(object):
    """Predict the genetic translation table (GTT) used in prokaryotic organisms."""

    def __init__(self, cpus=1,custom_model_path=None,debug=False):
        """Initialize."""

        self.logger = logging.getLogger('timestamp')
        self.warnings = logging.getLogger('warnings')

        self.genome_file_suffix = GENOME_FILE_SUFFIX
        self.protein_file_suffix = PROTEIN_FILE_SUFFIX
        self.nt_gene_file_suffix = NT_GENE_FILE_SUFFIX
        self.gff_file_suffix = GFF_FILE_SUFFIX
        self.checksum_suffix = CHECKSUM_SUFFIX

        self.custom_model_path = custom_model_path

        self.cpus = cpus
        self.debug = debug

    def predict(self, genomes, out_dir, prefix, force):
            """Call Prodigal with TT 4 and 11 to identify genome profiles. Run Jellyfish to calculate kmers.
            Calculate the probability of each genome to be translated with TT 4/11 and TT 4/25.

            Parameters
            ----------
            genomes : dict
                Genome IDs as the key, path to genome file as value.
            out_dir : str
                Path to the output directory.
            prefix : str
                Prefix to append to generated files.
            force : bool
                Overwrite any existing files.

            Returns
            -------
            bool
                True if the process completes successfully.
            """

            reports = {}
            self.called_gene_dir = os.path.join(out_dir, DIR_PREDICT_GENES)
            self.failed_genomes = os.path.join(
                out_dir, PATH_FAILS.format(prefix=prefix))
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

            genome_dictionary = prodigal.run(genomes,self.custom_model_path)


            gene_files = [(db_genome_id, genome_dictionary[db_genome_id]['aa_gene_path'])
                          for db_genome_id in genome_dictionary.keys()]
            # if gene_files is empty, we have no genomes to process
            # we still need to write all the genomes failing the identify step to the bac120_summary file
            if len(gene_files) == 0:
                self.logger.warning('All genomes failed the identify step.')
                return reports


            self.logger.log(CONFIG.LOG_TASK,
                            'Summarising translation tables.')
            self._report_identified_translation_table(genome_dictionary, out_dir, prefix,
                                                            reports)

            return True

    def _report_identified_translation_table(self, gene_dict, outdir, prefix, reports):
        """Write the identified translation tables to disk."""

        # Summarise the copy number of each AR53 and BAC120 markers.
        tln_summary_file = TranslationSummaryFile(outdir, prefix)
        feature_file = FeatureFile(outdir, prefix)

        for db_genome_id, info in tqdm_log(sorted(gene_dict.items()), unit='genome'):
            # Write the best translation table to disk for this genome.
            summary_row = TranslationSummaryFileRow(gid=db_genome_id)
            summary_row.best_tln_table = info.get("best_translation_table")
            summary_row.coding_density_4 = info.get("coding_density_4")
            summary_row.coding_density_11 = info.get("coding_density_11")
            summary_row.gc_percent = info.get("gc_percent")
            summary_row.n50 = info.get("n50")
            summary_row.genome_size = info.get("genome_size")
            summary_row.contig_count = info.get("contig_count")
            summary_row.confidence=info.get("confidence")
            # we add the warnings if there are any
            warnings_list = info.get("warnings", [])
            if len(warnings_list) > 0:
                summary_row.warnings = warnings_list

            ensemble_preds_dict = info.get("ensemble_preds", {})
            if ensemble_preds_dict:
                summary_row.ensemble_preds = ensemble_preds_dict

            tln_summary_file.add_row(summary_row)

            features = info.get("feature_vector", {})
            if features:
                feature_file.add_row(db_genome_id, info.get("best_translation_table"), features)

        tln_summary_file.write()
        feature_file.write()

        #we remove the the outdir from the path to get a relative path
        relative_path_tln_summary_file = os.path.relpath(tln_summary_file.path, outdir)
        relative_path_feature_file = os.path.relpath(feature_file.path, outdir)
        symlink_f(relative_path_tln_summary_file,os.path.join(outdir,os.path.basename(tln_summary_file.path)),force=True)
        symlink_f(relative_path_feature_file,os.path.join(outdir,os.path.basename(feature_file.path)),force=True)

        return True
