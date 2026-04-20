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
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from tqdm import tqdm

from gtranslate.biolib_lite.execute import check_dependencies
from gtranslate.biolib_lite.common import remove_extension, check_dir_exists, check_file_exists, make_sure_path_exists
from gtranslate.exceptions import GTranslateExit
from gtranslate.files.batchfile import Batchfile
from gtranslate.files.prodigal.tln_table_summary import TranslationSummaryFile
from gtranslate.misc import Misc
from gtranslate.plots.plotter import FeaturePlotter
from gtranslate.tbl_predictor import TablePredictor
from gtranslate.tools import remove_intermediate_files
from gtranslate.training_manager import TrainingManager


class OptionsParser(object):

    def __init__(self, version,output_dir=None):
        """Initialization.

        Parameters
        ----------
        version : str
            The current version number (e.g. 0.2.2).
        """
        self.logger = logging.getLogger('timestamp')
        self.warnings = logging.getLogger('warnings')
        self.version = version

        self.genomes_to_process = None


    @staticmethod
    def _verify_file_path(file_path: str) -> bool:
        if ' ' in file_path:
            raise GTranslateExit(f'Error: The genome file path contains spaces:\n{file_path}\n'
                                 f'Rename the file or move it to a directory without spaces.')
        return True

    def _verify_genome_id(self, genome_id: str) -> bool:
        """Ensure genome ID will be valid in Newick tree.

        Parameters
        ----------
        genome_id : str
            The string representing the genome identifier.

        Returns
        -------
        bool
            True if the genome identifier is legal.

        Raises
        ------
        GTDBTkExit
            If the genome identifier contains illegal characters.
        """
        if len(genome_id) == 0:
            raise GTranslateExit('Genome name cannot be blank, check for input files '
                             'without a name, or empty columns in the batchfile.')
        invalid_chars = frozenset('()[],;= ')
        if any((c in invalid_chars) for c in genome_id):
            self.logger.error(f'Invalid genome ID: {genome_id}')
            self.logger.error(f'The following characters are invalid: '
                              f'{" ".join(invalid_chars)}')
            raise GTranslateExit(f'Invalid genome ID: {genome_id}')
        return True

    def _genomes_to_process(self, genome_dir, batchfile, extension):
        """Get genomes to process.

        Parameters
        ----------
        genome_dir : str
            Directory containing genomes.
        batchfile : str
            File describing genomes.
        extension : str
            Extension of files to process.

        Returns
        -------
        genomic_files : d[genome_id] -> FASTA file
            Map of genomes to their genomic FASTA files.
        """

        genomic_files = dict()
        if genome_dir:
            for f in os.listdir(genome_dir):
                if f.endswith(extension):
                    genome_id = remove_extension(f, extension)
                    genomic_files[genome_id] = os.path.join(genome_dir, f)

        elif batchfile:
            if not os.path.exists(batchfile):
                raise GTranslateExit(f'Batchfile does not exist: {batchfile}')
            batchfile_fh = Batchfile(batchfile)
            genomic_files = batchfile_fh.genome_path

        # Check that all of the genome IDs are valid.
        for genome_key in genomic_files:
            self._verify_genome_id(genome_key)

        # Check that there are no illegal characters in the file path
        for file_path in genomic_files.values():
            self._verify_file_path(file_path)

        # Check that the prefix is valid and the path exists
        invalid_paths = list()
        for genome_key, genome_path in tqdm(genomic_files.items(), desc="Checking paths",ncols=100):
            if not Path(genome_path).exists():
                invalid_paths.append((genome_key, genome_path))

        # Report on any invalid paths
        if len(invalid_paths) > 0:
            self.warnings.info(f'Reading from batchfile: {batchfile}')
            self.warnings.error(f'The following {len(invalid_paths)} genomes '
                                f'have invalid paths specified in the batchfile:')
            for g_path, g_gid in invalid_paths:
                self.warnings.info(f'{g_gid}\t{g_path}')
            raise GTranslateExit(f'There are {len(invalid_paths)} paths in the '
                             f'batchfile which do not exist, see gtranslate.warnings.log')

        if len(genomic_files) == 0:
            if genome_dir:
                self.logger.error(f'No genomes found in directory: {genome_dir}. '
                                  f'Check the --extension flag or verify that the directory contains the expected files.')
            else:
                self.logger.error(f'No genomes found in batch file: {batchfile}. '
                                  f'Ensure the batch file is formatted correctly.')
            raise GTranslateExit

        return genomic_files

    def detect_table(self, options):
        """Detect the genetic translation table (GTT) used in prokaryotic organisms.

        Parameters
        ----------
        options : argparse.Namespace
            The CLI arguments input by the user.
        """
        self.logger.info('Running detect_table')

        check_dependencies(['prodigal'], exit_on_fail=True)

        if options.genome_dir:
            check_dir_exists(options.genome_dir)

        if options.batchfile:
            check_file_exists(options.batchfile)

        make_sure_path_exists(options.out_dir)

        genomes = self._genomes_to_process(options.genome_dir,
                                                       options.batchfile,
                                                       options.extension)

        table_predictor = TablePredictor(options.cpus,options.custom_model_path)
        reports = table_predictor.predict(genomes,
                         options.out_dir,
                         options.prefix,
                         options.force)

        if not options.keep_called_genes:
            self.logger.info('Cleaning up intermediate files.')
            self.remove_intermediate_files(options.out_dir)

        self.logger.info('Done.')

    def generate_plot(self, options):
        """Generate an interactive HTML dashboard to explore the features used for GTT prediction.

        Parameters
        ----------
        options : argparse.Namespace
            The CLI arguments input by the user.
        """
        self.logger.info('Running generate_plot')

        plotter = FeaturePlotter(options.feature_file, options.output_file,options.selected_genome_file)
        plotter.generate_html()


        self.logger.info('Done.')

    def ground_truth(self, options):
        """Determine the ground truth for genomes based on their taxonomic classification.

        Parameters
        ----------
        options : argparse.Namespace
            The CLI arguments input by the user.
        """
        self.logger.info('Selecting Ground Truth translation tables based on taxonomic classification.')

        training_manager=TrainingManager()
        training_manager.select_ground_truth(options.taxonomy_file, options.output_file, options.manual_gt_file)

        self.logger.info('Done.')

    def build_features(self, options):
        """Generate feature vectors for training models.

        Parameters
        ----------
        options : argparse.Namespace
            The CLI arguments input by the user.
        """
        self.logger.info('Generating feature vectors for training models.')

        check_dependencies(['prodigal'], exit_on_fail=True)

        genomes = self._genomes_to_process(options.genome_dir,
                                            options.batchfile,
                                            options.extension)

        training_manager=TrainingManager(cpus=options.cpus)
        training_manager.build_features(genomes, options.out_dir,options.force)

        self.logger.info('Done.')

    def fit_models(self, options):
        """Train the models based on training data.

        Parameters
        ----------
        options : argparse.Namespace
            The CLI arguments input by the user.
        """
        self.logger.info('Training models based on training data.')

        training_manager=TrainingManager(cpus=options.cpus,seed=options.seed)

        training_manager.fit_models(options.feature_file, options.tt_file,options.out_dir,options.split_data)

        self.logger.info('Done.')



    def run_test(self, options):
        """Run the test suite. To make sure the program is working as expected.
        """

        # Use a temporary directory if none is supplied.
        if options.out_dir:
            out_dir_fh = None
            make_sure_path_exists(options.out_dir)
        else:
            out_dir_fh = tempfile.TemporaryDirectory(prefix='gtranslate_tmp_')
            options.out_dir = out_dir_fh.name
            self.logger.info('Using a temporary directory as out_dir was not specified.')

        try:
            output_dir = os.path.join(options.out_dir, 'output')
            genome_test_dir = os.path.join(options.out_dir, 'genomes')
            if os.path.exists(genome_test_dir):
                self.logger.error(f'Test directory {genome_test_dir} already exists.')
                self.logger.error('Test must be run in a new directory.')
                sys.exit(1)

            current_path = os.path.dirname(os.path.realpath(__file__))
            input_dir = os.path.join(current_path, 'tests', 'data', 'genomes')

            shutil.copytree(input_dir, genome_test_dir)

            args = ['gtranslate', 'detect_table', '--genome_dir', genome_test_dir,
                    '--out_dir', output_dir, '--cpus', str(options.cpus)]
            self.logger.info('Command: {}'.format(' '.join(args)))

            # Pipe the output and write to disk.
            path_stdout = os.path.join(options.out_dir, 'test_execution.log')
            with subprocess.Popen(args, stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE, encoding='utf-8') as proc:
                with open(path_stdout, 'w') as fh_stdout:
                    bar_fmt = ' <TEST OUTPUT> '.center(22) + '{desc}'
                    with tqdm(bar_format=bar_fmt, leave=False) as p_bar:
                        while True:
                            line = proc.stdout.readline()
                            if not line:
                                break
                            fh_stdout.write(f'{line}')
                            p_bar.set_description_str(line.strip())
                proc.wait()
                exit_code = proc.returncode

            summary_fh = TranslationSummaryFile(output_dir, 'gtranslate')

            if exit_code != 0:
                self.logger.error('The test returned a non-zero exit code.')
                self.logger.error('A detailed summary of the execution log can be '
                                  'found here: {}'.format(path_stdout))
                self.logger.error('The test has failed.')
                sys.exit(1)
            if not os.path.exists(summary_fh.path):
                self.logger.error(f"{summary_fh.path} is missing.")
                self.logger.error('A detailed summary of the execution log can be '
                                  'found here: {}'.format(path_stdout))
                self.logger.error('The test has failed.')
                sys.exit(1)
        finally:
            if out_dir_fh:
                out_dir_fh.cleanup()

        self.logger.info('Test has successfully finished.')
        return True

    def check_install(self,options):
        """ Verify all GTDB-Tk data files are present.

        Raises
        ------
        ReferenceFileMalformed
            If one or more reference files are malformed.
        """
        self.logger.info("Running install verification")
        misc = Misc()
        misc.check_install()
        self.logger.info('Done.')



    def remove_intermediate_files(self,out_dir):
        """Remove intermediate files from the output directory.
        Parameters
        ----------
            out_dir : str
                The output directory.
            workflow_name : str
                The name of the workflow that generated the intermediate files.
        """
        if not os.path.exists(out_dir):
            self.logger.warning(f'Output directory does not exist: {out_dir}. Skipping cleanup.')
            return
        self.logger.info('Deleting generated gene files.')
        remove_intermediate_files(out_dir)
        self.logger.info('gene files deleted.')

    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)

        Parameters
        ----------
        options : argparse.Namespace
            The CLI arguments input by the user.
        """

        # Stop processing if python 2 is being used.
        if sys.version_info.major < 3:
            raise GTranslateExit('Python 2 is no longer supported.')

        # Correct user paths
        if hasattr(options, 'out_dir') and options.out_dir:
            options.out_dir = os.path.expanduser(options.out_dir)

        # Assert that the number of CPUs is a positive integer.
        if hasattr(options, 'cpus') and options.cpus < 1:
            self.logger.warning(
                'You cannot use less than 1 CPU, defaulting to 1.')
            options.cpus = 1

        if options.subparser_name == 'detect_table':
            self.detect_table(options)
        elif options.subparser_name == 'generate_plot':
            self.generate_plot(options)
        elif options.subparser_name == 'test':
            self.run_test(options)
        elif options.subparser_name == 'check_install':
            self.check_install(options)
        elif options.subparser_name == 'ground_truth':
            self.ground_truth(options)
        elif options.subparser_name == 'build_features':
            self.build_features(options)
        elif options.subparser_name == 'fit_models':
            self.fit_models(options)
        else:
            raise GTranslateExit(f'Unknown gTranslate command: {options.subparser_name}')

        return 0
