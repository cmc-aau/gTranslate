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
import multiprocessing as mp
import os
import shutil
import subprocess

from gtranslate.biolib_lite.prodigal_runner import Prodigal as ProdigalRunner
from gtranslate.config.output import CHECKSUM_SUFFIX
from gtranslate.exceptions import ProdigalException
from gtranslate.files.prodigal.tln_table import TlnTableFile
from gtranslate.tools import sha256, file_has_checksum, tqdm_log



class Prodigal(object):
    """Perform ab initio gene prediction using Prodigal."""

    def __init__(self,
                 threads,
                 failed_genomes_file,
                 called_genes_dir,
                 protein_file_suffix,
                 nt_gene_file_suffix,
                 gff_file_suffix,
                 force):
        """Initialize."""

        self.logger = logging.getLogger('timestamp')
        self.warnings = logging.getLogger('warnings')
        self.failed_genomes_file = failed_genomes_file
        self.called_genes_dir = called_genes_dir

        self.threads = threads

        self.protein_file_suffix = protein_file_suffix
        self.nt_gene_file_suffix = nt_gene_file_suffix
        self.gff_file_suffix = gff_file_suffix
        self.force = force
        self.version = self._get_version()

    def _get_version(self):
        try:
            env = os.environ.copy()
            proc = subprocess.Popen(['prodigal', '-v'], stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, env=env, encoding='utf-8')

            output, error = proc.communicate()

            return error.split('\n')[1].split()[1].replace(':', '')
        except:
            return "(version unavailable)"

    def _run_prodigal(self, genome_id, fasta_path, custom_model_path):
        """Run Prodigal.

        Parameters
        ----------
        fasta_path : str
            Path to FASTA file to process.
        :return
            False if an error occurred.
        """

        # Create objects to write the output data to.
        output_dir = os.path.join(self.called_genes_dir, genome_id)
        tln_table_file = TlnTableFile(output_dir, genome_id)

        # Set the paths for output files.
        aa_gene_file = os.path.join(output_dir, genome_id + self.protein_file_suffix)
        nt_gene_file = os.path.join(output_dir, genome_id + self.nt_gene_file_suffix)
        gff_file = os.path.join(output_dir, genome_id + self.gff_file_suffix)
        out_files = (nt_gene_file, gff_file, tln_table_file.path, aa_gene_file)

        # Check if this genome has already been processed (skip).
        if all([file_has_checksum(x) for x in out_files]):
            tln_table_file.read()
            self.warnings.info(f'Skipped Prodigal processing for: {genome_id}')
            return aa_gene_file, nt_gene_file, gff_file, tln_table_file, True , False

        # Run Prodigal
        pr = ProdigalRunner(1, False)
        summary_stats = pr.run([(genome_id, fasta_path, custom_model_path)], output_dir, called_genes=False)

        # An error occurred in BioLib Prodigal.
        if not summary_stats:
            if self.force:
                return None
            else:
                raise Exception("An error was encountered while running Prodigal.")

        summary_stats = summary_stats[list(summary_stats.keys())[0]]

        if summary_stats.is_empty:
            shutil.rmtree(output_dir)
            if self.force:
                return summary_stats
            else:
                raise Exception("An error was encountered while running Prodigal.")

        else:
            # rename output files to adhere to GTDB conventions and desired genome
            # ID
            shutil.move(summary_stats.aa_gene_file, aa_gene_file)
            shutil.move(summary_stats.nt_gene_file, nt_gene_file)
            shutil.move(summary_stats.gff_file, gff_file)

            # save translation table information
            tln_table_file.best_tln_table = summary_stats.best_tln_table
            tln_table_file.coding_density_4 = summary_stats.coding_density_4
            tln_table_file.coding_density_11 = summary_stats.coding_density_11
            tln_table_file.gc_percent = summary_stats.gc_percent
            tln_table_file.n50 = summary_stats.n50
            tln_table_file.genome_size = summary_stats.genome_size
            tln_table_file.contig_count = summary_stats.contig_count
            tln_table_file.warnings = summary_stats.pred_warnings
            tln_table_file.confidence = summary_stats.pred_confidence
            tln_table_file.ensemble_preds = summary_stats.ensemble_preds
            tln_table_file.feature_vector = summary_stats.feature_vector
            tln_table_file.write()

            # create a hash of each file
            for out_file in out_files:
                if out_file is not None:
                    with open(out_file + CHECKSUM_SUFFIX, 'w') as fh:
                        fh.write(sha256(out_file))

        return aa_gene_file, nt_gene_file, gff_file, tln_table_file, False , summary_stats.is_empty

    def _worker(self, out_dict, worker_queue, writer_queue, n_skipped):
        """This worker function is invoked in a process."""

        while True:
            data = worker_queue.get(block=True, timeout=None)
            if data is None:
                break

            genome_id, file_path,custom_model_path = data

            rtn_files = self._run_prodigal(genome_id, file_path,custom_model_path)
            was_skipped = False

            # Only proceed if an error didn't occur in BioLib Prodigal
            if rtn_files:
                aa_gene_file, nt_gene_file, gff_file, translation_table_file,was_skipped,is_empty = rtn_files

                if was_skipped:
                    with n_skipped.get_lock():
                        n_skipped.value += 1

                prodigal_infos = {"aa_gene_path": aa_gene_file,
                                    "nt_gene_path": nt_gene_file,
                                    "gff_path": gff_file,
                                    "translation_table_path": translation_table_file.path,
                                    "best_translation_table": translation_table_file.best_tln_table,
                                    "coding_density_4": translation_table_file.coding_density_4,
                                    "coding_density_11": translation_table_file.coding_density_11,
                                    "gc_percent": translation_table_file.gc_percent,
                                    "n50": translation_table_file.n50,
                                    "genome_size": translation_table_file.genome_size,
                                    "contig_count": translation_table_file.contig_count,
                                    "confidence": translation_table_file.confidence,
                                    "warnings": translation_table_file.warnings,
                                    "ensemble_preds": translation_table_file.ensemble_preds,
                                    "feature_vector": translation_table_file.feature_vector,
                                    "is_empty": is_empty}

                out_dict[genome_id] = prodigal_infos

            writer_queue.put((genome_id,was_skipped))

    def _writer(self, num_items, writer_queue):
        """Store or write results of worker threads in a single thread."""
        with tqdm_log(total=num_items, unit='genome') as p_bar:
            for item in iter(writer_queue.get, None):  # Receive each item from the queue
                p_bar.update()
                try:
                    genome_id, was_skipped = item  # Unpack the tuple (genome_id, was_skipped)
                    if was_skipped:
                        p_bar.set_postfix({"skipped": f"{genome_id}"}, refresh=True)  # Update progress bar with skipped message
                except ValueError:
                    pass  # Handle unexpected item structures gracefully

    def run(self, genomic_files, custom_model_path=None):
        """Run Prodigal across a set of genomes.

        Parameters
        ----------
        genomic_files : dict
            Dictionary indicating the genomic and gene file for each genome.
        tln_tables : Dict[str, int]
            Mapping of genome id to user-specified translation table.
        """

        # populate worker queue with data to process
        worker_queue = mp.Queue()
        writer_queue = mp.Queue()

        if custom_model_path:
            self.logger.warning(f"Custom model path provided. Overriding default models with custom models from {custom_model_path}.")

        for genome_id, file_path in genomic_files.items():
            worker_queue.put((genome_id, file_path, custom_model_path))

        for _ in range(self.threads):
            worker_queue.put(None)

        worker_proc = []
        writer_proc = None
        try:
            manager = mp.Manager()
            out_dict = manager.dict()
            n_skipped = mp.Value('i', 0)

            worker_proc = [mp.Process(target=self._worker, args=(out_dict,
                                                                 worker_queue,
                                                                 writer_queue,
                                                                 n_skipped))
                           for _ in range(self.threads)]
            writer_proc = mp.Process(target=self._writer, args=(len(genomic_files),
                                                                writer_queue))

            writer_proc.start()
            for p in worker_proc:
                p.start()

            for p in worker_proc:
                p.join()

                # Gracefully terminate the program.
                if p.exitcode != 0:
                    raise ProdigalException('Prodigal returned a non-zero exit code.')

            writer_queue.put(None)
            writer_proc.join()
        except Exception as e:
            for p in worker_proc:
                p.terminate()
            if writer_proc:
                writer_proc.terminate()
            raise ProdigalException(f'An exception was caught while running Prodigal: {e}')

        # Report if any genomes were skipped due to having already been processed.
        if n_skipped.value > 0:
            genome_s = 'genome' if n_skipped.value == 1 else 'genomes'
            self.logger.warning(f'Prodigal skipped {n_skipped.value} {genome_s} '
                                f'due to pre-existing data, see warnings.log')


        # Report on any genomes which failed to have any genes called
        result_dict = dict()
        lq_gids = list()
        empty_gids = list()
        fails = open(self.failed_genomes_file,'w')
        for gid, gid_dict in out_dict.items():
            if gid_dict['is_empty']:
                empty_gids.append(gid)
            elif os.path.getsize(gid_dict['aa_gene_path']) <= 1:
                lq_gids.append(gid)
            else:
                result_dict[gid] = gid_dict

        if len(lq_gids) > 0:
            self.logger.warning(f'Skipping {len(lq_gids+empty_gids)} of {len(genomic_files)} '
                                f'genomes as no genes were called by Prodigal. '
                                f'Check the genome quality (see gtdbtk.warnings.log).')
            self.warnings.warning(f'The following {len(lq_gids+empty_gids)} genomes have '
                                  f'been excluded from analysis due to Prodigal '
                                  f'failing to call any genes:')


            # If there are few low-quality genomes just output to console.
            if len(lq_gids+empty_gids) > 10:
                for lq_gid in lq_gids:
                    self.warnings.info(lq_gid)
                    fails.write(f'{lq_gid}\tNo genes were called by Prodigal\n')
                for empty_gid in empty_gids:
                    self.warnings.info(empty_gid)
                    fails.write(f'{empty_gid}\tEmpty file\n')
            else:
                for lq_gid in lq_gids:
                    self.logger.warning(f'Skipping: {lq_gid}')
                    self.warnings.info(lq_gid)
                    fails.write(f'{lq_gid}\tNo genes were called by Prodigal\n')
                for empty_gid in empty_gids:
                    self.logger.warning(f'Skipping: {empty_gid}')
                    self.warnings.info(empty_gid)
                    fails.write(f'{empty_gid}\tEmpty file\n')

        fails.close()

        return result_dict
