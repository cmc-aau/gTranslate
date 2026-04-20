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

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2014'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'

import logging
import ntpath
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .common import make_sure_path_exists, check_file_exists
from .execute import check_on_path
from .parallel import Parallel
from .seq_io import read_fasta, write_fasta
from .seq_tk import N50, calculate_gc_content
from ..classifiers.ensemble_predictor import TTPredictor


@dataclass
class ConsumerData:
    aa_gene_file: str
    nt_gene_file: str
    gff_file: str
    is_empty: bool
    pred_confidence: float
    pred_warnings: list
    ensemble_preds: dict
    feature_vector: dict
    metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        for key, value in self.metadata.items():
            setattr(self, key, value) # Dynamically add metadata keys as attributes


class Prodigal(object):
    """Wrapper for running Prodigal in parallel."""

    def __init__(self, cpus, verbose=True):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use.
        verbose : boolean
            Flag indicating if progress should be reported.
        """

        self.logger = logging.getLogger('timestamp')

        check_on_path('prodigal')

        self.cpus = cpus
        self.verbose = verbose

    def _producer(self, genome_file_tuple):
        """Apply prodigal to genome with most suitable translation table.

        Parameters
        ----------
        genome_file : str
            Fasta file for genome.
        """

        genome_id,genome_file,custom_model_path = genome_file_tuple

        aa_gene_file = os.path.join(self.output_dir, genome_id + '_genes.faa')
        nt_gene_file = os.path.join(self.output_dir, genome_id + '_genes.fna')
        gff_file = os.path.join(self.output_dir, genome_id + '.gff')

        best_translation_table = -1
        table_coding_density = {4: -1, 11: -1}
        gene_numbers = {4: 0, 11: 0, 25: 0}

        # New dictionary to store Tryptophan AND Glycine counts
        table_trp_counts = {4: {'UGA': 0, 'UGG': 0, 'GLY': 0},
                            11: {'UGA': 0, 'UGG': 0, 'GLY': 0},
                            25: {'UGA': 0, 'UGG': 0, 'GLY': 0}}

        if self.called_genes:
            os.system('cp %s %s' %
                      (os.path.abspath(genome_file), aa_gene_file))
        else:
            seqs = read_fasta(genome_file)

            if len(seqs) == 0:
                self.logger.warning('Cannot call Prodigal on an empty genome. '
                                    'Skipped: {}'.format(genome_file))
                return (genome_id, aa_gene_file, nt_gene_file, gff_file, {}, 0.0, ["Skipped: Empty genome"], True)

            with tempfile.TemporaryDirectory('gtranslate_prodigal_tmp_') as tmp_dir:
                # if this is a gzipped genome, re-write the uncompressed genome
                # file to disk
                prodigal_input = genome_file
                if genome_file.endswith('.gz'):
                    prodigal_input = os.path.join(tmp_dir, os.path.basename(genome_file[0:-3]) + '.fna')
                    write_fasta(seqs, prodigal_input)

                # there may be ^M character in the input file,
                # the following code is similar to dos2unix command to remove
                # those special characters.
                with open(prodigal_input, 'r') as fh:
                    text = fh.read().replace('\r\n', '\n')
                processed_prodigal_input = os.path.join(
                    tmp_dir, os.path.basename(prodigal_input))
                with open(processed_prodigal_input, 'w') as fh:
                    fh.write(text)

                # determine number of bases
                total_bases = 0
                for seq in seqs.values():
                    total_bases += len(seq)

                # call genes under different translation tables
                translation_tables = [4, 11]
                local_warnings = []
                for translation_table in translation_tables:
                    os.makedirs(os.path.join(tmp_dir, str(translation_table)))
                    _aa_gene_file_tmp, nt_gene_file_tmp, gff_file_tmp, fallback_warning = self.run_prodigal_command(translation_table,
                                                                                                 tmp_dir, 
                                                                                                 genome_id,
                                                                                                 processed_prodigal_input, 
                                                                                                 total_bases)

                    if fallback_warning:
                        local_warnings.append(fallback_warning)



                    # determine coding density
                    prodigalParser = ProdigalGeneFeatureParser(gff_file_tmp)

                    coding_bases = 0
                    for seq_id, _seq in seqs.items():
                        coding_bases += prodigalParser.coding_bases(seq_id)

                    coding_density = float(coding_bases) / total_bases
                    table_coding_density[translation_table] = coding_density * 100

                    # --- NEW LOGIC: Calculate In-Frame UGA/UGG/GLY Counts ---
                    uga, ugg, gly = self._count_codons_of_interest(nt_gene_file_tmp)
                    table_trp_counts[translation_table]['UGA'] = uga
                    table_trp_counts[translation_table]['UGG'] = ugg
                    table_trp_counts[translation_table]['GLY'] = gly

                genome_metadata_dict = {}
                genome_metadata_dict['gc_percent'] = calculate_gc_content(seqs)
                genome_metadata_dict['n50'] = N50(seqs)
                genome_metadata_dict['genome_size'] = total_bases
                genome_metadata_dict['contig_count'] = len(seqs)

                # we store the coding density in the genome metadata dictionary
                genome_metadata_dict['coding_density_4'] = table_coding_density[4]
                genome_metadata_dict['coding_density_11'] = table_coding_density[11]

                # create a dataframe to store the genome information with columns in the same order as the classifiers
                # get columns from the scaler
                temp_df = pd.DataFrame( columns = ['GC', 'Coding_density_4','Coding_density_11','Density_Diff'])
                temp_df['Coding_density_4'] = [table_coding_density[4]]
                temp_df['Coding_density_11'] = [table_coding_density[11]]
                temp_df['Density_Diff'] = [table_coding_density[4] - table_coding_density[11]]
                temp_df['GC'] = [genome_metadata_dict['gc_percent']]
                #we add glycine information
                raw_ratio= np.log((table_trp_counts[4]['UGA'] + 1) / ((table_trp_counts[4]['UGG'] + 1)))
                #lets clip the trp ratio to be between -6 and 5 as in the training data
                temp_df['Trp_ratio'] = [max(-6.0, min(5.0, raw_ratio))]

                temp_df['Trp_magnitude'] = [np.log(table_trp_counts[4]['UGA'] + table_trp_counts[4]['UGG'] + 1)]

                raw_gly_ratio = np.log((table_trp_counts[4]['UGA'] + 1) / (table_trp_counts[4]['GLY'] + 1))
                temp_df['Gly_ratio'] = [max(-10.0, min(0.0, raw_gly_ratio))]

                temp_df['UGG_density'] = table_trp_counts[4]['UGG'] / (table_trp_counts[4]['GLY'])

                # predict best translation table
                predictor = TTPredictor(custom_model_path)

                best_translation_table,pred_confidence,pred_warnings,ensemble_preds,feature_vector = predictor.predict_translation_table(temp_df)
                best_translation_table = int(best_translation_table)

                pred_warnings.extend(local_warnings)

                if best_translation_table not in [4, 11]:
                    # Create the missing temporary folder for this specific table (e.g., /25/)
                    os.makedirs(os.path.join(tmp_dir, str(best_translation_table)), exist_ok=True)

                    # Actually run Prodigal to generate the .faa, .fna, and .gff files
                    _, _, _, fallback_warning = self.run_prodigal_command(
                        best_translation_table,
                        tmp_dir,
                        genome_id,
                        genome_file,
                        total_bases
                    )

                    if fallback_warning:
                        pred_warnings.extend(local_warnings)

                genome_metadata_dict['best_tln_table'] = best_translation_table

                shutil.copyfile(os.path.join(tmp_dir, str(best_translation_table),
                                             genome_id + '_genes.faa'), aa_gene_file)
                shutil.copyfile(os.path.join(tmp_dir, str(best_translation_table),
                                             genome_id + '_genes.fna'), nt_gene_file)
                shutil.copyfile(os.path.join(tmp_dir, str(best_translation_table),
                                             genome_id + '.gff'), gff_file)

        return (genome_id, aa_gene_file, nt_gene_file, gff_file,genome_metadata_dict,
                pred_confidence,pred_warnings,ensemble_preds,feature_vector,False)

    def _consumer(self, produced_data, consumer_data):
        """Consume results from producer processes.

         Parameters
        ----------
        produced_data : tuple
            Summary statistics for called genes for a specific genome.
        consumer_data : list
            Summary statistics of called genes for each genome.

        Returns
        -------
        consumer_data: d[genome_id] -> namedtuple(aa_gene_file,
                                                    nt_gene_file,
                                                    gff_file,
                                                    best_translation_table,
                                                    coding_density_4,
                                                    coding_density_11)
            Summary statistics of called genes for each genome.
        """

        if consumer_data is None:
            consumer_data = {}

        genome_id, aa_gene_file, nt_gene_file, gff_file, metadata_dict, pred_confidence, pred_warnings,ensemble_preds,feature_vector, is_empty = produced_data


        for warning in pred_warnings:
            # Check for our specific string so we only log the Prodigal fallbacks
            # (or remove the if-statement to log ALL predictor warnings too)
            if "Used Prodigal 'meta' mode fallback" in warning:
                self.logger.warning(f"{warning} for genome {genome_id}.")


        # FIX: Use keyword arguments to guarantee the data goes to the correct dataclass fields
        consumer_data[genome_id] = ConsumerData(
            aa_gene_file=aa_gene_file,
            nt_gene_file=nt_gene_file,
            gff_file=gff_file,
            is_empty=is_empty,
            pred_confidence=pred_confidence,
            pred_warnings=pred_warnings,
            ensemble_preds=ensemble_preds,
            metadata=metadata_dict,
            feature_vector=feature_vector
        )

        return consumer_data

    def _progress(self, processed_items, total_items):
        """Report progress of consumer processes.

        Parameters
        ----------
        processed_items : int
            Number of genomes processed.
        total_items : int
            Total number of genomes to process.

        Returns
        -------
        str
            String indicating progress of data processing.
        """

        return self.progress_str % (processed_items, total_items, float(processed_items) * 100 / total_items)

    def run(self,
            genome_files,
            output_dir,
            called_genes=False,
            translation_table=None,
            meta=False,
            closed_ends=False):
        """Call genes with Prodigal.

        Call genes with prodigal and store the results in the
        specified output directory. For convenience, the
        called_gene flag can be used to indicate genes have
        previously been called and simply need to be copied
        to the specified output directory.

        Parameters
        ----------
        genome_files : list of str
            Nucleotide fasta files to call genes on.
        called_genes : boolean
            Flag indicating if genes are already called.
        translation_table : int
            Specifies desired translation table, use None to automatically
            select between tables 4 and 11.
        meta : boolean
            Flag indicating if prodigal should call genes with the metagenomics procedure.
        closed_ends : boolean
            If True, do not allow genes to run off edges (throws -c flag).
        output_dir : str
            Directory to store called genes.

        Returns
        -------
        d[genome_id] -> namedtuple(best_translation_table
                                            coding_density_4
                                            coding_density_11)
            Summary statistics of called genes for each genome.
        """

        self.called_genes = called_genes
        self.translation_table = translation_table
        self.meta = meta
        self.closed_ends = closed_ends
        self.output_dir = output_dir

        make_sure_path_exists(self.output_dir)

        progress_func = None
        if self.verbose:
            file_type = 'genomes'
            self.progress_str = '  Finished processing %d of %d (%.2f%%) genomes.'
            if meta:
                file_type = 'scaffolds'
                if len(genome_files):
                    file_type = ntpath.basename(genome_files[0][1])

                self.progress_str = '  Finished processing %d of %d (%.2f%%) files.'

            self.logger.info('Identifying genes within %s: ' % file_type)
            progress_func = self._progress

        parallel = Parallel(self.cpus)
        summary_stats = parallel.run(
            self._producer, self._consumer, genome_files, progress_func)

        # An error was encountered during Prodigal processing, clean up.
        if not summary_stats:
            shutil.rmtree(self.output_dir)

        return summary_stats


    def run_prodigal_command(self, translation_table, tmp_dir, genome_id, genome_file, total_bases):
        aa_gene_file_tmp = os.path.join(tmp_dir, str(
            translation_table), genome_id + '_genes.faa')
        nt_gene_file_tmp = os.path.join(tmp_dir, str(
            translation_table), genome_id + '_genes.fna')
        gff_file_tmp = os.path.join(tmp_dir, str(
            translation_table), genome_id + '.gff')

        # check if there is sufficient bases to calculate prodigal
        # parameters
        if total_bases < 100000 or self.meta:
            proc_str = 'meta'  # use best pre-calculated parameters
        else:
            proc_str = 'single'  # estimate parameters from data

        args = '-m'
        if self.closed_ends:
            args += ' -c'

        cmd = ['prodigal', args, '-p', proc_str, '-q',
               '-f', 'gff', '-g', str(translation_table),
               '-a', aa_gene_file_tmp, '-d', nt_gene_file_tmp,
               '-i', genome_file, '-o', gff_file_tmp]

        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, encoding='utf-8')

        stdout, stderr = proc.communicate()
        fallback_warning = None
        if proc.returncode != 0 and proc_str == 'single':
            fallback_warning = f"Used Prodigal 'meta' mode fallback for TT{translation_table}"
            #self.logger.warning(f"Prodigal 'single' mode failed for {genome_id}. Retrying with 'meta' mode...")
            cmd[cmd.index('single')] = 'meta'

            # Re-run with the modified command
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
            stdout, stderr = proc.communicate()

        if proc.returncode != 0:
            self.logger.warning('Error running Prodigal on genome: '
                                '{}'.format(genome_file))
            self.logger.warning('Error message:')
            for line in stderr.splitlines():
                print(line)
            self.logger.warning('This genome is skipped.')

        return aa_gene_file_tmp, nt_gene_file_tmp, gff_file_tmp, fallback_warning

    def _count_codons_of_interest(self, gene_file):
        """
        Count in-frame Trp (TGA, TGG) and Standard Glycine (GGN) codons.

        Parameters
        ----------
        gene_file : str
            Path to the nucleotide gene FASTA file produced by Prodigal.

        Returns
        -------
        tuple (int, int, int)
            (count_of_TGA, count_of_TGG, count_of_Std_Gly)
        """
        uga_count = 0
        ugg_count = 0
        gly_count = 0  # Sum of GGA, GGC, GGG, GGT

        # Parse the gene file (standard FASTA)
        try:
            genes = read_fasta(gene_file)
        except Exception:
            return 0, 0, 0

        for seq in genes.values():
            seq = seq.upper()
            # Iterate over the sequence in steps of 3 to ensure we are reading in-frame codons
            # Prodigal output (-d) is already aligned to the start codon.
            for i in range(0, len(seq), 3):
                codon = seq[i:i + 3]
                if len(codon) < 3:
                    continue

                # Check Tryptophan Candidates
                if codon == 'TGA':
                    uga_count += 1
                elif codon == 'TGG':
                    ugg_count += 1

                # Check Standard Glycine Candidates (GGN)
                elif codon in ['GGA', 'GGC', 'GGG', 'GGT']:
                    gly_count += 1

        return uga_count, ugg_count, gly_count


class ProdigalGeneFeatureParser(object):
    """Parses prodigal gene feature files (GFF) output."""

    def __init__(self, filename):
        """Initialization.

        Parameters
        ----------
        filename : str
            GFF file to parse.
        """
        check_file_exists(filename)

        self.genes = {}
        self.last_coding_base = {}

        self.__parseGFF(filename)

        self.coding_base_masks = {}
        for seq_id in self.genes:
            self.coding_base_masks[seq_id] = self.__build_coding_base_mask(
                seq_id)

    def __parseGFF(self, filename):
        """Parse genes from GFF file.

        Parameters
        ----------
        filename : str
            GFF file to parse.
        """
        bGetTranslationTable = True
        with open(filename, 'r') as fh:
            for line in fh:
                if bGetTranslationTable and line.startswith('# Model Data'):
                    data_model_info = line.split(':')[1].strip().split(';')
                    dict_data_model = {}
                    for item in data_model_info:
                        k = item.split('=')[0]
                        v = item.split('=')[1]
                        dict_data_model[k] = v

                    self.translationTable = int(
                        dict_data_model.get('transl_table'))
                    bGetTranslationTable = False

                if line[0] == '#':
                    continue

                line_split = line.split('\t')
                seq_id = line_split[0]
                if seq_id not in self.genes:
                    geneCounter = 0
                    self.genes[seq_id] = {}
                    self.last_coding_base[seq_id] = 0

                geneId = seq_id + '_' + str(geneCounter)
                geneCounter += 1
                start = int(line_split[3])
                end = int(line_split[4])

                self.genes[seq_id][geneId] = [start, end]
                self.last_coding_base[seq_id] = max(
                    self.last_coding_base[seq_id], end)

    def __build_coding_base_mask(self, seq_id):
        """Build mask indicating which bases in a sequences are coding.

        Parameters
        ----------
        seq_id : str
            Unique id of sequence.
        """

        # safe way to calculate coding bases as it accounts
        # for the potential of overlapping genes
        coding_base_mask = np.zeros(self.last_coding_base[seq_id], dtype=bool)
        for pos in self.genes[seq_id].values():
            coding_base_mask[pos[0]:pos[1] + 1] = True

        return coding_base_mask

    def coding_bases(self, seq_id, start=0, end=None):
        """Calculate number of coding bases in sequence between [start, end).

        To process the entire sequence set start to 0, and
        end to None.

        Parameters
        ----------
        seq_id : str
            Unique id of sequence.
        start : int
            Start calculation at this position in sequence.
        end : int
            End calculation just before this position in the sequence.
        """

        # check if sequence has any genes
        if seq_id not in self.genes:
            return 0

        # set end to last coding base if not specified
        if end is None:
            end = self.last_coding_base[seq_id]

        return np.sum(self.coding_base_masks[seq_id][start:end])
