#! /usr/bin/env python3

# Evaluate if misclassifications are likely due to an incorrect ground truth.

__prog_name__ = 'evaluate_misclassifications.py'
__prog_desc__ = 'Evaluate if misclassifications are likely due to an incorrect ground truth.'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2026'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.1.0'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'


import os
import sys
import logging
import argparse
import shutil
import subprocess
import concurrent.futures
from collections import defaultdict
from dataclasses import dataclass
from tqdm import tqdm

from gtranslate.biolib_lite.logger import logger_setup
from gtranslate.biolib_lite.execute import check_dependencies
from gtranslate.tools import canonical_gid


@dataclass
class GroundTruthData:
    translation_table: str
    gtdb_taxonomy: str
    ncbi_taxonomy: str


@dataclass
class gTranslateData:
    translation_table: str
    cd4: float
    cd11: float
    gc: float
    n50: int
    genome_size: int
    contig_count: int


@dataclass
class ProkkaResult:
    tRNA_Trp_tca: int
    RF2: int    


@dataclass
class CodettaResult:
    base: str
    columns: int


class EvaluateMisclassifications(object):
    """Evaluate if misclassifications are likely due to an incorrect ground truth."""

    def __init__(self):
        """Initialize."""

        check_dependencies(['prokka'])

        self.log = logging.getLogger('timestamp')

    def parse_ground_truth_file(self, ground_truth_file: str) -> dict:
        """Parse ground truth for genomes."""

        ground_truth = {}
        with open(ground_truth_file) as f:
            header = f.readline().strip().split('\t')

            gid_idx = header.index('Genome ID')
            gt_idx = header.index('Ground truth table')
            gtdb_idx = header.index('GTDB taxonomy')
            ncbi_idx = header.index('NCBI taxonomy')

            for line in f:
                tokens = line.strip().split('\t')

                gid = tokens[gid_idx]

                ground_truth[gid] = GroundTruthData(
                    tokens[gt_idx],
                    tokens[gtdb_idx],
                    tokens[ncbi_idx]
                )

        return ground_truth

    def parse_gtranslate_misclassifications(self, gtranslate_results_file: str, ground_truth: dict) -> dict:
        """Parse file with gTranslate results to find misclassifications relative to ground truth."""

        misclassifications = {}
        with open(gtranslate_results_file) as f:
            header = f.readline().strip().split('\t')

            gid_idx = header.index('user_genome')
            tt_idx = header.index('best_tln_table')
            cd4_ind = header.index('coding_density_4')
            cd11_idx = header.index('coding_density_11')
            gc_idx = header.index('gc_percent')
            n50_idx = header.index('n50')
            gs_idx = header.index('genome_size')
            contig_count_idx = header.index('contig_count')

            for line in f:
                tokens = line.strip().split('\t')

                gid = canonical_gid(tokens[gid_idx])
                tt = tokens[tt_idx]

                if gid not in ground_truth:
                    self.log.error(f'Genome not in ground truth: {gid}')
                    sys.exit(1)

                if ground_truth[gid].translation_table == tt:
                    # correctly predicted by gTranslate
                    continue

                misclassifications[gid] = gTranslateData(
                    tt,
                    float(tokens[cd4_ind]),
                    float(tokens[cd11_idx]),
                    float(tokens[gc_idx]),
                    int(tokens[n50_idx]),
                    int(tokens[gs_idx]),
                    int(tokens[contig_count_idx])
                ) 

        return misclassifications

    def parse_genome_path_file(self, genome_path_file: str, misclassified_gids: dict) -> dict:
        """Parse path to genomic FASTA file of putatively misclassified genomes."""

        genome_paths = {}
        with open(genome_path_file) as f:
            for line in f:
                tokens = line.strip().split('\t')

                gp = tokens[0]
                gid = tokens[1]
                if gid not in misclassified_gids:
                    continue

                genome_paths[gid] = gp

        return genome_paths

    def parse_prokka_results(self, gff_file: str):
        """Parse Prokka annotations for tRNA-Trp, RF1, and RF2, genes."""

        tRNA_Trp_tca = 0
        rf2 = 0
        with open(gff_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue

                tokens = line.strip().split('\t')
                stat = tokens[-1].split(';')


                if 'product=Peptide chain release factor 2' in stat or 'product=Peptide chain release factor RF2' in stat:
                    rf2 += 1
                elif 'product=tRNA-Trp(tca)' in stat:
                    tRNA_Trp_tca += 1
                  
        return ProkkaResult(tRNA_Trp_tca, rf2)

    def decompress_genomes(self, genome_paths: dict, misclassified_gids: dict, out_dir: str) -> dict:
        """Decompress putatively misclassified genomes."""

        genome_dir = os.path.join(out_dir, 'genomes')
        os.makedirs(genome_dir, exist_ok=True)

        genome_paths_decompressed = {}
        for gid in misclassified_gids:
            # copy and decompress genomic FASTA file
            decompressed_genome_file = os.path.join(genome_dir, f'{gid}.fna')
            if not os.path.exists(decompressed_genome_file):
                cur_genome_file = decompressed_genome_file + '.gz'
                shutil.copyfile(genome_paths[gid], cur_genome_file)
                subprocess.run(['pigz', '-d', cur_genome_file], check=True)
                genome_paths_decompressed[gid] = decompressed_genome_file

        return genome_paths_decompressed

    def _run_prokka_for_genome(self, gid: str, genome_path: str, prokka_dir: str) -> tuple:
        """Run Prokka for a single genome."""

        results = {}
        for tt in [4, 11]:
            # run prokka
            cmd = [
                'prokka', '--kingdom', 'Bacteria', '--outdir', prokka_dir,
                '--locustag', f'{gid}-tt{tt}', '--prefix', f'{gid}-tt{tt}',
                '--gcode', str(tt), genome_path,
                '--cpus', '8', '--force',
            ]
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            results[tt] = self.parse_prokka_results(os.path.join(prokka_dir, f'{gid}-tt{tt}.gff'))
            
        return gid, results

    def identify_uga_reassignment_properties(self, 
                                             genome_paths_decompressed: dict, 
                                             cpus: int,
                                             out_dir: str) -> dict:
        """Identify genomic properties of UGA reassignment using Prokka."""

        prokka_dir = os.path.join(out_dir, 'prokka')
        os.makedirs(prokka_dir, exist_ok=True)

        genomic_properties = defaultdict(dict)
        num_workers = max(1, cpus // 8)
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
            futures = [
                executor.submit(self._run_prokka_for_genome, gid, genome_paths_decompressed[gid], prokka_dir)
                for gid in genome_paths_decompressed
            ]
            
            for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures), desc="Running Prokka"):
                gid, results = future.result()
                genomic_properties[gid] = results

        return genomic_properties

    def _run_codetta_for_genome(self, gid: str, genome_path: str, codetta_dir: str) -> tuple:
        """Run Codetta for a single genome."""

        # TBD: replace with call to Codetta installed system wide
        inference_file = f'{codetta_dir}/{gid}-inference'
        cmd = f'~/git/codetta/codetta.py {genome_path} --resource_directory ~/git/codetta/resources'
        cmd += f' -s {codetta_dir}/{gid}-summary.tsv'
        cmd += f' --align_output {codetta_dir}/{gid}-align'
        cmd += f' --inference_output {inference_file}'
        os.system(cmd)

        # parse Codetta results
        codetta_result = None
        with open(inference_file) as f:
            for line in f:
                if line.startswith('TGA'):
                    tokens = line.strip().split()
                    codetta_result = (tokens[1], tokens[5])
                    break

        return gid, CodettaResult(codetta_result[0], int(codetta_result[1]))

    def run_codetta(self, genome_paths_decompressed: dict, cpus: int, out_dir: str) -> dict:
        """Determine support for UGA reassignment according to Codetta."""

        codetta_dir = os.path.join(out_dir, 'codetta')
        os.makedirs(codetta_dir, exist_ok=True)

        codetta_results = {}
        num_workers = max(1, cpus // 3)
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
            futures = [
                executor.submit(self._run_codetta_for_genome, gid, genome_paths_decompressed[gid], codetta_dir)
                for gid in genome_paths_decompressed
            ]
            
            for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures), desc="Running Codetta"):
                gid, result = future.result()
                codetta_results[gid] = result

        return codetta_results

    def write_results(self, 
                        ground_truth: dict,
                        misclassified_gids: dict, 
                        genomic_properties: dict, 
                        codetta_results: dict,
                        out_dir: str) -> None:
        """Write file with evidence for or against a UGA reassignment."""

        fout = open(os.path.join(out_dir, 'misclassification_evidence.tsv'), 'wt')
        fout.write('Genome ID\tGround truth\tgTranslate prediction\tCheckM prediction\tCodetta prediction\tCodetta columns')
        fout.write('\tgTranslate correct\tCheckM correct\tCodetta correct')
        fout.write('\tTable 11: tRNA-Trp/Gly(tca)\tTable 11: RF2')
        fout.write('\tTable 4: tRNA-Trp/Gly(tca)\tTable 4: RF2')
        fout.write('\tGTDB taxonomy\tNCBI taxonomy')
        fout.write('\tCoding density 4 or 25\tCoding density 11)')
        fout.write('\tGC content\tN50\tGenome size\tContig count\n')

        for gid in misclassified_gids:
            # determine CheckM prediction
            cd4 = float(misclassified_gids[gid].cd4)
            cd11 = float(misclassified_gids[gid].cd11)
            checkm_prediction = '11'
            if (cd4 - cd11) > 5 and cd4 > 70:
                checkm_prediction = '4 or 25'

            checkm_correct = ground_truth[gid].translation_table == checkm_prediction or (ground_truth[gid].translation_table in ['4', '25'] and checkm_prediction == '4 or 25')

            # determine Codetta prediction
            amino_acid, codetta_cols = codetta_results[gid].base, codetta_results[gid].columns
            codetta_tt = '11'
            if amino_acid == 'W':
                codetta_tt = '4'
            elif amino_acid == 'G':
                codetta_tt = '25'
            elif amino_acid != '?':
                codetta_tt = amino_acid

            if int(codetta_cols) < 34:
                codetta_correct = ground_truth[gid].translation_table == '11'
            else:
                codetta_correct = ground_truth[gid].translation_table == codetta_tt

            # write out results for genome
            fout.write(gid)
            fout.write('\t{}\t{}\t{}\t{}\t{}'.format(
                ground_truth[gid].translation_table,
                misclassified_gids[gid].translation_table,
                checkm_prediction,
                codetta_tt,
                codetta_cols
            ))

            fout.write('\t{}\t{}\t{}'.format(
                misclassified_gids[gid].translation_table == ground_truth[gid].translation_table,
                checkm_correct,
                codetta_correct
            ))

            fout.write('\t{}\t{}\t{}\t{}'.format(
                genomic_properties[gid][4].tRNA_Trp_tca,
                genomic_properties[gid][4].RF2,
                genomic_properties[gid][11].tRNA_Trp_tca,
                genomic_properties[gid][11].RF2
            ))

            fout.write('\t{}\t{}'.format(
                ground_truth[gid].gtdb_taxonomy,
                ground_truth[gid].ncbi_taxonomy
            ))

            fout.write('\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                misclassified_gids[gid].cd4,
                misclassified_gids[gid].cd11,
                misclassified_gids[gid].gc,
                misclassified_gids[gid].n50,
                misclassified_gids[gid].genome_size,
                misclassified_gids[gid].contig_count
            ))

            fout.write('\n')

        fout.close()

    def run(self, 
            ground_truth_file: str,
            gtranslate_results_file: str, 
            genome_path_file: str, 
            cpus: int,
            out_dir: str) -> None:
        """Evaluate if misclassifications are likely due to an incorrect ground truth."""

        # parse ground trurth for genomes
        self.log.info('Parsing ground truth for each genome:')
        ground_truth = self.parse_ground_truth_file(ground_truth_file)
        self.log.info(f' - identified ground truth for {len(ground_truth):,} genomes')

        # identify genomes with putative misclassification by gTranslate
        self.log.info('Identifying genomes with putative misclassification by gTranslate:')
        misclassified_gids = self.parse_gtranslate_misclassifications(gtranslate_results_file, ground_truth)
        self.log.info(f' - identified {len(misclassified_gids):,} putative misclassifications')

        # read path to putatively misclassified genomes
        self.log.info('Parsing path to genomic FASTA file of putatively misclassified genomes:')
        genome_paths = self.parse_genome_path_file(genome_path_file, misclassified_gids)
        self.log.info(f' - identified path for {len(misclassified_gids):,} genomes')
        assert len(genome_paths) == len(misclassified_gids)

        # decompress putatively misclassified genomes
        self.log.info('Decompressing putatively misclassified genomes:')
        genome_paths_decompressed = self.decompress_genomes(genome_paths, misclassified_gids, out_dir)
        self.log.info(f' - decompressed {len(genome_paths_decompressed):,} genomes')

        # identify genomic properties of UGA reassignment using Prokka
        self.log.info('Identifying genomic properties of UGA reassignment using Prokka:')
        genomic_properties = self.identify_uga_reassignment_properties(genome_paths_decompressed,
                                                                        cpus, 
                                                                        out_dir)
        self.log.info(f' - identified genomic properties for {len(genomic_properties):,} genomes')

        # determine support for UGA reassignment according to Codetta
        self.log.info('Identifying support for UGA reassignment according to Codetta:')
        codetta_results = self.run_codetta(genome_paths_decompressed, cpus, out_dir)
        self.log.info(f' - identified genomic properties for {len(genomic_properties):,} genomes')

        # write file with evidence for or against a UGA reassignment
        self.log.info('Writing file with evidence for or against a UGA reassignment.')
        self.write_results(ground_truth,
                           misclassified_gids, 
                           genomic_properties, 
                           codetta_results,
                           out_dir)

        self.log.info('Done.')


def main():
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--ground_truth_file', required=True, help='File with ground truth for each genome.')
    parser.add_argument('--gtranslate_results_file', required=True, help='File with gTranslate prediction for each genome.')
    parser.add_argument('--genome_path_file', required=True,  help='File indicating path to genomic FASTA files.')
    parser.add_argument('--out_dir', required=True, help='Output directory.')
    parser.add_argument('--cpus', type=int, default=1, help='Number of CPUs to use.')

    args = parser.parse_args()

    # setup logger
    logger_setup(args.out_dir, 
                 __prog_name__.replace('.py', '.log'), 
                 __prog_name__.replace('.py', ''),
                 __version__, 
                 False, False)

    # run program
    p = EvaluateMisclassifications()
    p.run(args.ground_truth_file,
            args.gtranslate_results_file, 
            args.genome_path_file, 
            args.cpus,
            args.out_dir)


if __name__ == '__main__':
    main()
