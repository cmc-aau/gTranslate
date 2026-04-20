#! /usr/bin/env python3

# Creating files required to train gTranslate on a new GTDB release.
#
# Training should be performed on:
#  - all GTDB genomes passing standard QC
#  - additional genomes with reassigned stop codons that pass a less stringent set of QC criteria

__prog_name__ = 'gtdb_training_set.py'
__prog_desc__ = 'Creating files required to train gTranslate on a new GTDB release.'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2026'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.1.0'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'


import os
import logging
import argparse
import gzip
import random
from collections import defaultdict
from dataclasses import dataclass

from gtranslate.biolib_lite.logger import logger_setup
from gtranslate.tools import canonical_gid


@dataclass
class GenomeData:
    gtdb_taxonomy: str
    ncbi_taxonomy: str


class GtdbTrainingSet(object):
    """Creating files required to train gTranslate on a new GTDB release."""

    def __init__(self):
        """Initialize."""

        # QC criteria for above taxa
        self.MIN_COMP = 50
        self.MAX_CONT = 10
        self.MIN_QUAL = 50

        # genomes to retain without needing to pass QC
        self.QC_EXCEPTIONS = set(['G015697925']) # s__Pinguicoccus supinus representative genome

        # maximum number of genomes to select per species
        self.MAX_PER_SPECIES = 100

        self.log = logging.getLogger('timestamp')

    def parse_metadata_file(self, gtdb_metadata_file: str) -> tuple[dict, dict]:
        """Parse manual ground truth file."""

        genome_data = {}
        gtdb_sp_rid = {}
        open_file = gzip.open if gtdb_metadata_file.endswith('.gz') else open 
        with open_file(gtdb_metadata_file, 'rt') as f:
            header = f.readline().strip().split('\t')

            gid_idx = header.index("accession")
            gtdb_taxonomy_idx = header.index('gtdb_taxonomy')
            ncbi_taxonomy_idx = header.index('ncbi_taxonomy')

            gtdb_rep_idx = header.index('gtdb_representative')

            for line in f:
                tokens = line.strip().split('\t')
                gid = canonical_gid(tokens[gid_idx])
                genome_data[gid] = GenomeData(tokens[gtdb_taxonomy_idx], tokens[ncbi_taxonomy_idx])

                is_rep = tokens[gtdb_rep_idx].lower().startswith('t')
                if is_rep:
                    gtdb_sp = [t.strip() for t in tokens[gtdb_taxonomy_idx].split(';')][-1]
                    gtdb_sp_rid[gtdb_sp] = gid

        return genome_data, gtdb_sp_rid

    def parse_taxonomy_file(self, taxonomy_file: str) -> dict:
        """Parse taxonomy file."""

        taxonomy = {}
        open_file = gzip.open if taxonomy_file.endswith('.gz') else open 
        with open_file(taxonomy_file, 'rt') as f:
            for line in f:
                tokens = line.strip().split('\t')
                gid = canonical_gid(tokens[0])
                taxonomy[gid] = tokens[1]

        return taxonomy

    def parse_qc_failed_file(self, gtdb_qc_failed_file: str, ncbi_taxonomy: dict) -> tuple[dict, dict]:
        """Parse genomes with reassigned stop codon that failed the GTDB QC, but should still be included."""

        genome_data = {}
        open_file = gzip.open if gtdb_qc_failed_file.endswith('.gz') else open 
        with open_file(gtdb_qc_failed_file, 'rt') as f:
            header = f.readline().strip().split('\t')

            gid_idx = header.index("Accession")
            gtdb_taxonomy_idx = header.index('GTDB taxonomy')
            ncbi_sp_idx = header.index('NCBI species')

            cm2_comp_idx = header.index("CM2 completeness (%)")
            cm2_cont_idx = header.index("CM2 contamination (%)")
            cm2_qual_idx = header.index("CM2 quality")

            for line in f:
                tokens = line.strip().split('\t')
                
                gid = canonical_gid(tokens[gid_idx])

                comp = float(tokens[cm2_comp_idx])
                cont = float(tokens[cm2_cont_idx])
                qual = float(tokens[cm2_qual_idx])

                pass_qc = comp >= self.MIN_COMP and cont <= self.MAX_CONT and qual >= self.MIN_QUAL

                if gid not in self.QC_EXCEPTIONS and not pass_qc:
                    # genome failed even more lenient QC
                    continue

                gtdb_taxonomy = tokens[gtdb_taxonomy_idx]
                genome_data[gid] = GenomeData(gtdb_taxonomy, ncbi_taxonomy[gid])

        return genome_data
        
    def parse_genome_dir_file(self, genome_dir_file: str, training_gids: dict) -> dict:
        """Parse path to genomic FASTA files for training genomes."""

        genome_paths = {}
        open_file = gzip.open if genome_dir_file.endswith('.gz') else open 
        with open_file(genome_dir_file, 'rt') as f:
            for line in f:
                tokens = line.strip().split('\t')

                gid = canonical_gid(tokens[0])
                if gid not in training_gids:
                    continue

                genome_dir = tokens[1]
                accn = genome_dir.split('/')[-1]
                genome_path = os.path.join(genome_dir, f'{accn}_genomic.fna.gz')

                genome_paths[gid] = genome_path

        return genome_paths

    def run(self, gtdb_bac_metadata_file: str, 
            gtdb_ar_metadata_file: str, 
            gtdb_qc_failed: str, 
            ncbi_taxonomy_file: str,
            genome_dir_file: str, 
            out_dir: str) -> None:
        """Creating files required to train gTranslate on a new GTDB release."""

        # parse genomes from GTDB metadata files
        self.log.info('Parsing GTDB metadata files:')
        gtdb_bac_gids, gtdb_bac_sp_rid = self.parse_metadata_file(gtdb_bac_metadata_file)
        gtdb_ar_gids, gtdb_ar_sp_rid = self.parse_metadata_file(gtdb_ar_metadata_file)
        gtdb_gids = {**gtdb_bac_gids, **gtdb_ar_gids}
        gtdb_sp_rid = {**gtdb_bac_sp_rid, **gtdb_ar_sp_rid}
        self.log.info(f' - identified {len(gtdb_bac_gids):,} bacterial genomes from {len(gtdb_bac_sp_rid):,} species')
        self.log.info(f' - identified {len(gtdb_ar_gids):,} archaeal genomes from {len(gtdb_ar_sp_rid):,} species')
        self.log.info(f' - identified {len(gtdb_gids):,} total genomes from {len(gtdb_sp_rid):,} species')

        # subsample to maximum number of genomes per species
        self.log.info(f'Sampling to {self.MAX_PER_SPECIES} genomes per species:')

        sp_gids = defaultdict(set)
        for gid, genome_data in gtdb_gids.items():
            gtdb_sp = [t.strip() for t in genome_data.gtdb_taxonomy.split(';')][-1]
            sp_gids[gtdb_sp].add(gid)
        
        gtdb_gids_sampled = {}
        for sp, gids in sp_gids.items():
            if len(gids) > self.MAX_PER_SPECIES:
                rid = gtdb_sp_rid[sp]
                gids.remove(rid)
                gids = [rid] + random.sample(list(gids), self.MAX_PER_SPECIES - 1)

            for gid in gids:
                gtdb_gids_sampled[gid] = gtdb_gids[gid]

        self.log.info(f' - retained {len(gtdb_gids_sampled):,} genomes for training')

        # read NCBI taxonomy for genomes
        self.log.info('Parsing NCBI taxonomy file:')
        ncbi_taxonomy = self.parse_taxonomy_file(ncbi_taxonomy_file)
        self.log.info(f' - identified taxonomy for {len(ncbi_taxonomy):,} genomes')

        # parse genomes from taxa with reassigned stop codon that failed the GTDB QC,
        # but should still be included
        self.log.info('Parsing genomes that failed GTDB QC:')
        gtdb_reassigned_gids = self.parse_qc_failed_file(gtdb_qc_failed, ncbi_taxonomy)
        self.log.info(f' - identified {len(gtdb_reassigned_gids):,} genomes to retain for training')

        # combine all training genomes
        self.log.info('Combining all genomes to use for training:')
        training_gids = {**gtdb_gids_sampled, **gtdb_reassigned_gids}
        self.log.info(f' - identified {len(training_gids):,} total genomes for training')

        # sanity check that all genomes that are a QC expection are accounted for
        for gid in self.QC_EXCEPTIONS:
            assert gid in training_gids

        # get final count of reassigned species
        self.log.info('Tabulating number of genomes in reassigned species:')
        reassigned_sp_count = defaultdict(int)
        for gid, genome_data in training_gids.items():
            ncbi_sp = [t.strip() for t in genome_data.ncbi_taxonomy.split(';')][-1]
            gtdb_sp = [t.strip() for t in genome_data.gtdb_taxonomy.split(';')][-1]

            if ncbi_sp in self.NCBI_REASSIGNED_SP:
                reassigned_sp_count[ncbi_sp] += 1
            elif gtdb_sp in self.GTDB_REASSIGNED_SP:
                reassigned_sp_count[gtdb_sp] += 1

        for sp, count in sorted(reassigned_sp_count.items()):
            self.log.info(f' - {sp}: {count}')
                
        # read path to genomic FASTA files
        self.log.info('Parsing path to genomic FASTA files for training genomes:')
        genome_paths = self.parse_genome_dir_file(genome_dir_file, training_gids)
        self.log.info(f' - identified path for {len(training_gids):,} genomes')
        assert len(genome_paths) == len(training_gids)

        # create taxonomy file for all training genomes
        self.log.info('Creating taxonomy file.')
        fout = open(os.path.join(out_dir, 'training_taxonomy.tsv'), 'wt')
        fout.write('Genome ID\tGTDB taxonomy\tNCBI taxonomy\n')
        for gid, genome_data in training_gids.items():
            fout.write('{}\t{}\t{}\n'.format(
                gid, 
                genome_data.gtdb_taxonomy, 
                genome_data.ncbi_taxonomy))
        fout.close()

        # create file with path to genomic FASTA file for all training genomes
        self.log.info('Creating genome path file.')
        fout = open(os.path.join(out_dir, 'training_genome_path.tsv'), 'wt')
        for gid in training_gids:
            fout.write(f'{genome_paths[gid]}\t{gid}\n')
        fout.close()

        self.log.info('Done.')


def main():
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--gtdb_bac_metadata_file', required=True, help='GTDB metadata file for bacterial genomes.')
    parser.add_argument('--gtdb_ar_metadata_file', required=True,  help='GTDB metadata file for bacterial genomes.')
    parser.add_argument('--gtdb_qc_failed', required=True, help='File indicating genomes that failed standard GTDB QC criteria.')
    parser.add_argument('--ncbi_taxonomy_file', required=True,  help='File indicating standardized NCBI taxonomic classification of each genome.')
    parser.add_argument('--genome_dir_file', required=True,  help='File indicating path to data files for genome assemblies.')
    parser.add_argument('--out_dir', help='Output directory.')

    args = parser.parse_args()

    # setup logger
    logger_setup(args.out_dir, "gtdb_training_set.log", "gTranslate", __version__, False, False)

    # run program
    p = GtdbTrainingSet()
    p.run(args.gtdb_bac_metadata_file, 
          args.gtdb_ar_metadata_file, 
          args.gtdb_qc_failed,
          args.ncbi_taxonomy_file,
          args.genome_dir_file,
          args.out_dir)


if __name__ == '__main__':
    main()
