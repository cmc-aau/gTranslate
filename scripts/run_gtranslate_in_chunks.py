#! /usr/bin/env python3

# Run gTranslate across genomes in chunks.

__prog_name__ = 'run_gtranslate_in_chunks.py'
__prog_desc__ = 'Run gTranslate across genomes in chunks.'

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
import subprocess

from gtranslate.biolib_lite.logger import logger_setup
from gtranslate.biolib_lite.execute import check_dependencies


class gTranslateInChunks(object):
    """Run gTranslate across genomes in chunks."""

    def __init__(self):
        """Initialize."""

        check_dependencies(['prodigal'])

        self.log = logging.getLogger('timestamp')

    def run(self, 
            genome_path_file: str, 
            custom_model_path: str,
            chunk_size: int,
            cpus: int,
            out_dir: str) -> None:
        """Run gTranslate across genomes in chunks."""

        # read path to genomic FASTA files
        self.log.info('Reading path to genomic FASTA files:')
        genome_paths = []
        with open(genome_path_file) as f:
            for line in f:
                tokens = line.strip().split('\t')
                genome_paths.append(tokens[0])

        self.log.info(f' - identified {len(genome_paths):} genomes')

        # process genomes in chunks
        num_chunks = (len(genome_paths) + chunk_size - 1) // chunk_size
        self.log.info(f'Processing genomes in {num_chunks:,} chunks:')
        for chunk_num in range(num_chunks):
            chunk_start = chunk_num * chunk_size
            chunk_end = min((chunk_num + 1) * chunk_size, len(genome_paths))

            genome_chunk_dir = os.path.join(out_dir, f'genome_chunk{chunk_num}')
            if os.path.exists(genome_chunk_dir):
                self.log.warning(f'Skipping: {genome_chunk_dir}')
                continue

            os.makedirs(genome_chunk_dir)
            for genome_path in genome_paths[chunk_start:chunk_end]:
                os.symlink(genome_path, os.path.join(genome_chunk_dir, os.path.basename(genome_path)))

            gtranslate_out_dir = os.path.join(out_dir, f'gtranslate{chunk_num}')

            if False:
                cmd = [
                    'gtranslate', 'detect_table',
                    '--genome_dir', genome_chunk_dir,
                    '--out_dir', gtranslate_out_dir,
                    '-x', 'gz',
                    '--custom_model_path', custom_model_path,
                    '--cpus', str(cpus)
                ]
                subprocess.run(cmd, check=True)
            else:
                # useful for local execution during development
                cmd = 'python3 ~/git/gtranslate/gtranslate detect_table'
                cmd += f' --genome_dir {genome_chunk_dir}'
                cmd += f' --out_dir {gtranslate_out_dir}'
                cmd += ' -x gz'
                cmd += f' --cpus {cpus}'
                cmd += f' --custom_model_path {custom_model_path}'
                os.system(cmd)

            self.log.info(f' - processed chunk {chunk_num+1}')

        # combine gTranslate output for chunks into single file
        self.log.info('Combining results into single file.')
        fout = open(os.path.join(out_dir, 'gtranslate.translation_table_summary.tsv'), 'w')
        write_header = True
        for chunk_num in range(num_chunks):
            gtranslate_out_dir = os.path.join(out_dir, f'gtranslate{chunk_num}')
            if not os.path.exists(gtranslate_out_dir):
                self.log.error(f'Missing expected results directory: {gtranslate_out_dir}.')
                continue
            
            gtranslate_results_file = os.path.join(gtranslate_out_dir, 'gtranslate.translation_table_summary.tsv')
            with open(gtranslate_results_file, 'r') as f:
                header = f.readline().strip().split('\t')
                if write_header:
                    write_header = False
                    fout.write('\t'.join(header) + '\n')

                for line in f:
                    tokens = line.strip().split('\t')
                    fout.write(line)

        fout.close()

        self.log.info('Done.')


def main():
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--genome_path_file', required=True,  help='File indicating path to genomic FASTA files.')
    parser.add_argument('--custom_model_path', required=True, help='Path to directory containing model files.')
    parser.add_argument('--out_dir', required=True, help='Output directory.')
    parser.add_argument('--chunk_size', default=10000, type=int, help='Chunk size.')
    parser.add_argument('--cpus', default=1, type=int, help='Number of CPUs to use.')

    args = parser.parse_args()

    # setup logger
    logger_setup(args.out_dir, 
                 __prog_name__.replace('.py', '.log'), 
                 __prog_name__.replace('.py', ''),
                  __version__, 
                  False, False)

    # run program
    p = gTranslateInChunks()
    p.run(args.genome_path_file, 
            args.custom_model_path,
            args.chunk_size,
            args.cpus,
            args.out_dir)

    
if __name__ == '__main__':
    main()
