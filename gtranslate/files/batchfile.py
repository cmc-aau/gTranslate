import logging

from gtranslate.exceptions import GTranslateExit


class Batchfile(object):

    def __init__(self, path: str):
        self.path = path
        self.genome_path= self.read(path)

    @staticmethod
    def read(path):
        logger = logging.getLogger('timestamp')
        genomes = dict()
        seen_paths = set()
        warnings = list()
        with open(path) as fh:
            for line_no, line in enumerate(fh):
                line_split = line.strip().split("\t")
                if line_split[0] == '':
                    continue  # blank line

                if len(line_split) not in {2}:
                    raise GTranslateExit('Batch file must contain 2 '
                                     'columns (genome_path, genome_id).')

                if len(line_split) == 2:
                    genome_file, genome_id = line_split


                if genome_file is None or genome_file == '':
                    warnings.append(f'Missing genome path on line {line_no + 1}.')
                elif genome_id is None or genome_id == '':
                    warnings.append(f'Missing genome ID on line {line_no + 1}.')
                elif genome_id in genomes:
                    warnings.append(f'Genome ID {genome_id} appears multiple times.')
                if genome_file in seen_paths:
                    logger.warning(f'Genome file appears multiple times: {genome_file}')

                # All good, record the value.
                genomes[genome_id] = genome_file
                seen_paths.add(genome_file)

        # Check if any warnings were raised.
        if len(warnings) > 0:
            warning_str = '\n'.join(warnings)
            raise GTranslateExit(f'Please check the format of your batchfile, '
                             f'the following errors were found: {warning_str}')

        return genomes
