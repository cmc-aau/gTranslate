import hashlib
import logging
import shutil

import math
import os
import random
import re
import time
from itertools import islice

from tqdm import tqdm


from gtranslate.config.output import CHECKSUM_SUFFIX, DIR_IDENTIFY_INTERMEDIATE
from gtranslate.biolib_lite.seq_io import read_fasta

order_rank = ["d__", "p__", "c__", "o__", 'f__', 'g__', 's__']

##################################################
############MISC UTILITIES########################
##################################################

RE_CANONICAL = re.compile(r'^(?:GB_)?(?:RS_)?(?:GCF_)?(?:GCA_)?(\d{9})\.\d')


def canonical_gid(gid: str) -> str:
    """Get canonical form of NCBI genome accession.

    Example:
        G005435135 -> G005435135
        GCF_005435135.1 -> G005435135
        GCF_005435135.1_ASM543513v1_genomic -> G005435135
        RS_GCF_005435135.1 -> G005435135
        GB_GCA_005435135.1 ->

    :param gid: Genome accesion to conver to canonical form.
    :return: Canonical form of accession.
    """

    match = RE_CANONICAL.match(gid)
    if match:
        return f'G{match[1]}'
    else:
        return gid


def get_genomes_size(genome_path):
    """Returns the size of a specific genome file."""
    seqs = read_fasta(genome_path)
    return sum([len(seq) for seq in seqs.values()])


def splitchunks(d, n):
    chunksize = int(math.ceil(len(d) / float(n)))
    it = iter(d)
    for _ in range(0, len(d), chunksize):
        yield {k: d[k] for k in islice(it, chunksize)}

def splitchunks_list(l, n):
    """Yield successive n-sized chunks from l."""
    chunksize = int(math.ceil(len(l) / float(n)))
    for i in range(0, len(l), chunksize):
        yield l[i:i + chunksize]


def generateTempTableName():
    rng = random.SystemRandom()
    suffix = ''
    for _ in range(0, 10):
        suffix += rng.choice(
            'abcefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')
    return "TEMP" + suffix + str(int(time.time()))


def merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z

def confirm(msg):
    raw = input(msg + " (y/N): ")
    if raw.upper() in ["Y", "YES"]:
        return True
    return False

def sha256(input_file):
    """Determine SHA256 hash for file.

    Parameters
    ----------
    input_file : str
        Name of file.
    Returns
    -------
    str
        SHA256 hash.
    """
    block_size = 65536
    hasher = hashlib.sha256()
    with open(input_file, 'rb') as afile:
        buf = afile.read(block_size)
        while len(buf) > 0:
            hasher.update(buf)
            buf = afile.read(block_size)
    return hasher.hexdigest()


def file_has_checksum(file_path, checksum_suffix=CHECKSUM_SUFFIX):
    """Check that the file contents match the checksum.

    Parameters
    ----------
    file_path : str
        Name of the file to check.
    checksum_suffix : str
        Suffix used to denote the file checksum.

    Returns
    -------
    bool
        True if the file has a checksum and it matches the original contents,
        False otherwise.

    """
    check_path = file_path + checksum_suffix
    if os.path.isfile(file_path) and os.path.isfile(check_path):
        with open(check_path, 'r') as check_f:
            return sha256(file_path) == check_f.read()
    return False


def symlink_f(src, dst, force=True):
    """Create a symbolic link pointing to src named dst.

    Parameters
    ----------
    src : str
        The source file.
    dst : str
        The destination file.
    force : bool
        Overwrite any file found with the same name as dst.

    """
    if force and os.path.isfile(dst):
        os.remove(dst)
    os.symlink(src, dst)


def remove_intermediate_files(output_dir):
    """Remove intermediate files.

    Parameters
    ----------
    output_dir : str
        The path to the output directory.
    wf_name : str
        The name of the workflow to delete intermediate files.
    """
    #Remove predict step intermediate files
    intermediate_identify = os.path.join(output_dir, DIR_IDENTIFY_INTERMEDIATE)
    if os.path.exists(intermediate_identify) and os.path.isdir(intermediate_identify):
        shutil.rmtree(intermediate_identify)



class tqdm_log(object):
    """A basic wrapper for the tqdm progress bar. Automatically reports the
    runtime statistics after exit.
    """

    def __init__(self, iterable=None, **kwargs):
        # Setup reporting information.
        self.logger = logging.getLogger('timestamp')
        self.start_ts = None

        # Set default parameters.
        default = {'leave': False,
                   'smoothing': 0.1,
                   'bar_format': '==> Processed {n_fmt}/{total_fmt} {unit}s '
                                 '({percentage:.0f}%) |{bar:15}| [{rate_fmt}, ETA {remaining}] {postfix}'}
        merged = {**default, **kwargs}
        self.args = merged

        # Instantiate tqdm
        self.tqdm = tqdm(iterable, **merged)

    def log(self):
        try:
            # Collect values.
            delta = time.time() - self.start_ts
            n = self.tqdm.n
            unit = self.args.get('unit', 'item')

            # Determine the scale.
            if delta > 60:
                time_unit = 'minute'
                scale = 1 / 60
            elif delta > 60 * 60:
                time_unit = 'hour'
                scale = 1 / (60 * 60)
            elif delta > 60 * 60 * 24:
                time_unit = 'day'
                scale = 1 / (60 * 60 * 24)
            else:
                time_unit = 'second'
                scale = 1

            # Scale units.
            value = delta * scale
            per = n / value

            # Determine what scale to use for the output.
            if per > 1:
                per_msg = f'{per:,.2f} {unit}s/{time_unit}'
            else:
                per_msg = f'{1 / per:,.2f} {time_unit}s/{unit}'

            # Output the message.
            s = 's' if n > 1 else ''
            msg = f'Completed {n:,} {unit}{s} in {value:,.2f} {time_unit}s ({per_msg}).'
            self.logger.info(msg)
        except Exception:
            pass


    def __enter__(self):
        self.start_ts = time.time()
        return self.tqdm

    def __iter__(self):
        try:
            self.start_ts = time.time()
            for item in self.tqdm:
                yield item
            self.log()
        finally:
            self.tqdm.close()

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.log()
        self.tqdm.close()

    def __del__(self):
        self.tqdm.close()
