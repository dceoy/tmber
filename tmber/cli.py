#!/usr/bin/env python
"""
Tumor Mutational Burden Analyzer

Usage:
    tmber -h|--help
    tmber --version
    tmber bed [--debug|--info] [--cpus=<int>] [--human-autosome]
        [--target-letters=<str>] [--dest-dir=<path>] <fa_path>
    tmber tmb [--debug|--info] [--cpus=<int>] [--include-filtered]
        [--dest-dir=<path>] <bed_path> <vcf_path>...

Commands:
    bed                 Identify regions consisting of target letters in FASTA
    tmb                 Calculate variant counts on a BED resion

Options:
    -h, --help          Print help and exit
    --version           Print version and exit
    --debug, --info     Execute a command with debug|info messages
    --cpus=<int>        Limit CPU cores to use
    --human-autosome    Extract only human autosomes (chr1-22)
    --target-letters=<str>
                        Specify nucleic acid codes to include [default: ACGT]
    --dest-dir=<path>   Specify a path to an output TSV file [default: .]
    --include-filtered  Include filtered variants (`PASS` or `.`)

Args:
    <fa_path>           Path to a genome FASTA file
    <bed_path>          Path to a BED file for TMB
    <vcf_path>          Path to a VCF file
"""

import logging
import os

from docopt import docopt
from psutil import cpu_count

from . import __version__
from .bed import create_bed_from_fa
from .tmb import calculate_tmb
from .util import fetch_executable


def main():
    args = docopt(__doc__, version=__version__)
    if args['--debug']:
        lv = logging.DEBUG
    elif args['--info']:
        lv = logging.INFO
    else:
        lv = logging.WARNING
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S', level=lv
    )
    logger = logging.getLogger(__name__)
    logger.debug(f'args:{os.linesep}{args}')
    n_cpu = int(args['--cpus'] or cpu_count())
    if args['bed']:
        create_bed_from_fa(
            fa_path=args['<fa_path>'], dest_dir_path=args['--dest-dir'],
            bgzip=fetch_executable('bgzip'),
            human_autosome=args['--human-autosome'],
            target_letters=args['--target-letters'],
            n_cpu=n_cpu
        )
    if args['tmb']:
        calculate_tmb(
            vcf_paths=args['<vcf_path>'], bed_path=args['<bed_path>'],
            dest_dir_path=args['--dest-dir'],
            bedtools=fetch_executable('bedtools'),
            bgzip=fetch_executable('bgzip'),
            include_filtered=args['--include-filtered'], n_cpu=n_cpu
        )
