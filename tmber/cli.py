#!/usr/bin/env python
"""
Tumor Mutational Burden Analyzer

Usage:
    tmber -h|--help
    tmber --version
    tmber tmb [--debug|--info] [--cpus=<int>] [--include-filtered]
        [--dest-dir=<path>] <bed_path> <vcf_path>...
    tmber fa2bed [--debug|--info] [--cpus=<int>] [--uppercase]
        [--human-autosome] <fa_path> <bed_path>

Commands:
    tmb                     Calculate variant counts on a BED resion
    fa2bed                  Create an accessible BED file from a FASTA file
                            (excluding Ns in sequences)

Options:
    -h, --help              Print help and exit
    --version               Print version and exit
    --debug, --info         Execute a command with debug|info messages
    --cpus=<int>            Limit CPU cores used
    --include-filtered      Include filtered variants (`PASS` or `.`)
    --dest-dir=<path>       Specify a path to an output TSV file [default: .]
    --uppercase             Extract only uppercases in sequences
    --human-autosome        Extract only human autosomes (chr1-22)

Args:
    <bed_path>              Path to a BED file for TMB
    <vcf_path>              Path to a VCF file
    <fa_path>               Path to a genome FASTA file
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
    if args['fa2bed']:
        create_bed_from_fa(
            fa_path=args['<fa_path>'], bed_path=args['<bed_path>'],
            bgzip=fetch_executable('bgzip'), uppercase=args['--uppercase'],
            human_autosome=args['--human-autosome'],
            n_cpu=int(args['--cpus'] or cpu_count())
        )
    if args['tmb']:
        calculate_tmb(
            vcf_paths=args['<vcf_path>'], bed_path=args['<bed_path>'],
            dest_dir_path=args['--dest-dir'], bgzip=fetch_executable('bgzip'),
            include_filtered=args['--include-filtered'],
            n_cpu=int(args['--cpus'] or cpu_count())
        )
