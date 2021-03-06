#!/usr/bin/env python
"""
Tumor Mutational Burden Analyzer

Usage:
    tmber -h|--help
    tmber --version
    tmber bed [--debug|--info] [--cpus=<int>] [--human-autosome]
        [--target-letters=<str>] [--dest-dir=<path>] <fa_path>
    tmber tmb [--debug|--info] [--cpus=<int>] [--include-filtered]
        [--sample=<name>] [--min-af=<float>] [--max-af=<float>]
        [--dest-dir=<path>] <vcf_path> <bed_path>...

Commands:
    bed                 Identify regions consisting of target letters in FASTA
    tmb                 Tally variants and calculate TMB on a BED resion

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
    --sample=<name>     Specify a sample column including AF in a VCF file
    --min-af=<float>    Specify the min AF limit
    --max-af=<float>    Specify the max AF limit

Args:
    <fa_path>           Path to a genome FASTA file
    <vcf_path>          Path to a VCF file
    <bed_path>          Path to a BED file for TMB
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
    logger.info(f'n_cpu: {n_cpu}')
    if args['bed']:
        create_bed_from_fa(
            fa_path=args['<fa_path>'], dest_dir_path=args['--dest-dir'],
            bgzip=fetch_executable('bgzip'),
            human_autosome=args['--human-autosome'],
            target_letters=args['--target-letters'], n_cpu=n_cpu
        )
    if args['tmb']:
        if args['--min-af'] or args['--max-af']:
            assert bool(args['--sample']), '--sample requred for AF limits'
        calculate_tmb(
            vcf_path=args['<vcf_path>'], bed_paths=args['<bed_path>'],
            dest_dir_path=args['--dest-dir'],
            bedtools=fetch_executable('bedtools'),
            bgzip=fetch_executable('bgzip'),
            include_filtered=args['--include-filtered'],
            sample_name=args['--sample'],
            min_af=(float(args['--min-af']) if args['--min-af'] else None),
            max_af=(float(args['--max-af']) if args['--max-af'] else None),
            n_cpu=n_cpu
        )
