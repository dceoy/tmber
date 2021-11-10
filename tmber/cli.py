#!/usr/bin/env python
"""
Tumor Mutational Burden Analyzer

Usage:
    tmber -h|--help
    tmber --version
    tmber tmb [--debug|--info] [--cpus=<int>] [--dest-dir=<path>] <bed_path>
        <vcf_path>...

Options:
    -h, --help          Print help and exit
    --version           Print version and exit
    --debug, --info     Execute a command with debug|info messages
    --cpus=<int>        Limit CPU cores used
    --dest-dir=<path>   Specify a destination directory path [default: .]

Args:
    <bed_path>          Path to a BED file for TMB
    <vcf_path>          Path to a VCF file
"""

import logging
import os

from docopt import docopt
from psutil import cpu_count

from . import __version__
from .tmb import calculate_tmb
from .util import fetch_executable


def main():
    args = docopt(__doc__, version=__version__)
    _set_log_config(debug=args['--debug'], info=args['--info'])
    logger = logging.getLogger(__name__)
    logger.debug(f'args:{os.linesep}{args}')
    print(args)
    if args['tmb']:
        calculate_tmb(
            vcf_paths=args['<vcf_path>'], bed_path=args['<bed_path>'],
            dest_dir_path=args['--dest-dir'], bgzip=fetch_executable('bgzip'),
            n_cpu=int(args['--cpus'] or cpu_count()),
            executable=fetch_executable('bash')
        )


def _set_log_config(debug=None, info=None):
    if debug:
        lv = logging.DEBUG
    elif info:
        lv = logging.INFO
    else:
        lv = logging.WARNING
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S', level=lv
    )
