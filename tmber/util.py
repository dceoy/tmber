#!/usr/bin/env python

import bz2
import gzip
import logging
import os
import subprocess
from collections import OrderedDict
from pathlib import Path
from pprint import pformat

import pandas as pd
import yaml


def fetch_executable(cmd, ignore_errors=False):
    executables = [
        cp for cp in [
            str(Path(p).joinpath(cmd))
            for p in os.environ['PATH'].split(os.pathsep)
        ] if os.access(cp, os.X_OK)
    ]
    if executables:
        return executables[0]
    elif ignore_errors:
        return None
    else:
        raise RuntimeError(f'command not found: {cmd}')


def stream_vcf_lines(path, exclude_filtered=True, bgzip='bgzip', pigz=None,
                     pbzip2=None, n_cpu=1):
    dtype = OrderedDict([
        ('CHROM', str), ('POS', int), ('ID', str), ('REF', str), ('ALT', str),
        ('QUAL', str), ('FILTER', str), ('INFO', str), ('FORMAT', str),
        ('ADDITIONAL', str)
    ])
    columns = list(dtype.keys())
    maxsplit = len(columns) - 1
    for s in _open_and_stream_file(path=path, bgzip=bgzip, pigz=pigz,
                                   pbzip2=pbzip2, n_cpu=n_cpu):
        if not s.startswith('#'):
            v = s.strip().split('\t', maxsplit=maxsplit)
            df = pd.DataFrame([v], columns=columns[:len(v)])
            if not (exclude_filtered
                    and df.iloc[0]['FILTER'] not in {'PASS', '.'}):
                yield df.astype(dtype=dtype)


def stream_bed_lines(path, bgzip='bgzip', pigz=None, pbzip2=None, n_cpu=1):
    dtype = OrderedDict([
        ('chrom', str), ('chromStart', int), ('chromEnd', int), ('name', str),
        ('score', int), ('strand', str), ('thickStart', int),
        ('thickEnd', int), ('itemRgb', str), ('blockCount', int),
        ('blockSizes', int), ('blockStarts', int), ('ADDITIONAL', str)
    ])
    columns = list(dtype.keys())
    maxsplit = len(columns) - 1
    for s in _open_and_stream_file(path=path, bgzip=bgzip, pigz=pigz,
                                   pbzip2=pbzip2, n_cpu=n_cpu):
        if not s.startswith(('browser', 'track')):
            v = s.strip().split('\t', maxsplit=maxsplit)
            yield pd.DataFrame(
                [v], columns=columns[:len(v)]
            ).astype(dtype=dtype)


def _open_and_stream_file(path, bgzip='bgzip', pigz=None, pbzip2=None,
                          n_cpu=1):
    tsv_path = str(path)
    if tsv_path.endswith(('.gz', '.bgz')):
        if bgzip:
            f = _run_and_parse_subprocess(
                args=[bgzip, '-@', str(n_cpu), '-dc', tsv_path]
            )
        elif pigz:
            f = _run_and_parse_subprocess(
                args=[pigz, '-p', str(n_cpu), '-dc', tsv_path]
            )
        else:
            f = gzip.open(tsv_path, mode='rt')
    elif tsv_path.endswith('.bz2'):
        if pbzip2:
            f = _run_and_parse_subprocess(
                args=[pbzip2, f'-p{n_cpu}', '-dc', tsv_path]
            )
        else:
            f = bz2.open(tsv_path, mode='rt')
    else:
        f = open(tsv_path, mode='r')
    for s in f:
        yield s


def _run_and_parse_subprocess(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, **kwargs):
    logger = logging.getLogger(__name__)
    logger.debug(f'args: {args}')
    with subprocess.Popen(args=args, stdout=stdout, stderr=stderr,
                          **kwargs) as p:
        for line in p.stdout:
            yield line.decode('utf-8')
        o, e = p.communicate()
        if p.returncode == 0:
            pass
        else:
            logger.error(
                f'STDERR from subprocess `{p.args}`:{os.linesep}'
                + e.decode('utf-8')
            )
            raise subprocess.CalledProcessError(
                returncode=p.returncode, cmd=p.args, output=o, stderr=e
            )


def print_log(message):
    logger = logging.getLogger(__name__)
    logger.debug(message)
    print(f'>>\t{message}', flush=True)


def read_yml(path):
    logger = logging.getLogger(__name__)
    with open(str(path), 'r') as f:
        d = yaml.load(f, Loader=yaml.FullLoader)
    logger.debug('YAML data:' + os.linesep + pformat(d))
    return d


def print_yml(data):
    print(yaml.dump(data))
