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
from Bio import SeqIO


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


def read_fasta_and_generate_seq(path):
    print_log(f'Read a FASTA file: {path}')
    if path.endswith('.gz'):
        f = gzip.open(path, 'rt')
    elif path.endswith('.bz2'):
        f = bz2.open(path, 'rt')
    else:
        f = open(path, 'r')
    for s in SeqIO.parse(f, 'fasta'):
        yield s.id, s.seq
    f.close()


def print_yml(data):
    logger = logging.getLogger(__name__)
    logger.debug(data)
    print(yaml.dump(data))


def read_bed(path, columns=None, **kwargs):
    print_log(f'Read a BED file: {path}')
    dtype = {
        'chrom': str, 'chromStart': int, 'chromEnd': int, 'name': str,
        'score': int, 'strand': str, 'thickStart': int, 'thickEnd': int,
        'itemRgb': str, 'blockCount': int, 'blockSizes': int,
        'blockStarts': int, 'ADDITIONAL': str
    }
    bed_lines = (
        [
            OrderedDict([(k, v) for k, v in d.items() if k in columns])
            for d in _stream_bed_lines(path=path, **kwargs)
        ] if columns else [d for d in _stream_bed_lines(path=path, **kwargs)]
    )
    return (
        pd.DataFrame(bed_lines) if bed_lines else pd.DataFrame(columns=columns)
    ).pipe(
        lambda d:
        d.astype(dtype={k: v for k, v in dtype.items() if k in d.columns})
    )


def read_vcf(path, columns=None, **kwargs):
    print_log(f'Read a VCF file: {path}')
    dtype = {
        'CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
        'QUAL': str, 'FILTER': str, 'INFO': str, 'FORMAT': str,
        'ADDITIONAL': str
    }
    vcf_lines = (
        [
            OrderedDict([(k, v) for k, v in d.items() if k in columns])
            for d in _stream_vcf_lines(path=path, **kwargs)
        ] if columns else [d for d in _stream_vcf_lines(path=path, **kwargs)]
    )
    return (
        pd.DataFrame(vcf_lines) if vcf_lines else pd.DataFrame(columns=columns)
    ).pipe(
        lambda d:
        d.astype(dtype={k: v for k, v in dtype.items() if k in d.columns})
    )


def _stream_vcf_lines(path, include_filtered=False, bgzip='bgzip', pigz=None,
                      pbzip2=None, n_cpu=1):
    columns = [
        'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
        'ADDITIONAL'
    ]
    maxsplit = len(columns) - 1
    for s in _open_and_stream_file(path=path, bgzip=bgzip, pigz=pigz,
                                   pbzip2=pbzip2, n_cpu=n_cpu):
        if not s.startswith('#'):
            values = s.strip().split('\t', maxsplit=maxsplit)
            od = OrderedDict(zip(columns[:len(values)], values))
            if include_filtered or od['FILTER'] in {'PASS', '.'}:
                yield od


def _stream_bed_lines(path, bgzip='bgzip', pigz=None, pbzip2=None, n_cpu=1):
    columns = [
        'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
        'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes',
        'blockStarts', 'ADDITIONAL'
    ]
    maxsplit = len(columns) - 1
    for s in _open_and_stream_file(path=path, bgzip=bgzip, pigz=pigz,
                                   pbzip2=pbzip2, n_cpu=n_cpu):
        if not s.startswith(('browser', 'track')):
            values = s.strip().split('\t', maxsplit=maxsplit)
            yield OrderedDict(zip(columns[:len(values)], values))


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
