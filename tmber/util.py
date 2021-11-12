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
    print_log(f'Read a FASTA file:\t{path}')
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


def read_bed(path, merge=True, **kwargs):
    print_log(
        'Read a BED file{0}:\t{1}'.format(
            (' (merging intervals)' if merge else ''), path
        )
    )
    dtype = {
        'chrom': str, 'chromStart': int, 'chromEnd': int, 'name': str,
        'score': int, 'strand': str, 'thickStart': int, 'thickEnd': int,
        'itemRgb': str, 'blockCount': int, 'blockSizes': int,
        'blockStarts': int, 'ADDITIONAL': str
    }
    bed_lines = [d for d in _stream_bed_lines(path=path, **kwargs)]
    return (
        pd.DataFrame(bed_lines) if bed_lines else pd.DataFrame()
    ).pipe(
        lambda d:
        d.astype(dtype={k: v for k, v in dtype.items() if k in d.columns})
    )


def read_vcf(path, sample_name=None, min_af=None, max_af=None,
             include_filtered=False, **kwargs):
    print_log(
        'Read a VCF file ({0} filtered variants):\t{1}'.format(
            ('including' if include_filtered else 'excluding'), path
        )
    )
    dtype = {
        'CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
        'QUAL': str, 'FILTER': str, 'INFO': str, 'FORMAT': str
    }
    vcf_lines = [d for d in _stream_vcf_lines(path=path, **kwargs)]
    if not vcf_lines:
        return pd.DataFrame()
    else:
        df_vcf = pd.DataFrame(vcf_lines).pipe(
            lambda d:
            d.astype(dtype={k: v for k, v in dtype.items() if k in d.columns})
        )
        if min_af is None and max_af is None:
            return df_vcf
        else:
            if min_af is None:
                condition_str = f'AF <= {max_af}'
            elif max_af is None:
                condition_str = f'AF >= {min_af}'
            else:
                condition_str = f'{min_af} <= AF <= {max_af}'
            print_log(f'Extract variants with {condition_str}:\t{sample_name}')
            return df_vcf.assign(
                AF=lambda d: _extract_sample_af(df=d, sample_name=sample_name)
            ).pipe(
                lambda d: (d if min_af is None else d[d['AF'] >= min_af])
            ).pipe(
                lambda d: (d if max_af is None else d[d['AF'] <= max_af])
            )


def _extract_sample_af(df, sample_name):
    assert sample_name in df.columns, f'column not found: {sample_name}'
    return df[['FORMAT', sample_name]].apply(
        lambda r: dict(zip(r[0].split(':'), r[1].split(':'))).get('AF'),
        axis=1
    ).astype(float)


def _stream_vcf_lines(path, include_filtered=False, **kwargs):
    columns = None
    maxsplit = None
    for s in _open_and_stream_file(path=path, **kwargs):
        if not s.startswith('#'):
            assert (columns and maxsplit), 'columns not found'
            values = s.strip().split('\t', maxsplit=maxsplit)
            od = OrderedDict(zip(columns, values))
            if include_filtered or od['FILTER'] in {'PASS', '.'}:
                yield od
        elif s.startswith('#CHROM'):
            columns = s[1:].strip().split('\t')
            maxsplit = len(columns) - 1


def _stream_bed_lines(path, merge=True, bedtools='bedtools', **kwargs):
    columns = [
        'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
        'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes',
        'blockStarts', 'ADDITIONAL'
    ]
    maxsplit = len(columns) - 1
    iters = (
        _run_and_parse_subprocess(
            args=f'{bedtools} sort -i {path} | {bedtools} merge -i -',
            shell=True
        ) if merge else _open_and_stream_file(path=path, **kwargs)
    )
    for s in iters:
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
