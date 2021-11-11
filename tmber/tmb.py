#!/usr/bin/env python

import logging
import os
import re
from pathlib import Path

import numpy as np

from .util import print_yml, read_bed, read_vcf


def calculate_tmb(vcf_paths, bed_path, dest_dir_path='.', bgzip='bgzip',
                  include_filtered=False, n_cpu=1):
    logger = logging.getLogger(__name__)
    vcfs = [Path(p).resolve() for p in vcf_paths]
    bed = Path(bed_path).resolve()
    dest_dir = Path(dest_dir_path).resolve()
    print_yml([
        {'n_cpu': n_cpu}, {'bed': bed}, {'vcfs': vcfs},
        {'dest_dir': dest_dir}
    ])
    df_bed = read_bed(
        path=str(bed), columns=['chrom', 'chromStart', 'chromEnd'],
        bgzip=bgzip, n_cpu=n_cpu
    ).assign(
        chrom=lambda d: _normalize_chrom_name(series=d['chrom'])
    )
    logger.debug(f'df_bed:{os.linesep}{df_bed}')
    bed_stem = Path(bed_path).stem
    output_count_tsvs = [
        dest_dir.joinpath(
            re.sub(r'(\.vcf|)$', f'.vcf.{bed_stem}.count.tsv', Path(p).stem)
        ) for p in vcf_paths
    ]
    for vcf, tsv in zip(vcfs, output_count_tsvs):
        df_v = _extract_alteration(
            vcf_path=str(vcf), df_bed=df_bed, bgzip=bgzip,
            include_filtered=include_filtered, n_cpu=n_cpu
        )
        logger.debug(f'df_v:{os.linesep}{df_v}')
        df_vc = df_v.group(['ref', 'alt']).size().to_frame(
            name='observed_count'
        )
        logger.debug(f'df_vc:{os.linesep}{df_vc}')
        df_vc.to_csv(str(tsv), sep='\t')


def _extract_alteration(vcf_path, df_bed, **kwargs):
    return read_vcf(
        path=vcf_path, columns=['CHROM', 'POS', 'REF', 'ALT', 'INFO'], **kwargs
    ).assign(
        chrom=lambda d: _normalize_chrom_name(series=d['CHROM']),
        pos_start=lambda d: d['POS'],
        pos_end=lambda d: np.where(
            d['ALT'].str.startswith('<'),
            d['INFO'].str.replace(
                r'^(.+;|)END=([0-9]+);.*$', r'\2'
            ).pipe(
                lambda i: i.where(i.str.isdigit(), np.nan)
            ).astype(float),
            (d['POS'] + d['REF'].str.len() - 1).astype(float)
        ).astype(int)
    ).rename(
        columns={'REF': 'ref', 'ALT': 'alt'}
    )[['chrom', 'pos_start', 'pos_end', 'ref', 'alt']].pipe(
        lambda d: d[
            d[['chrom', 'pos_start', 'pos_end']].apply(
                lambda r: (
                    (df_bed['chrom'] == r[0])
                    & (df_bed['chromStart'] < r[1])
                    & (df_bed['chromEnd'] >= r[2])
                ).any(),
                axis=1
            )
        ]
    )


def _normalize_chrom_name(series):
    return series.str.replace(r'^(chr|)', 'chr', flags=re.IGNORECASE)
