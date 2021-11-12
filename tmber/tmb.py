#!/usr/bin/env python

import logging
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

from .util import print_log, read_bed, read_vcf


def calculate_tmb(vcf_paths, bed_path, dest_dir_path='.', bedtools='bedtools',
                  bgzip='bgzip', include_filtered=False, n_cpu=1):
    logger = logging.getLogger(__name__)
    vcfs = [Path(p).resolve() for p in vcf_paths]
    bed = Path(bed_path).resolve()
    dest_dir = Path(dest_dir_path).resolve()
    df_bed = read_bed(
        path=str(bed), columns=['chrom', 'chromStart', 'chromEnd'],
        bedtools=bedtools
    ).assign(
        chrom=lambda d: _normalize_chrom_name(series=d['chrom'])
    )
    logger.debug(f'df_bed:{os.linesep}{df_bed}')
    bed_stem = Path(bed_path).stem
    df_vc = pd.concat(
        [
            _extract_alteration(
                vcf_path=str(v), df_bed=df_bed, bgzip=bgzip,
                include_filtered=include_filtered, n_cpu=n_cpu
            ).groupby(['ref', 'alt']).size().to_frame(
                name='observed_count'
            ).reset_index().assign(
                vcf=v.name,
                variant_type=lambda d: d[['ref', 'alt']].apply(
                    lambda r: _determine_sequence_ontology(ref=r[0], alt=r[1]),
                    axis=1
                )
            ) for v in vcfs
        ],
        ignore_index=True, sort=False
    ).set_index(['vcf', 'variant_type', 'ref', 'alt']).sort_index()
    logger.debug(f'df_vc:{os.linesep}{df_vc}')
    output_tsv = dest_dir.joinpath(f'{bed_stem}.tmber.tsv')
    print_log(f'Write a TSV file:\t{output_tsv}')
    df_vc.to_csv(output_tsv, sep='\t')


def _determine_sequence_ontology(ref, alt):
    if {'[', ']'} & set(alt):
        return 'structural_variant'
    elif alt.startswith(('<DEL>', '<DEL:')):
        return 'deletion'
    elif alt.startswith(('<INS>', '<INS:')):
        return 'insertion'
    elif alt.startswith(('<DUP>', '<DUP:')):
        return 'duplication'
    elif alt.startswith(('<INV>', '<INV:')):
        return 'inversion'
    elif alt.startswith(('<CNV>', '<CNV:')):
        return 'copy_number_variation'
    elif alt == '.':
        return 'no_sequence_alteration'
    elif alt == '*':
        return 'deletion'
    else:
        alt0 = alt.split(',', maxsplit=1)[0]
        if len(ref) == 1 and len(alt0) == 1:
            return 'SNV'
        elif ref[0] == alt0[0] and len(ref) > 1 and len(alt0) == 1:
            return 'deletion'
        elif ref[0] == alt0[0] and len(ref) == 1 and len(alt0) > 1:
            return 'insertion'
        elif len(ref) > 1 and len(alt0) > 1:
            return 'delins'
        else:
            return ''


def _extract_alteration(vcf_path, df_bed, **kwargs):
    return read_vcf(
        path=vcf_path, columns=['CHROM', 'POS', 'REF', 'ALT', 'INFO'], **kwargs
    ).assign(
        chrom=lambda d: _normalize_chrom_name(series=d['CHROM']),
        pos_start=lambda d: d['POS'],
        pos_end=lambda d: np.where(
            d['ALT'].str.startswith('<'),
            d['INFO'].str.replace(
                r'^(.+;|)END=([0-9]+);.*$', r'\2', regex=True
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
    return series.str.replace(
        r'^(chr|)', 'chr', regex=True, flags=re.IGNORECASE
    )
