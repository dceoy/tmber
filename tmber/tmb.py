#!/usr/bin/env python

import logging
import os
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

from .util import print_log, read_bed, read_vcf


def calculate_tmb(vcf_path, bed_paths, dest_dir_path='.', bedtools='bedtools',
                  bgzip='bgzip', include_filtered=False, sample_name=None,
                  min_af=None, max_af=None, n_cpu=1):
    logger = logging.getLogger(__name__)
    vcf = Path(vcf_path).resolve()
    beds = [Path(p).resolve() for p in bed_paths]
    output_tsv = Path(dest_dir_path).resolve().joinpath(
        '.'.join([
            re.sub(r'\.(gz|bz2|bgz)', '', Path(vcf_path).name),
            ('' if min_af is None and max_af is None else str(sample_name)),
            '_n_'.join([
                '{0}{1:0>2}'.format(k, int(v * 100)) for k, v
                in [('min_af', min_af), ('max_af', max_af)] if v is not None
            ]),
            'tmber.tsv'
        ])
    )
    output_tsv_columns = [
        'bed_name', 'bed_size', 'variant_type', 'ref', 'alt', 'observed_count'
    ]
    df_vcf = read_vcf(
        path=str(vcf), sample_name=sample_name, min_af=min_af, max_af=max_af,
        include_filtered=include_filtered, bgzip=bgzip, n_cpu=n_cpu
    )
    logger.debug(f'df_vcf:{os.linesep}{df_vcf}')
    if df_vcf.shape[0] == 0:
        df_tmb = pd.DataFrame(column=output_tsv_columns).set_index(
            output_tsv_columns[:-1]
        )
    else:
        df_alt = df_vcf.assign(
            chrom=lambda d: d['CHROM'].str.replace(
                r'^(chr|)', 'chr', regex=True, flags=re.IGNORECASE
            ),
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
            columns={'POS': 'pos_start', 'REF': 'ref', 'ALT': 'alt'}
        )[['chrom', 'pos_start', 'pos_end', 'ref', 'alt']].drop_duplicates()
        logger.debug(f'df_alt:{os.linesep}{df_alt}')
        fs = list()
        with ProcessPoolExecutor(max_workers=n_cpu) as x:
            for bed in beds:
                df_bed = read_bed(path=str(bed), merge=True, bedtools=bedtools)
                logger.debug(f'df_bed:{os.linesep}{df_bed}')
                assert df_bed.shape[0] > 0
                fs.append(x.submit(_tally_alt, df_alt, df_bed, bed.name))
            f_results = [f.result() for f in as_completed(fs)]
        df_tmb = pd.concat(
            f_results, ignore_index=True, sort=False
        ).assign(
            variant_type=lambda d: d[['ref', 'alt']].apply(
                lambda r: _determine_sequence_ontology(ref=r[0], alt=r[1]),
                axis=1
            )
        ).set_index([
            'bed_name', 'bed_size', 'variant_type', 'ref', 'alt'
        ]).sort_index()
    logger.debug(f'df_tmb:{os.linesep}{df_tmb}')
    print_log(f'Write a TSV file:\t{output_tsv}')
    df_tmb.to_csv(output_tsv, sep='\t')


def _tally_alt(df_alt, df_bed, bed_name):
    return df_alt[
        df_alt[['chrom', 'pos_start', 'pos_end']].apply(
            lambda r: (
                (df_bed['chrom'] == r[0])
                & (df_bed['chromStart'] < r[1])
                & (df_bed['chromEnd'] >= r[2])
            ).any(),
            axis=1
        )
    ].groupby(['ref', 'alt']).size().to_frame(
        name='observed_count'
    ).reset_index().assign(
        bed_name=bed_name,
        bed_size=(df_bed['chromEnd'] - df_bed['chromStart']).sum()
    )


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
