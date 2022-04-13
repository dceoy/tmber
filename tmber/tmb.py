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
    assert vcf.is_file(), f'file not found: {vcf}'
    beds = [Path(p).resolve() for p in bed_paths]
    for bed in beds:
        assert bed.is_file(), f'file not found: {bed}'
    output_alt_tsv = Path(dest_dir_path).resolve().joinpath(
        '.'.join([
            re.sub(r'\.(gz|bz2|bgz)', '', Path(vcf_path).name),
            ('' if min_af is None and max_af is None else str(sample_name)),
            '_n_'.join([
                '{0}{1:0>3}'.format(k, int(v * 1000)) for k, v
                in [('min_af', min_af), ('max_af', max_af)] if v is not None
            ]),
            'alt.tsv'
        ])
    )
    output_tmb_tsv = output_alt_tsv.parent.joinpath(
        Path(output_alt_tsv.stem).stem + '.tmb.tsv'
    )
    df_vcf = read_vcf(
        path=str(vcf), sample_name=sample_name, min_af=min_af, max_af=max_af,
        include_filtered=include_filtered, bgzip=bgzip, n_cpu=n_cpu
    )
    logger.debug(f'df_vcf:{os.linesep}{df_vcf}')
    bed_dfs = dict()
    df_size = pd.DataFrame(columns=['bed_name', 'bed_size'])
    for b in beds:
        df_bed = read_bed(path=str(b), merge=True, bedtools=bedtools)
        logger.debug(f'df_bed:{os.linesep}{df_bed}')
        assert df_bed.shape[0] > 0
        bed_dfs[b.name] = df_bed
        bed_size = (df_bed['chromEnd'] - df_bed['chromStart']).sum()
        logger.debug(f'bed_size: {bed_size}')
        assert bed_size > 0
        df_size = pd.concat([
            df_size, pd.DataFrame([{'bed_name': b.name, 'bed_size': bed_size}])
        ])
    logger.debug(f'df_size:{os.linesep}{df_size}')
    if df_vcf.shape[0] == 0:
        logger.info('No variant detected.')
        df_var = pd.DataFrame(
            column=[
                'bed_name', 'bed_size', 'variant_type', 'ref', 'alt',
                'observed_alt_count'
            ]
        ).set_index(['bed_name', 'bed_size', 'variant_type', 'ref', 'alt'])
    else:
        df_var = df_vcf.assign(
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
        logger.debug(f'df_var:{os.linesep}{df_var}')
        fs = list()
        with ProcessPoolExecutor(max_workers=n_cpu) as x:
            fs = [
                x.submit(_tally_alt, df_var, v, k) for k, v in bed_dfs.items()
            ]
            f_results = [f.result() for f in as_completed(fs)]
        df_alt = pd.concat(
            f_results, ignore_index=True, sort=False
        ).assign(
            variant_type=lambda d: d[['ref', 'alt']].apply(
                lambda r: _determine_sequence_ontology(ref=r[0], alt=r[1]),
                axis=1
            )
        ).pipe(
            lambda d: d.merge(df_size, on='bed_name', how='left')
        ).set_index([
            'bed_name', 'bed_size', 'variant_type', 'ref', 'alt'
        ]).sort_index()
    logger.debug(f'df_alt:{os.linesep}{df_alt}')
    print_log(f'Write a TSV file:\t{output_alt_tsv}')
    df_alt.to_csv(output_alt_tsv, sep='\t')
    df_tmb = df_alt.reset_index().pipe(
        lambda d: pd.concat([
            d[
                d['variant_type'] != 'no_sequence_alteration'
            ][['bed_name', 'bed_size', 'variant_type', 'observed_alt_count']],
            *[
                df_size.assign(variant_type=t, observed_alt_count=0) for t in [
                    'SNV', 'deletion', 'insertion', 'delins',
                    'structural_variant', 'duplication', 'inversion',
                    'copy_number_variation'
                ]
            ]
        ])
    ).groupby([
        'bed_name', 'bed_size', 'variant_type'
    ])['observed_alt_count'].sum().to_frame().reset_index().pipe(
        lambda d: pd.concat([
            d,
            d.groupby([
                'bed_name', 'bed_size'
            ])['observed_alt_count'].sum().to_frame().assign(
                variant_type='total'
            ).reset_index()[d.columns]
        ])
    ).assign(
        mutations_per_mb=lambda d:
        (d['observed_alt_count'] / d['bed_size'] * 1000000)
    ).set_index(['bed_name', 'bed_size', 'variant_type'])
    logger.debug(f'df_tmb:{os.linesep}{df_tmb}')
    print_log(f'Write a TSV file:\t{output_tmb_tsv}')
    df_tmb.to_csv(output_tmb_tsv, sep='\t')


def _tally_alt(df_var, df_bed, bed_name):
    return df_var[
        df_var[['chrom', 'pos_start', 'pos_end']].apply(
            lambda r: (
                (df_bed['chrom'] == r[0])
                & (df_bed['chromStart'] < r[1])
                & (df_bed['chromEnd'] >= r[2])
            ).any(),
            axis=1
        )
    ].groupby(['ref', 'alt']).size().to_frame(
        name='observed_alt_count'
    ).reset_index().assign(bed_name=bed_name)


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
        elif len(ref) >= 1 and len(alt0) >= 1:
            return 'delins'
        else:
            raise ValueError(f'unsupported REF/ALT: "{ref}" / "{alt}"')
