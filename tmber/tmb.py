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
                  bgzip='bgzip', include_filtered=False, n_cpu=1):
    logger = logging.getLogger(__name__)
    vcf = Path(vcf_path).resolve()
    beds = [Path(p).resolve() for p in bed_paths]
    output_tsv = Path(dest_dir_path).resolve().joinpath(
        re.sub(r'\.(gz|bz2|bgz)', '', Path(vcf_path).name) + '.tmber.tsv'
    )
    df_vcf = read_vcf(
        path=str(vcf), bgzip=bgzip, include_filtered=include_filtered,
        n_cpu=n_cpu
    )[['CHROM', 'POS', 'REF', 'ALT', 'INFO']]
    logger.debug(f'df_vcf:{os.linesep}{df_vcf}')
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
    )[['chrom', 'POS', 'pos_end', 'REF', 'ALT']].rename(
        columns={'POS': 'pos_start', 'REF': 'ref', 'ALT': 'alt'}
    ).drop_duplicates()
    logger.debug(f'df_alt:{os.linesep}{df_alt}')
    fs = list()
    with ProcessPoolExecutor(max_workers=n_cpu) as x:
        for bed in beds:
            df_bed = read_bed(path=str(bed), bedtools=bedtools)
            logger.debug(f'df_bed:{os.linesep}{df_bed}')
            fs.append(x.submit(_tally_alterations, df_alt, df_bed, bed.name))
        df_vc = pd.concat(
            [f.result() for f in as_completed(fs)], ignore_index=True,
            sort=False
        ).pipe(
            lambda d: d.merge(
                d[['ref', 'alt']].drop_duplicates().assign(
                    variant_type=lambda e: e.apply(
                        lambda r:
                        _determine_sequence_ontology(ref=r[0], alt=r[1]),
                        axis=1
                    )
                ),
                on=['ref', 'alt'], how='left'
            )
        ).set_index(['bed_name', 'variant_type', 'ref', 'alt']).sort_index()
    logger.debug(f'df_vc:{os.linesep}{df_vc}')
    print_log(f'Write a TSV file:\t{output_tsv}')
    df_vc.to_csv(output_tsv, sep='\t')


def _tally_alterations(df_alt, df_bed, bed_name):
    return df_alt.pipe(
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
    ).groupby(['ref', 'alt']).size().to_frame(
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
