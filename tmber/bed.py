#!/usr/bin/env python

import logging
import os
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from .util import print_log, read_fasta


def create_bed_from_fa(fa_path, dest_dir_path, bgzip='bgzip',
                       human_autosome=False, target_letters='ACGT', n_cpu=1):
    logger = logging.getLogger(__name__)
    fa = Path(fa_path).resolve()
    bed = Path(dest_dir_path).resolve().joinpath(
        re.sub(r'\.(gz|bz2|bgz)', '', Path(fa_path).name)
        + '.' + ''.join(set(target_letters)) + '.bed'
    )
    autosomes = {f'chr{i}' for i in range(1, 23)}
    chrom_seqs = {
        k: v for k, v in read_fasta(path=str(fa)).items()
        if human_autosome and k in autosomes
    }
    print_log('Identify target letters: {}'.format(set(target_letters)))
    with ProcessPoolExecutor(max_workers=n_cpu) as x:
        fs = [
            x.submit(_identify_target_region, k, v.seq, target_letters)
            for k, v in chrom_seqs.items()
        ]
        df_bed = pd.concat(
            [f.result() for f in as_completed(fs)], ignore_index=True,
            sort=False
        ).sort_values(['chrom', 'chromStart', 'chromEnd'])
    logger.debug(f'df_bed:{os.linesep}{df_bed}')
    print_log(f'Write a BED file: {bed}')
    df_bed.to_csv(bed, sep='\t', header=False, index=False)


def _identify_target_region(chrom, sequence, target_letters='ACGT'):
    logger = logging.getLogger(__name__)
    bseq = pd.Series(list(sequence)).isin(set(target_letters)).astype(int)
    if bseq.sum() > 0:
        logger.info(f'Identify regions to extract: {chrom}')
        return bseq.pipe(
            lambda s: pd.DataFrame({
                'chrom': chrom,
                'chromStart': [
                    *([0] if s.iloc[0] == 1 else list()),
                    *s[s.diff() == 1].index
                ],
                'chromEnd': [
                    *s[s.diff() == -1].index,
                    *([len(s)] if s.iloc[-1] == 1 else list())
                ]
            })
        )
    else:
        logger.info(f'No region to extract: {chrom}')
