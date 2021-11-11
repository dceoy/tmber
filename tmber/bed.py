#!/usr/bin/env python

import logging
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from .util import print_log, print_yml, read_fasta


def create_bed_from_fa(fa_path, bed_path, bgzip='bgzip', uppercase=False,
                       human_autosome=False, n_cpu=1):
    logger = logging.getLogger(__name__)
    fa = Path(fa_path).resolve()
    bed = Path(bed_path).resolve()
    autosomes = {f'chr{i}' for i in range(1, 23)}
    print_yml([{'n_cpu': n_cpu}, {'fa': fa}, {'bed': bed}])
    with ProcessPoolExecutor(max_workers=n_cpu) as x:
        fs = [
            x.submit(
                _identify_target_region,
                k, v.seq, ('AGCT' if uppercase else 'AGCTagct')
            ) for k, v in read_fasta(path=fa_path).items()
            if human_autosome and k in autosomes
        ]
        df_bed = pd.concat(
            [f.result() for f in as_completed(fs)], ignore_index=True,
            sort=False
        ).sort_values(['chrom', 'chromStart', 'chromEnd'])
    logger.debug(f'df_bed:{os.linesep}{df_bed}')
    print_log(f'Write a BED file: {bed}')
    df_bed.to_csv(bed, sep='\t', header=False, index=False)


def _identify_target_region(chrom, sequence, target_cases='AGCT'):
    logger = logging.getLogger(__name__)
    bseq = pd.Series(list(sequence)).isin(set(target_cases)).astype(int)
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
                    *(s[s.diff() == 1].index + 1),
                    *([len(s)] if s.iloc[-1] == 1 else list())
                ]
            })
        )
    else:
        logger.info(f'No region to extract: {chrom}')
