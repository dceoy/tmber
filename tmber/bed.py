#!/usr/bin/env python

import logging
import os
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from .util import print_log, read_fasta_and_generate_seq


def create_bed_from_fa(fa_path, dest_dir_path, bgzip='bgzip',
                       human_autosome=False, target_letters='ACGT', n_cpu=1):
    logger = logging.getLogger(__name__)
    target_letter_set = set(target_letters)
    print_log('Set target letters:\t{}'.format(target_letter_set))
    fa = Path(fa_path).resolve()
    bed = Path(dest_dir_path).resolve().joinpath(
        re.sub(r'\.(gz|bz2|bgz)', '', Path(fa_path).name)
        + '.' + ''.join(target_letter_set) + '.bed'
    )
    autosomes = {f'chr{i}' for i in range(1, 23)}
    fs = list()
    with ProcessPoolExecutor(max_workers=n_cpu) as x:
        for chrom, seq in read_fasta_and_generate_seq(path=str(fa)):
            seq_len = len(seq)
            if human_autosome and chrom in autosomes:
                print_log(
                    f'Detect the target letters:\t{chrom}\t({seq_len} bp)'
                )
                fs.append(
                    x.submit(
                        _identify_target_region, chrom, seq, target_letter_set
                    )
                )
            else:
                logger.info(f'Skip detection: {chrom} ({seq_len} bp)')
        f_results = [f.result() for f in as_completed(fs)]
    df_bed = pd.concat(
        f_results, ignore_index=True, sort=False
    ).sort_values(['chrom', 'chromStart', 'chromEnd'])
    logger.debug(f'df_bed:{os.linesep}{df_bed}')
    print_log(f'Write a BED file:\t{bed}')
    df_bed.to_csv(bed, sep='\t', header=False, index=False)


def _identify_target_region(chrom, sequence, target_letter_set):
    bseq = pd.Series(list(sequence)).isin(target_letter_set).astype(int)
    if bseq.sum() > 0:
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
        logger = logging.getLogger(__name__)
        logger.info(f'Target letters not detected: {chrom}')
