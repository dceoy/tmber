#!/usr/bin/env python

from pathlib import Path

from .util import print_yml


def calculate_tmb(vcf_paths, bed_path, dest_dir_path='.', bgzip='bgzip',
                  n_cpu=1, log_txt_path=None, quiet=False,
                  executable='/bin/bash'):
    print_yml([{'n_cpu': n_cpu}, {'bed': bed_path}, {'vcf': vcf_paths}])
    vcfs = [Path(p).resolve() for p in vcf_paths]
    bed = Path(bed_path).resolve()
    print(bed, vcfs)
