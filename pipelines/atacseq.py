__author__ = 'Dominic Fitzgerald'
import os
import subprocess
from dive.components import Software, Parameter, Redirect


def run_pipeline(reads, options):
    # Instantiate options
    output_dir = options['output_dir']
    logs_dir = options['logs_dir']
    lib_prefix = options['lib_prefix']
    step = options['step']
    config = options['config']
    run_is_paired_end = options['run_is_paired_end']

    cutadapt = Software('cutadapt', config['cutadapt']['path'])
    fastqc = Software('FastQC', config['fastqc']['path'])
    bwa_aln = Software('BWA aln', config['bwa']['path'] + ' aln')
    bwa_sampe = Software('BWA sampe', config['bwa']['path'] + ' sampe')
