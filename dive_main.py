#!/bin/python
import sys
import os
# import logging
import json
import importlib
import time
import subprocess
# from components.pipeline_software_base import PipelineSoftwareBase, SoftwareConfigService
# from components.pipeline_software import *
# import datetime
from dive.components import Software, Parameter, Redirect


def main():
    import argparse
    parser = argparse.ArgumentParser(description='RNAseq Unified Analysis Pipeline: Chicago Flavor')
    parser.add_argument('-r', '--read', action='append', help='File paths to raw data. Separate with colon for paired-end.')
    parser.add_argument('-o', '--output-dir', dest='output_dir', default='.', help='Destination for processed files.')
    parser.add_argument('-l', '--lib-prefix', dest='lib_prefix', help='Library prefix name.')
    parser.add_argument('-s', '--step', type=int, default=0, help='Start step in the pipeline. Defaults to 0 for beginning.')
    parser.add_argument('-c', '--config', help='Path to configuration JSON file.')
    parser.add_argument('-x', '--extra', action='append', help='Extra info designated for a specific pipeline. '
                                                               + 'Specify multiple times for multiple values. Format: \"{key}:{value}\"')
    # parser.add_argument('-h', '--help', help='Print this help message, and optionally print help for given pipeline module.')
    parser.add_argument('pipeline_module', help='Python module to pipeline logic.')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    user_args = parser.parse_args()
    reads = user_args.read
    output_dir = user_args.output_dir
    lib_prefix = user_args.lib_prefix
    step = user_args.step
    extra = user_args.extra if user_args.extra else []
    pipeline_module = user_args.pipeline_module

    # Print help if asked for
    # TODO

    # Make sure some reads were provided
    if len(reads) < 1:
        print 'Must provide at least one raw data file.'
        sys.exit(1)

    # Determine type of reads by user input
    run_is_paired_end = False
    if len(reads[0].split(':')) == 2:
        run_is_paired_end = True

    # Make sure there's a trailing slash on output dir
    # Create it if it doesn't exist
    output_dir = output_dir.rstrip('/')
    output_dir += '/'
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # Create logs directory
    logs_dir = os.path.join(output_dir, 'logs')
    if not os.path.isdir(logs_dir):
        os.makedirs(logs_dir)

    # Create tmp directory
    tmp_dir = os.path.join(output_dir, 'tmp')
    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir)

    # Read in configuration file
    config = json.loads(open(user_args.config, 'r').read())

    # Build extra_info list
    extra_info = {e_key.strip(): e_val.strip() for e_key, e_val in [e.split(':') for e in extra]}

    # Change to output directory
    os.chdir(output_dir)

    # Gather options into options dictionary
    options = {
        'output_dir': output_dir,
        'logs_dir': logs_dir,
        'lib_prefix': lib_prefix,
        'step': step,
        'config': config,
        'run_is_paired_end': run_is_paired_end,
        'extra_info': extra_info
    }

    # Import pipeline and run
    try:
        pipeline = importlib.import_module(pipeline_module)
        pipeline.run_pipeline(reads, options)
        print(' '.join(['>', time.strftime('%d %b %Y %H:%M:%S'), 'Pipeline ran successfully']))
    except Exception, e:
        print e

    # Remove temporary files and directory
    subprocess.call(['rm', '-rf', tmp_dir])

main()
