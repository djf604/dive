__author__ = 'Dominic Fitzgerald'
import os
import subprocess
import re
import sys

# Get user input parameters
try:
    start = int(sys.argv[1])
except:
    start = 1
    stop = 9999
try:
    stop = int(sys.argv[2])
except:
    stop = 9999

# Set up paths
ucla_dir = '/lustre/beagle2/djf604/synapse/UCLA_R01MH094714/RAW'
ucla_out = '/lustre/beagle2/djf604/PsychENCODE/Unified/UCLA_R01MH094714/chicago19'
ucla_pbs = '/lustre/beagle2/djf604/software/PEC/dive/pbs/chicago-19.pbs'

# Make output directory if it doesn't exist
subprocess.call(['mkdir', '-p', ucla_out])

# Get all input files
ucla_files = os.listdir(ucla_dir)

for i, ucla_pair in enumerate(ucla_files):
    # Start at user specified start
    if (i + 1) < start:
        continue
    # Stop at user specified stop
    if (i + 1) > stop:
        break
    # Change into sample directory, get filenames
    os.chdir(os.path.join(ucla_dir, ucla_pair))
    read1, read2 = os.listdir('.')[:2]
    if re.search(r'\.R1\.', read2) is not None:
        read1, read2 = read2, read1

    # Create output directory for sample, and logs directory
    ucla_pair_out = os.path.join(ucla_out, ucla_pair)
    subprocess.call(['mkdir', '-p', os.path.join(ucla_pair_out, 'logs')])

    # Set up qsub arguments
    ucla_args = []
    ucla_args.append('READS="{}"'.format(':'.join([
        os.path.join(ucla_dir, ucla_pair, read1),
        os.path.join(ucla_dir, ucla_pair, read2)]))
    )
    ucla_args.append('OUTDIR="{}"'.format(ucla_pair_out))
    ucla_args.append('LIBNAME="{}"'.format(ucla_pair))
    ucla_args.append('FORWARD_ADAPTER=""')
    ucla_args.append('REVERSE_ADAPTER=""')
    ucla_args.append('RUN_IS_STRANDED="true"')
    ucla_args.append('CUFFLINKS_LIB_TYPE="fr-firststrand"')
    ucla_args.append('HTSEQ_STRANDED="yes"')

    # Change to logs directory so qsub output will go here
    os.chdir(os.path.join(ucla_pair_out, 'logs'))
    subprocess.call(['echo', '"{} of {}"'.format(i + 1, len(ucla_files))],
                    stdout=open('order.log', 'w'))

    subprocess.call(['qsub', '-v', ','.join(ucla_args), ucla_pbs])
    print 'Submitted job {} of {}'.format(i + 1, len(ucla_files))
