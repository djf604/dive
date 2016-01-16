__author__ = 'Dominic Fitzgerald'
import os
import subprocess
import re
import sys

try:
    start = int(sys.argv[1])
except:
    start = 1
    stop = 9999
try:
    stop = int(sys.argv[2])
except:
    stop = 9999

ucla_dir = '/lustre/beagle2/djf604/synapse/UCLA_R01MH094714/RAW'
ucla_out = '/lustre/beagle2/djf604/workspace/analysis/UCLA_R01MH094714'
ucla_pbs = '/lustre/beagle2/djf604/software/PEC/dive/pbs/trendy-24.pbs'

ucla_files = os.listdir(ucla_dir)

for i, ucla_pair in enumerate(ucla_files):
    if (i + 1) < start:
        continue
    if (i + 1) > stop:
        break
    os.chdir(os.path.join(ucla_dir, ucla_pair))
    read1, read2 = os.listdir('.')[:2]
    if re.search(r'\.R1\.', read2) is not None:
        read1, read2 = read2, read1

    ucla_pair_out = os.path.join(ucla_out, ucla_pair)
    # subprocess.call(['mkdir', '-p', ucla_pair_out])

    ucla_args = []
    ucla_args.append('READS="{}"'.format(':'.join([
        os.path.join(ucla_dir, ucla_pair, read1),
        os.path.join(ucla_dir, ucla_pair, read2)]))
    )
    ucla_args.append('OUTDIR="{}"'.format(ucla_pair_out))
    ucla_args.append('LIBNAME="{}"'.format(ucla_pair))
    ucla_args.append('FORWARD_ADAPTER=""')
    ucla_args.append('REVERSE_ADAPTER=""')
    ucla_args.append('SAILFISH_LIBTYPE="{}"'.format('ISF'))

    # subprocess.call(['qsub', '-v', ','.join(ucla_args), ucla_pbs])
    # print ' '.join(['qsub', '-v', ','.join(ucla_args), ucla_pbs])
