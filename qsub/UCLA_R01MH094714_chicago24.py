__author__ = 'Dominic Fitzgerald'
import os
import subprocess
import re

ucla_dir = '/lustre/beagle2/djf604/synapse/UCLA_R01MH094714/RAW'
ucla_out = '/lustre/beagle2/djf604/workspace/analysis/UCLA_R01MH094714'
ucla_pbs = '/lustre/beagle2/djf604/software/PEC/dive/pbs/chicago-24.pbs'

ucla_files = os.listdir(ucla_dir)

for ucla_pair in ucla_files:
    os.chdir(os.path.join(ucla_dir, ucla_pair))
    read1, read2 = os.listdir('.')[:2]
    if re.search(r'\.R1\.', read2) is not None:
        read1, read2 = read2, read1

    ucla_pair_out = os.path.join(ucla_out, ucla_pair)
    ucla_pbs_logs = os.path.join(ucla_pair_out, 'logs')
    # subprocess.call(['mkdir', '-p', ucla_pair_out])
    subprocess.call(['mkdir', '-p', ucla_pbs_logs])

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

    os.chdir(ucla_pbs_logs)
    # subprocess.call(['qsub', '-v', ','.join(ucla_args), ucla_pbs])
    print ' '.join(['qsub', '-v', ','.join(ucla_args), ucla_pbs])

