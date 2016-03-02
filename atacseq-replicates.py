import os
import itertools
import subprocess
import sys
import uuid
import multiprocessing as mp
from dive.components import Software, Parameter, Redirect


def run_pipeline(bam1_path, bam2_path, options):
    # Instantiate options
    minimum_filter = options['minimum_filter']
    output_dir = options['output_dir']
    bam1_bid = options['bam1_bid']
    bam2_bid = options['bam2_bid']
    output_stub = options['output_stub']
    count_path = options['count_window_py_path']
    plot_path = options['plot_replicates_R_path']
    bam1_path = options['bam1_path']
    bam2_path = options['bam2_path']
    bam1_mapped_reads = options['bam1_mapped_reads']
    bam2_mapped_reads = options['bam2_mapped_reads']

    # Make output dir
    subprocess.call(['mkdir', '-p', output_dir])

    # Instantiate software
    count_window_py = Software('count_window.py', count_path)
    plot_replicates_R = Software('plot_replicates.R', plot_path)

    # Create downsampled BAM
    downsampled_bam_filename = str(uuid.uuid4()) + '.bam'
    if bam1_mapped_reads < bam2_mapped_reads:
        # bam2 is bigger, needs to be downsampled
        downsample_ratio = float(bam1_mapped_reads)/float(bam2_mapped_reads)
        print 'Downsample ratio: {}'.format(str(downsample_ratio))
        print 'Running samtools:'
        print 'samtools view -bs {} {} > {}'.format(str(downsample_ratio), bam2_path, os.path.join(os.path.dirname(bam2_path), downsampled_bam_filename))
        subprocess.call('samtools view -bs {} {} > {}'.format(str(downsample_ratio), bam2_path, os.path.join(os.path.dirname(bam2_path), downsampled_bam_filename)), shell=True)
        bam2_path = os.path.join(os.path.dirname(bam2_path), downsampled_bam_filename)
        print 'Creating index for bam2'
        subprocess.call('samtools index {}'.format(bam2_path), shell=True)
        print 'Reassigned bam2 to {}'.format(bam2_path)
    else:
        # bid1 is bigger, needs to be downsampled
        downsample_ratio = float(bam2_mapped_reads)/float(bam1_mapped_reads)
        print 'Downsample ratio: {}'.format(str(downsample_ratio))
        print 'Running samtools:'
        print 'samtools view -bs {} {} > {}'.format(str(downsample_ratio), bam1_path, os.path.join(os.path.dirname(bam1_path), downsampled_bam_filename))
        subprocess.call('samtools view -bs {} {} > {}'.format(str(downsample_ratio), bam1_path, os.path.join(os.path.dirname(bam1_path), downsampled_bam_filename)), shell=True)
        bam1_path = os.path.join(os.path.dirname(bam1_path), downsampled_bam_filename)
        print 'Creating index for bam1'
        subprocess.call('samtools index {}'.format(bam1_path), shell=True)
        print 'Reassigned bam1 to {}'.format(bam1_path)

    count_window_py.run(
        Parameter(bam1_path),
        Parameter(bam2_path),
        Parameter(minimum_filter),
        Redirect(type='>', dest=os.path.join(output_dir, 'counts.tsv'))
    )

    main_lab = '\'Alignment Counts Genome-wide {}-{}\''
    x_lab = '\'{} reads, log2\''
    y_lab = '\'{} reads, log2\''
    plot_replicates_R.run(
        Parameter(os.path.join(output_dir, 'counts.tsv')),
        Parameter(main_lab.format(bam1_bid, bam2_bid)),
        Parameter(x_lab.format(bam1_bid)),
        Parameter(y_lab.format(bam2_bid)),
        Parameter(os.path.join(output_dir, output_stub))
    )

    subprocess.call('rm {}'.format(os.path.join(os.path.dirname(bam1_path), downsampled_bam_filename + '*')), shell=True)

replicates = {
    '1709': {
        'bam_path': '/mnt/cinder/PEC/techreps/pairwise/bams/2015-1709_150814_SN1070_0405_BHVKJTADXX_1.bl.unmappedrm.duprm.sorted.unique.bam',
        'mapped_reads': 51482023
    },
    '2662': {
        'bam_path': '/mnt/cinder/PEC/techreps/pairwise/bams/2015-2662_151103_SN1070_0432.bl.duprm.unmappedrm.unique.sorted.bam',
        'mapped_reads': 70491523
    },
    '2047': {
        'bam_path': '/mnt/cinder/PEC/techreps/pairwise/bams/2015-2047_150814_SN1070_0405_BHVKJTADXX_2.bl.unmappedrm.duprm.sorted.unique.bam',
        'mapped_reads': 51238843
    },
    '2850': {
        'bam_path': '/mnt/cinder/PEC/techreps/pairwise/bams/2015-2850_151023_SN1070_0425.bl.duprm.unmappedrm.unique.sorted.bam',
        'mapped_reads': 105923319
    },
    '2657': {
        'bam_path': '/mnt/cinder/PEC/techreps/pairwise/bams/2015-2657_151023_SN1070_0425.bl.duprm.unmappedrm.unique.sorted.bam',
        'mapped_reads': 100827412
    },
    '2658': {
        'bam_path': '/mnt/cinder/PEC/techreps/pairwise/bams/2015-2658_151103_SN1070_0432.bl.duprm.unmappedrm.unique.sorted.bam',
        'mapped_reads': 54688203
    },
    '2848': {
        'bam_path': '/mnt/cinder/PEC/techreps/pairwise/bams/2015-2848_151023_SN1070_0425.bl.duprm.unmappedrm.unique.sorted.bam',
        'mapped_reads': 65027346
    },
    '2849': {
        'bam_path': '/mnt/cinder/PEC/techreps/pairwise/bams/2015-2849_151103_SN1070_0432.bl.duprm.unmappedrm.unique.sorted.bam',
        'mapped_reads': 27604805
    }
}

replicates_bams = {
    '1709':'/mnt/cinder/dfitzgeraldSCRATCH/analysis/PsychENCODE/RAW/2015-1709/2015-1709_150814_SN1070_0405_BHVKJTADXX_1.bl.unmappedrm.duprm.sorted.unique.bam',
    '2662':'/mnt/cinder/PEC/PsychENCODE/RAW/2015-2662/2015-2662_151103_SN1070_0432.bl.duprm.unmappedrm.unique.sorted.bam',
    '2047':'/mnt/cinder/dfitzgeraldSCRATCH/analysis/PsychENCODE/RAW/2015-2047/2015-2047_150814_SN1070_0405_BHVKJTADXX_2.bl.unmappedrm.duprm.sorted.unique.bam',
    '2850':'/mnt/cinder/PEC/PsychENCODE/RAW/2015-2850_save/2015-2850_151023_SN1070_0425.bl.duprm.unmappedrm.unique.sorted.bam',
    '2657':'/mnt/cinder/PEC/PsychENCODE/RAW/2015-2657/2015-2657_151023_SN1070_0425.bl.duprm.unmappedrm.unique.sorted.bam',
    '2658':'/mnt/cinder/PEC/PsychENCODE/RAW/2015-2658/2015-2658_151103_SN1070_0432.bl.duprm.unmappedrm.unique.sorted.bam',
    '2848':'/mnt/cinder/PEC/PsychENCODE/RAW/2015-2848/2015-2848_151023_SN1070_0425.bl.duprm.unmappedrm.unique.sorted.bam',
    '2849':'/mnt/cinder/PEC/PsychENCODE/RAW/2015-2849/2015-2849_151103_SN1070_0432.bl.duprm.unmappedrm.unique.sorted.bam'
}

replicate_bids = [
    '1709',
    '2662',
    '2047',
    '2850',
    '2657',
    '2658',
    '2848',
    '2849'
]
# Determine all possible permutations of replicates, without replacement
replicate_pairs = list(itertools.combinations(replicate_bids, 2))

output_dir = '/mnt/cinder/PEC/techreps/pairwise/without_filter'
count_path = 'python /home/ubuntu/scripts/count_window.py'
plot_path = 'Rscript /home/ubuntu/scripts/plot_replicates.R'

# Sort replicate pairs, lower BID always comes first
pool = mp.Pool()
for pair in replicate_pairs:
    bid1 = pair[0]
    bid2 = pair[1]
    if bid1 > bid2:
        bid1, bid2 = bid2, bid1

    bam1_path = replicates[bid1]['bam_path']
    bam2_path = replicates[bid2]['bam_path']
    bam1_mapped_reads = int(replicates[bid1]['mapped_reads'])
    bam2_mapped_reads = int(replicates[bid2]['mapped_reads'])

    print 'BAM1: {}'.format(bam1_path)
    print 'BAM1 Mapped Reads: {}'.format(str(bam1_mapped_reads))
    print 'BAM2: {}'.format(bam2_path)
    print 'BAM2 Mapped Reads: {}'.format(str(bam2_mapped_reads))

    # if bam1_mapped_reads < bam2_mapped_reads:
    #     # bam2 is bigger, needs to be downsampled
    #     downsample_ratio = float(bam1_mapped_reads)/float(bam2_mapped_reads)
    #     print 'Downsample ratio: {}'.format(str(downsample_ratio))
    #     print 'Running samtools:'
    #     print 'samtools view -bs {} {} > {}'.format(str(downsample_ratio), bam2_path, os.path.join(os.path.dirname(bam2_path), 'downsampled.bam'))
    #     subprocess.call('samtools view -bs {} {} > {}'.format(str(downsample_ratio), bam2_path, os.path.join(os.path.dirname(bam2_path), 'downsampled.bam')), shell=True)
    #     bam2_path = os.path.join(os.path.dirname(bam2_path), 'downsampled.bam')
    #     print 'Creating index for bam2'
    #     subprocess.call('samtools index {}'.format(bam2_path), shell=True)
    #     print 'Reassigned bam2 to {}'.format(bam2_path)
    # else:
    #     # bid1 is bigger, needs to be downsampled
    #     downsample_ratio = float(bam2_mapped_reads)/float(bam1_mapped_reads)
    #     print 'Downsample ratio: {}'.format(str(downsample_ratio))
    #     print 'Running samtools:'
    #     print 'samtools view -bs {} {} > {}'.format(str(downsample_ratio), bam1_path, os.path.join(os.path.dirname(bam1_path), 'downsampled.bam'))
    #     subprocess.call('samtools view -bs {} {} > {}'.format(str(downsample_ratio), bam1_path, os.path.join(os.path.dirname(bam1_path), 'downsampled.bam')), shell=True)
    #     bam1_path = os.path.join(os.path.dirname(bam1_path), 'downsampled.bam')
    #     print 'Creating index for bam1'
    #     subprocess.call('samtools index {}'.format(bam1_path), shell=True)
    #     print 'Reassigned bam1 to {}'.format(bam1_path)

    pool.apply_async(run_pipeline, args=(bam1_path, bam2_path, {
            'minimum_filter': sys.argv[1],
            'output_dir': os.path.join(output_dir, '{}-{}'.format(bid1, bid2)),
            'bam1_bid': bid1,
            'bam2_bid': bid2,
            'output_stub': '{}-{}.replicates'.format(bid1, bid2),
            'count_window_py_path': count_path,
            'plot_replicates_R_path': plot_path,
            'bam1_path': bam1_path,
            'bam2_path': bam2_path,
            'bam1_mapped_reads': bam1_mapped_reads,
            'bam2_mapped_reads': bam2_mapped_reads
    }))


    # run_pipeline(bam1_path,
    #              bam2_path,
    #              {
    #                  'minimum_filter': sys.argv[1],
    #                  'output_dir': os.path.join(output_dir, '{}-{}'.format(bid1, bid2)),
    #                  'bam1_bid': bid1,
    #                  'bam2_bid': bid2,
    #                  'output_stub': '{}-{}.replicates'.format(bid1, bid2),
    #                  'count_window_py_path': count_path,
    #                  'plot_replicates_R_path': plot_path
    #              }
    # )
pool.close()
pool.join()
