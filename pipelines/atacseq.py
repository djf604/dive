__author__ = 'Dominic Fitzgerald'
import os
import subprocess
import re
from dive.components import Software, Parameter, Redirect, Pipe

READ1 = 0
FIRST_CHAR = 0


def count_gzipped_lines(filepath):
    zcat = subprocess.Popen(['zcat', filepath], stdout=subprocess.PIPE)
    num_lines = subprocess.check_output(['wc', '-l'], stdin=zcat.stdout)
    return num_lines.strip()


def run_pipeline(read_pairs, options):
    # Instantiate options
    output_dir = options['output_dir']
    logs_dir = options['logs_dir']
    lib_prefix = options['lib_prefix']
    step = options['step']
    config = options['config']
    run_is_paired_end = options['run_is_paired_end']

    try:
        forward_adapter = options['extra_info']['forward_adapter']
        reverse_adapter = options['extra_info']['reverse_adapter']
    except KeyError, e:
        # TODO Make better exception message
        raise KeyError('Some needed extra_info not given.')

    # Keep list of items to delete
    staging_delete= ['tmp']
    bwa_bam_outs = []
    qc_data = {
        'total_raw_reads_counts': [],
        'trimmed_reads_counts': [],
        # TODO Find a better way to store FastQC results
        'num_reads_mapped': [],
        'percent_duplicate_reads': '0',
        'num_unique_reads_mapped': [],  # TODO This isn't implemented
        'num_mtDNA_reads_mapped': [],  # TODO This isn't implemented
        'num_reads_mapped_after_filtering': '0',  # TODO This isn't implemented
        'num_peaks_called': '-1',
        # TODO Get number of peaks in annotation sites
    }

    # Instantiate software instances
    cutadapt = Software('cutadapt', config['cutadapt']['path'])
    fastqc = Software('FastQC', config['fastqc']['path'])
    bwa_aln = Software('BWA aln', config['bwa']['path'] + ' aln')
    bwa_sampe = Software('BWA sampe', config['bwa']['path'] + ' sampe')
    samtools_view = Software('samtools view',
                             config['samtools']['path'] + ' view')
    samtools_flagstat = Software('samtools flagstat',
                                 config['samtools']['path'] + ' flagstat')
    novosort = Software('novosort', config['novosort']['path'])
    picard_mark_dup = Software('Picard MarkDuplicates',
                               config['picard']['path'] + ' MarkDuplicates')
    picard_insert_metrics = Software('Picard CollectInsertSizeMetrics',
                                     config['picard']['path'] + ' CollectInsertSizeMetrics')
    bedtools_bamtobed = Software('bedtools bamtobed',
                        config['bedtools']['path'] + ' bamtobed')
    sortbed = Software('sortBed', config['sortBed']['path'])
    bedtools_intersect = Software('bedtools intersect',
                                  config['bedtools']['path'] + ' intersect')
    homer_maketagdir = Software('HOMER makeTagDirectory',
                                config['makeTagDirectory']['path'])
    homer_findpeaks = Software('HOMER findPeaks', config['findPeaks']['path'])
    homer_pos2bed = Software('HOMER pos2bed', config['pos2bed']['path'])

    if step <= 1:
        for i, read_pair in enumerate(read_pairs):
            read1, read2 = read_pair.split(':')

            # Get raw fastq read counts
            qc_data['total_raw_reads_counts'].append([
                str(int(count_gzipped_lines(read1))/4),
                str(int(count_gzipped_lines(read2))/4)
            ])

            trimmed_read1_filename = os.path.join(output_dir,
                                                  lib_prefix + '_{}_read1.trimmed.fastq.gz'.format(i))
            trimmed_read2_filename = os.path.join(output_dir,
                                                  lib_prefix + '_{}_read2.trimmed.fastq.gz'.format(i))

            cutadapt.run(
                Parameter('--quality-base=33'),
                Parameter('--minimum-length=5'),
                Parameter('-q', '30'),  # Minimum quality score
                Parameter('--output={}'.format(trimmed_read1_filename)),
                Parameter('--paired-output={}'.format(trimmed_read2_filename)),
                Parameter('-a', forward_adapter if forward_adapter else 'ZZZ'),
                Parameter('-A', reverse_adapter if reverse_adapter else 'ZZZ'),
                Parameter(read1),
                Parameter(read2),
                Redirect(type='1>', dest=os.path.join(logs_dir, 'cutadapt.summary.log'))
            )

            # Get trimmed fastq read counts
            qc_data['trimmed_reads_counts'].append([
                str(int(count_gzipped_lines(trimmed_read1_filename))/4),
                str(int(count_gzipped_lines(trimmed_read2_filename))/4)
            ])

            staging_delete.extend([trimmed_read1_filename, trimmed_read2_filename])
            read_pairs[i] = ':'.join([trimmed_read1_filename, trimmed_read2_filename])

    if step <= 2:
        # Make FastQC directory
        fastqc_output_dir = os.path.join(output_dir, 'fastqc')
        subprocess.call(['mkdir', '-p', fastqc_output_dir])
        for i, read_pair in enumerate(read_pairs):
            for read in read_pair.split(':'):
                fastqc.run(
                    Parameter('--outdir={}'.format(fastqc_output_dir)),
                    Parameter(read)
                )

                bwa_aln.run(
                    Parameter('-t', config['bwa']['threads']),
                    Parameter(config['bwa']['index-dir']),
                    Parameter(read),
                    Redirect(type='>', dest='{}.sai'.format(read))
                )

                staging_delete.append('{}.sai'.format(read))

    if step <= 3:
        for i, read_pair in enumerate(read_pairs):
            read1, read2 = read_pair.split(':')

            bwa_sampe.run(
                Parameter('-a', '2000'),  # Maximum insert size
                Parameter('-n', '1'),
                Parameter(config['bwa']['index-dir']),
                Parameter('{}.sai'.format(read1)),
                Parameter('{}.sai'.format(read2)),
                Parameter(read1),
                Parameter(read2),
                Redirect(type='2>', dest=os.path.join(logs_dir, 'bwa_sampe.log')),
                Pipe(
                    samtools_view.cmd(
                        Parameter('-hSb'),
                        Parameter('-o', '{}.{}.bam'.format(lib_prefix, i)),
                        Parameter('-')  # Get input from stdin
                    )
                )
            )

            bwa_bam_outs.append('{}.{}.bam'.format(lib_prefix, i))
            # staging_delete.append('{}.{}.bam'.format(lib_prefix, i))

    if step <= 4:
        for i, bwa_bam in enumerate(bwa_bam_outs):
            samtools_flagstat.run(
                Parameter(bwa_bam),
                Redirect(type='>', dest=bwa_bam + '.flagstat')
            )

            # Get number of mapped reads from this BAM
            with open(bwa_bam + '.flagstat') as flagstats:
                flagstats_contents = flagstats.read()
                target_line = re.search(r'(\d+) \+ \d+ mapped', flagstats_contents)
                if target_line is not None:
                    qc_data['num_reads_mapped'].append(str(int(target_line.group(1))/2))
                else:
                    qc_data['num_reads_mapped'].append('0')

        novosort.run(
            Parameter('--threads', config['novosort']['threads']),
            Parameter('--tmpcompression', '6'),
            Parameter('--tmpdir', 'tmp'),
            Parameter(' '.join([bam for bam in bwa_bam_outs])),
            Redirect(type='>', dest='{}.sortmerged.bam'.format(lib_prefix)),
            Redirect(type='2>', dest=os.path.join(logs_dir, 'novosort.log'))
        )

        # TODO Do we still want to filter to only uniquely mapped reads?

        # This would be the combined uniquely mapped reads and unmapped reads
        # samtools_view.run(
        #     Parameter('-b'),
        #     Parameter('-F', '268'),
        #     Parameter('-q', '10'),
        #     Parameter('-o', '{}.unique.unmappedrm.bam'.format(lib_prefix)),
        #     Parameter('{}.sortmerged.bam'.format(lib_prefix))
        # )

        # Mark and remove duplicates
        markduplicates_metrics_filepath = os.path.join(logs_dir,
                                                       'mark_dup.metrics')
        picard_mark_dup.run(
            Parameter('INPUT={}.sortmerged.bam'.format(lib_prefix)),
            Parameter('OUTPUT={}.duprm.bam'.format(lib_prefix)),
            Parameter('TMP_DIR=tmp'),
            Parameter('METRICS_FILE={}'.format(markduplicates_metrics_filepath)),
            Parameter('REMOVE_DUPLICATES=true'),
            Redirect(type='&>', dest=os.path.join(logs_dir, 'mark_dup.log'))
        )

        # Get percent duplicates
        with open(markduplicates_metrics_filepath) as markdup_metrics:
            for line in markdup_metrics:
                if line[FIRST_CHAR] == '#':
                    continue
                record = line.strip().split('\t')
                if len(record) == 9:
                    if re.match(r'\d+', record[7]) is not None:
                        qc_data['percent_duplicate_reads'] = record[7]

        # Filter down to uniquely mapped reads
        samtools_view.run(
            Parameter('-b'),
            Parameter('-F', '256'),
            Parameter('-q', '10'),
            Parameter('-o', '{}.unique.bam'.format(lib_prefix)),
            Parameter('{}.duprm.bam'.format(lib_prefix))
        )

        # Remove unmapped reads
        samtools_view.run(
            Parameter('-b'),
            Parameter('-F', '12'),
            Parameter('-o', '{}.unmappedrm.bam'.format(lib_prefix)),
            Parameter('{}.unique.bam'.format(lib_prefix))
        )

        # Stage delete for temporary files
        staging_delete.extend([
            '{}.sortmerged.bam'.format(lib_prefix),
            '{}.unique.bam'.format(lib_prefix),
            '{}.duprm.bam'.format(lib_prefix),
            '{}.unmappedrm.bam'.format(lib_prefix)
        ])

    if step <= 5:
        bedtools_intersect.run(
            Parameter('-v'),
            Parameter('-abam', '{}.unmappedrm.bam'.format(lib_prefix)),
            Parameter('-b', config['bedtools']['blacklist-bed']),
            Parameter('-f', '0.5'),
            Redirect(type='>', dest='{}.processed.bam'.format(lib_prefix))
        )

        picard_insert_metrics.run(
            Parameter('INPUT={}.processed.bam'.format(lib_prefix)),
            Parameter('OUTPUT={}'.format(os.path.join(logs_dir, lib_prefix + '.insert_size.metrics'))),
            Parameter('HISTOGRAM_FILE={}'.format(os.path.join(logs_dir, lib_prefix + '.insert_size.pdf')))
        )

        bedtools_bamtobed.run(
            Parameter('-i', '{}.processed.bam'.format(lib_prefix)),
            Redirect(type='>', dest='{}.processed.bed'.format(lib_prefix))
        )

    if step <= 6:
        homer_tagdir = '{}_tagdir'.format(lib_prefix)

        # Populate HOMER tag directory
        homer_maketagdir.run(
            Parameter(homer_tagdir),
            Parameter('-format', 'bed'),
            Parameter('{}.processed.bed'.format(lib_prefix)),
            Redirect(type='&>', dest=os.path.join(logs_dir, 'maketagdir.log'))
        )

        # Run HOMER peak calling program
        homer_findpeaks.run(
            Parameter(homer_tagdir),
            Parameter('-fragLength', '0'),
            Parameter('-fdr', '0.01'),
            Parameter('-localSize', '50000'),
            Parameter('-o', 'auto'),
            Parameter('-style', 'dnase'),
            Parameter('-size', '150'),
            Parameter('-minDist', '50'),
            Redirect(type='&>', dest=os.path.join(logs_dir, 'findpeaks.log'))
        )

        # Convert HOMER peaks file to bed format
        homer_pos2bed.run(
            Parameter(os.path.join(homer_tagdir, 'peaks.txt')),
            Redirect(type='>', dest='{}.unsorted.peaks.bed'.format(lib_prefix)),
            Redirect(type='2>', dest=os.path.join(logs_dir, 'pos2bed.log'))
        )

        # Sort called peaks bed file
        sortbed.run(
            Parameter('-i', '{}.unsorted.peaks.bed'.format(lib_prefix)),
            Redirect(type='>', dest='{}.peaks.bed'.format(lib_prefix))
        )

        staging_delete.append('{}.unsorted.peaks.bed'.format(lib_prefix))

        with open(os.path.join(logs_dir, 'qc_metrics.txt'), 'w') as qc_data_file:
            qc_data_file.write(str(qc_data) + '\n')

        # Get number of called peaks
        # meta_data['num_peaks'] = subprocess.check_output(['wc', '-l',
        #                                                   '{}.peaks.bed'.format(lib_prefix)])

    for delete_file in staging_delete:
        subprocess.call(['rm', '-rf', delete_file])
        # Commit




