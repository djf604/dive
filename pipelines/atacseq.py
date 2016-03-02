__author__ = 'Dominic Fitzgerald'
import os
import subprocess
import re
from dive.components import Software, Parameter, Redirect, Pipe

READ1 = 0


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
    meta_data = {
        'raw_counts': [],
        'trimmed_counts': [],
        'mapped_reads': [],
        'num_peaks': '-1'
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
            meta_data['raw_counts'].append([
                subprocess.check_output(['wc', '-l', read1]),
                subprocess.check_output(['wc', '-l', read2])
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
            meta_data['trimmed_counts'].append([
                subprocess.check_output(['wc', '-l', trimmed_read1_filename]),
                subprocess.check_output(['wc', '-l', trimmed_read2_filename])
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
                Redirect(type='>', dest=os.path.join(logs_dir, 'bwa_sampe.log')),
                Pipe(
                    samtools_view.cmd(
                        Parameter('-hSb'),
                        Parameter('-o', '{}.{}.bam'.format(lib_prefix, i)),
                        Parameter('-')  # Get input from stdin
                    )
                )
            )

            bwa_bam_outs.append('{}.{}.bam'.format(lib_prefix, i))
            staging_delete.append('{}.{}.bam'.format(lib_prefix, i))

    if step <= 4:
        for i, bwa_bam in enumerate(bwa_bam_outs):
            samtools_flagstat.run(
                Parameter(bwa_bam),
                Redirect(type='>', dest=bwa_bam + '.flagstat')
            )

            # Get number of mapped reads
            # Calculate as trimmed_reads * flagstats_mapped_percentage
            try:
                flagstats_percent = (re.search(r'\d+ \+ \d+ mapped \((\d+\.\d+)%:.*?\)',
                                               subprocess.check_output(['cat', bwa_bam + '.flagstat']))
                                     .group(1))
                meta_data['mapped_reads'].append(
                    str(int(int(meta_data['trimmed_counts'][i][READ1]) * float(flagstats_percent)))
                )
            except:
                pass

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
        samtools_view.run(
            Parameter('-b'),
            Parameter('-F', '268'),
            Parameter('-q', '10'),
            Parameter('-o', '{}.unique.unmappedrm.bam'.format(lib_prefix)),
            Parameter('{}.sortmerged.bam'.format(lib_prefix))
        )

        # Filter down to uniquely mapped reads
        samtools_view.run(
            Parameter('-b'),
            Parameter('-F', '256'),
            Parameter('-q', '10'),
            Parameter('-o', '{}.unique.bam'.format(lib_prefix)),
            Parameter('{}.sortmerged.bam'.format(lib_prefix))
        )

        picard_mark_dup.run(
            Parameter('INPUT={}.unique.bam'.format(lib_prefix)),
            Parameter('OUTPUT={}.duprm.bam'.format(lib_prefix)),
            Parameter('TMP_DIR=tmp'),
            Parameter('METRICS_FILE={}'.format(os.path.join(logs_dir,
                                                            'mark_dup.metrics'))),
            Parameter('REMOVE_DUPLICATES=true'),
            Redirect(type='&>', dest=os.path.join(logs_dir, 'mark_dup.log'))
        )

        # Remove unmapped reads
        samtools_view.run(
            Parameter('-b'),
            Parameter('-F', '12'),
            Parameter('-o', '{}.unmappedrm.bam'.format(lib_prefix)),
            Parameter('{}.duprm.bam'.format(lib_prefix))
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

        # Get number of called peaks
        meta_data['num_peaks'] = subprocess.check_output(['wc', '-l',
                                                          '{}.peaks.bed'.format(lib_prefix)])




