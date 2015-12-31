__author__ = 'Dominic Fitzgerald'
import os
from dive.components import Software, Parameter, Redirect


def run_pipeline(reads, options):
    # Instantiate options
    output_dir = options['output_dir']
    logs_dir = options['logs_dir']
    lib_prefix = options['lib_prefix']
    step = options['step']
    config = options['config']
    run_is_paired_end = options['run_is_paired_end']

    try:
        run_is_stranded = True if options['extra_info']['run_is_stranded'].lower() == 'true' else False
        cufflinks_lib_type = options['extra_info']['cufflinks_lib_type']
        forward_adapter = options['extra_info']['forward_adapter']
        reverse_adapter = options['extra_info']['reverse_adapter']
    except KeyError, e:
        # TODO Make better exception message
        raise KeyError('Some needed extra_info not given.')


    # Establish Software instances
    cutadapt = Software('Cutadapt', config['cutadapt']['path'])
    star = Software('STAR Two-Pass', config['STAR']['path'])
    novosort = Software('Novosort', config['novosort']['path'])
    samtools_flagstat = Software('Samtools Flagstat', config['samtools_flagstat']['path'])
    cufflinks = Software('Cufflinks', config['cufflinks']['path'])
    htseq = Software('HTSeq', config['htseq']['path'])

    # Step 1: Trimming | Cutadapt
    if step <= 1:
        for i, read in enumerate(reads):
            if run_is_paired_end:
                # Get paired-end reads, construct new filenames
                read1, read2 = read.split(',')
                trimmed_read1_filename = os.path.join(output_dir,
                                                      lib_prefix + '_{}_read1.trimmed.fastq.gz'.format(i))
                trimmed_read2_filename = os.path.join(output_dir,
                                                      lib_prefix + '_{}_read2.trimmed.fastq.gz'.format(i))

                # Run cutadapt
                cutadapt.run(
                    Parameter('--quality-base={}'.format(config['cutadapt']['quality-base'])),
                    Parameter('--minimum-length=5'),
                    Parameter('--output={}'.format(trimmed_read1_filename)),
                    Parameter('--paired-output={}'.format(trimmed_read2_filename)),
                    Parameter('-a', forward_adapter),
                    Parameter('-A', reverse_adapter),
                    Parameter('-q', '30'),
                    Parameter(read1),
                    Parameter(read2),
                    Redirect(type='1>', dest=os.path.join(logs_dir, 'cutadapt.summary'))
                )

                # Update reads list
                reads[i] = ','.join([trimmed_read1_filename, trimmed_read2_filename])
            else:
                # Construct new filename
                trimmed_read_filename = os.path.join(output_dir,
                                                     lib_prefix + '_{}.trimmed.fastq.gz'.format(i))

                # Run cutadapt
                cutadapt.run(
                    Parameter('--quality-base={}'.format(config['cutadapt']['quality-base'])),
                    Parameter('--minimum-length=5'),
                    Parameter('--output={}'.format(trimmed_read_filename)),
                    Parameter('-a', forward_adapter),
                    Parameter('-q', '30'),
                    Parameter(read),
                    Redirect(type='1>', dest=os.path.join(logs_dir, 'cutadapt.summary'))
                )

                # Update reads list
                reads[i] = trimmed_read_filename

    # Step 2: Alignment | STAR 2-pass
    if step <= 2:
        for i, read in enumerate(reads):
            if run_is_paired_end:
                read1, read2 = read.split(',')
                star.run(
                    Parameter('--runMode', 'alignReads'),
                    Parameter('--twopassMode', 'Basic'),
                    Parameter('--outFileNamePrefix', lib_prefix if lib_prefix[-1] == '.' else lib_prefix + '.'),
                    Parameter('--runThreadN', '8'),
                    Parameter('--genomeDir', config['STAR']['genome-dir']),
                    Parameter('--readFilesIn', read1, read2),
                    Parameter('--readFilesCommand', 'zcat'),
                    Parameter('--quantMode', 'TranscriptomeSAM', 'GeneCounts'),
                    Parameter('--outSAMtype', 'BAM', 'Unsorted'),
                    Parameter('--outFilterType', 'BySJout'),
                    Parameter('--outFilterMultimapNmax', '20'),
                    Parameter('--alignSJoverhangMin', '8'),
                    Parameter('--alignSJDBoverhangMin', '1'),
                    Parameter('--outFilterMismatchNmax', '8'),  # Should this be 2?
                    Parameter('--alignIntronMin', '20'),
                    Parameter('--alignIntronMax', '1000000'),
                    Parameter('--alignMatesGapMax', '1000000'),
                    Parameter('--outFilterIntronMotifs', 'RemoveNoncanonical') if run_is_stranded else Parameter()
                )
            else:
                # TODO Implement single-end STAR
                pass

    # Step 3a: Quantification | Cufflinks
    if step <= 3:
        pass

    # Step 3b: Quantification | HTSeq
