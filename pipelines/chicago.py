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

    try:
        run_is_stranded = True if options['extra_info']['run_is_stranded'].lower() == 'true' else False
        cufflinks_lib_type = options['extra_info']['cufflinks_lib_type']
        htseq_stranded = options['extra_info']['htseq_stranded']  # yes|no|reverse
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

    # Housekeeping
    star_output = []
    novosort_outfile = ''
    if htseq_stranded not in ['yes', 'no', 'reverse']:
        raise ValueError('htseq_stranded is not yes|no|reverse')
    if cufflinks_lib_type not in ['ff-firststrand',
                                  'ff-secondstrand',
                                  'ff-unstranded',
                                  'fr-firststrand',
                                  'fr-secondstrand',
                                  'fr-unstranded',
                                  'transfrags']:
        raise ValueError('cufflinks_lib_type is not ff-firststrand|ff-secondstrand|ff-unstranded|' +
                         'fr-firststrand|fr-secondstrand|fr-unstranded|transfrags')

    # Step 1: Trimming | Cutadapt
    if step <= 1:
        for i, read in enumerate(reads):
            if run_is_paired_end:
                # Get paired-end reads, construct new filenames
                read1, read2 = read.split(':')
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
                    Parameter('-a', forward_adapter if forward_adapter else 'ZZZ'),
                    Parameter('-A', reverse_adapter if reverse_adapter else 'ZZZ'),
                    Parameter('-q', '30'),
                    Parameter(read1),
                    Parameter(read2),
                    Redirect(type='1>', dest=os.path.join(logs_dir, 'cutadapt.chicago.summary'))
                )

                # Update reads list
                reads[i] = ':'.join([trimmed_read1_filename, trimmed_read2_filename])
            else:
                # Construct new filename
                trimmed_read_filename = os.path.join(output_dir,
                                                     lib_prefix + '_{}.trimmed.fastq.gz'.format(i))

                # Run cutadapt
                cutadapt.run(
                    Parameter('--quality-base={}'.format(config['cutadapt']['quality-base'])),
                    Parameter('--minimum-length=5'),
                    Parameter('--output={}'.format(trimmed_read_filename)),
                    Parameter('-a', forward_adapter if forward_adapter else 'ZZZ'),
                    Parameter('-q', '30'),
                    Parameter(read),
                    Redirect(type='1>', dest=os.path.join(logs_dir, 'cutadapt.chicago.summary'))
                )

                # Update reads list
                reads[i] = trimmed_read_filename

    # Step 2: Alignment | STAR 2-pass, Alignment Stats | samtools flagstat
    if step <= 2:
        for i, read in enumerate(reads):
            if run_is_paired_end:
                read1, read2 = read.split(':')
                star_outfile_prefix = os.path.join(output_dir,
                                                   lib_prefix + ('_' if lib_prefix[-1] != '.' else '') + '{}.')
                star.run(
                    Parameter('--runMode', 'alignReads'),
                    Parameter('--twopassMode', 'Basic'),
                    Parameter('--outFileNamePrefix', star_outfile_prefix.format(i)),
                    Parameter('--runThreadN', config['STAR']['threads']),
                    Parameter('--genomeDir', config['STAR']['genome-dir']),
                    Parameter('--readFilesIn', read1, read2),  # TODO Account for single-end
                    Parameter('--readFilesCommand', 'zcat'),
                    Parameter('--quantMode', 'TranscriptomeSAM', 'GeneCounts'),
                    Parameter('--outSAMtype', 'BAM', 'Unsorted'),
                    Parameter('--outFilterType', 'BySJout'),
                    Parameter('--outFilterMultimapNmax', '20'),
                    Parameter('--alignSJoverhangMin', '8'),
                    Parameter('--alignSJDBoverhangMin', '1'),
                    Parameter('--outFilterMismatchNmax', '2'),
                    Parameter('--alignIntronMin', '20'),
                    Parameter('--alignIntronMax', '1000000'),
                    Parameter('--alignMatesGapMax', '1000000'),
                    (
                        Parameter('--outFilterIntronMotifs', 'RemoveNoncanonical') if run_is_stranded
                        else Parameter('--outSAMstrandField', 'intronMotif')
                    )
                )

                samtools_flagstat.run(
                    Parameter(star_outfile_prefix + 'Aligned.out.bam'),
                    Redirect(type='>', dest=os.path.join(logs_dir, lib_prefix + '_{}.stat'.format(i)))
                )

                star_output.append(star_outfile_prefix.format(i) + 'Aligned.out.bam')
            else:
                # TODO Implement single-end STAR
                pass

    # Step 3: BAM Merge | Novosort
    if step <= 3:
        novosort_outfile = os.path.join(output_dir,
                                        lib_prefix + ('.' if lib_prefix[-1] != '.' else '') +
                                        'merged.Aligned.out.bam')
        novosort.run(
            Parameter('--tmpdir', os.path.join(output_dir, 'tmp')),
            Parameter(*[bam for bam in star_output]),
            Redirect(type='>', dest=novosort_outfile)
        )

    # Step 4a: Quantification | Cufflinks
    if step <= 4:
        cufflinks_output_dir = os.path.join(output_dir, 'cufflinks')
        subprocess.call(['mkdir', '-p', cufflinks_output_dir])
        cufflinks.run(
            Parameter('--GTF', config['cufflinks']['transcriptome-gtf']),
            Parameter('-p', config['cufflinks']['threads']),
            Parameter('--library-type', cufflinks_lib_type),
            Parameter('--upper-quartile-norm'),
            Parameter('-o', cufflinks_output_dir),
            Parameter('--max-bundle-frags', '1000000000'),
            Parameter(novosort_outfile)
        )

    # Step 4b: Quantification | HTSeq
    if step <= 5:
        htseq_output_dir = os.path.join(output_dir, 'htseq')
        subprocess.call(['mkdir', '-p', htseq_output_dir])
        for id_attr in ['gene_id', 'gene_name']:
            for feature_type in ['gene', 'transcript', 'exon']:
                htseq.run(
                    Parameter('-f', 'bam'),
                    Parameter('-r', 'name'),
                    Parameter('-s', htseq_stranded),
                    Parameter('-t', feature_type),
                    Parameter('-i', id_attr),
                    Parameter(novosort_outfile),
                    Parameter(config['htseq']['transcriptome-gtf']),
                    Redirect(type='>', dest=os.path.join(htseq_output_dir,
                                                         '{}.{}.counts'.format(feature_type,
                                                                               id_attr)))
                )
