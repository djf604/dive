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

    # Housekeeping
    star_output = []

    # Keep list of items to delete
    staging_delete = ['tmp']

    # Get extra options
    try:
        forward_adapter = options['extra_info']['forward_adapter']
        reverse_adapter = options['extra_info']['reverse_adapter']
        run_is_stranded = True if options['extra_info']['run_is_stranded'].lower() == 'true' else False
    except KeyError, e:
        # TODO Make better exception message
        raise KeyError('Some needed extra_info not given.')

    # Establish software instances
    cat = Software('cat', '/bin/cat')
    cutadapt = Software('cutadapt', config['cutadapt']['path'])
    star = Software('STAR', config['STAR']['path'])
    rsem = Software('RSEM', config['RSEM']['path'])
    bedGraph_to_bw = Software('bedGraphToBigWig', config['bedGraphToBigWig']['path'])

    if step <= 1 and len(reads) >= 2:
        if run_is_paired_end:
            # Aggregate read1s and read2s
            read1s, read2s = [], []
            for reads_set in reads:
                read1, read2 = reads_set.split(':')
                read1s.append(read1)
                read2s.append(read2)

            # Combine reads groups
            combined_reads = []
            for name, reads_group in [('read1', read1s), ('read2', read2s)]:
                combined_read_filename = os.path.join(output_dir, '{}.combined.{}.fastq.gz'.format(lib_prefix, name))
                combined_reads.append(combined_read_filename)
                staging_delete.append(combined_read_filename)
                cat.run(
                    Parameter(*[read for read in reads_group]),
                    Redirect(type='1>', dest=combined_read_filename)
                )

            # Update reads list
            reads = [':'.join(combined_reads)]
        else:
            # Combine reads
            combined_read_filename = os.path.join(output_dir, '{}.combined.fastq.gz'.format(lib_prefix))
            staging_delete.append(combined_read_filename)
            cat.run(
                Parameter(*[read for read in reads]),
                Redirect(type='1>', dest=combined_read_filename)
            )

            # Update reads list
            reads = [combined_read_filename]

    # Trim adapters with cutadapt
    if step <= 2:
        reads_set = reads[0]
        if run_is_paired_end:
            # Get paired-end reads, construct new filenames
            read1, read2 = reads_set.split(':')
            trimmed_read1_filename = os.path.join(output_dir, lib_prefix + '_read1.trimmed.fastq.gz')
            trimmed_read2_filename = os.path.join(output_dir, lib_prefix + '_read2.trimmed.fastq.gz')

            staging_delete.append(trimmed_read1_filename)
            staging_delete.append(trimmed_read2_filename)

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
                Redirect(type='1>', dest=os.path.join(output_dir, 'logs', 'cutadapt.trendy.summary'))
            )

            # Update reads list
            reads = ':'.join([trimmed_read1_filename, trimmed_read2_filename])

        else:
            # Construct new filename
            trimmed_read_filename = os.path.join(output_dir, lib_prefix + '.trimmed.fastq.gz')

            staging_delete.append(trimmed_read_filename)

            # Run cutadapt
            cutadapt.run(
                Parameter('--quality-base={}'.format(config['cutadapt']['quality-base'])),
                Parameter('--minimum-length=5'),
                Parameter('--output={}'.format(trimmed_read_filename)),
                Parameter('-a', forward_adapter if forward_adapter else 'ZZZ'),
                Parameter('-q', '30'),
                Parameter(reads[0]),
                Redirect(type='1>', dest=os.path.join(output_dir, 'logs', 'cutadapt.trendy.summary'))
            )

            # Update reads list
            reads = [trimmed_read_filename]

    # Step 3: Alignment
    if step <= 3:
        # Gets reads for paired-end and single-end
        if run_is_paired_end:
            read1, read2 = reads.split(':')
        else:
            read1 = reads[0]
            read2 = ''

        # Set up STAR parameters
        star_outfile_prefix = os.path.join(output_dir,
                                           lib_prefix + ('.' if lib_prefix[-1] != '.' else ''))
        star_common = [
            Parameter('--outFileNamePrefix', star_outfile_prefix),
            Parameter('--genomeDir', config['STAR']['genome-dir']),
            Parameter('--readFilesIn', read1, read2),
            Parameter('--readFilesCommand', 'zcat'),
            Parameter('--outFilterType', 'BySJout'),
            Parameter('--outFilterMultimapNmax', '20'),
            Parameter('--alignSJoverhangMin', '8'),
            Parameter('--alignSJDBoverhangMin', '1'),
            Parameter('--outFilterMismatchNmax', '999'),
            Parameter('--alignIntronMin', '20'),
            Parameter('--alignIntronMax', '1000000'),
            Parameter('--alignMatesGapMax', '1000000'),
            Parameter('--outSAMunmapped', 'Within'),
            Parameter('--outSAMattributes', 'NH', 'HI', 'AS', 'NM', 'MD'),
            Parameter('--outFilterMismatchNoverReadLmax', '0.04'),
            Parameter('--sjdbScore', '1')
        ]
        star_run = [
            Parameter('--runThreadN', config['STAR']['threads']),
            Parameter('--genomeLoad', 'LoadAndKeep'),
            Parameter('--limitBAMsortRAM', '10000000000')
        ]
        star_bam = [
            Parameter('--outSAMtype', 'BAM', 'SortedByCoordinate'),
            Parameter('--quantMode', 'TranscriptomeSAM')
        ]
        star_strand, star_wig = [], []

        if run_is_stranded:
            star_wig.append(Parameter('--outWigStrand', 'Stranded'))
        else:
            star_strand.append(Parameter('--outSAMstrandField', 'intronMotif'))
            star_wig.append(Parameter('--outWigStrand', 'Unstranded'))

        # TODO Encode has SAM Header metadata here, but I'm going to skip it for now
        star_meta = []

        # Run STAR alignment step
        star.run(*(star_common + star_run + star_bam + star_strand + star_meta))

        # Store STAR output files
        star_output_bam = star_outfile_prefix + 'Aligned.sortedByCoord.out.bam'

        # Generate bedGraph
        subprocess.call(['mkdir', 'p', os.path.join(output_dir, 'signal')])
        signal_output_dir = os.path.join(output_dir, 'signal',
                                         lib_prefix + ('.' if lib_prefix[-1] != '.' else ''))

        # Run STAR for signal generation
        star.run(
            Parameter('--runMode', 'inputAlignmentsFromBAM'),
            Parameter('--inputBAMfile', star_output_bam),
            Parameter('--outWigtype', 'bedGraph'),
            Parameter('--outFileNamePrefix', signal_output_dir),
            Parameter('--outWigReferencesPrefix', 'chr'),
            *star_wig
        )

        # Convert bedGraph to bigWig
        chrNL_txt = os.path.join(output_dir, 'chrNL.txt')
        with open(chrNL_txt, 'w') as chrNL_filehandle:
            subprocess.call(['grep', '^chr', os.path.join(config['STAR']['genome-dir'], 'chrNameLength.txt')],
                            stdout=chrNL_filehandle)

        # Generate temporary signal file path
        sig_tmp = os.path.join(output_dir, 'sig.tmp')
        staging_delete.append(sig_tmp)
        if run_is_stranded:
            strand = [None, '-', '+']
            for i_strand in [1, 2]:
                for i_mult in ['Unique', 'UniqueMultiple']:
                    # Get signal file for this iteration
                    signal_file = 'Signal.{}.str{}.out.bg'.format(i_mult, str(i_strand))
                    # Write to temporary signal file
                    with open(sig_tmp, 'w') as sig_tmp_filehandle:
                        subprocess.call(['grep', '^chr', os.path.join(signal_output_dir, signal_file)],
                                        stdout=sig_tmp_filehandle)
                    # Run bedGraph to bigWig conversion
                    bedGraph_to_bw.run(
                        Parameter(sig_tmp),
                        Parameter(chrNL_txt),
                        Parameter(os.path.join(signal_output_dir,
                                               'Signal.{}.strand{}.bw'.format(i_mult, strand[i_strand])))
                    )
        else:
            for i_mult in ['Unique', 'UniqueMultiple']:
                # Get signal file for this iteration
                signal_file = 'Signal.{}.str1.out.bg'.format(i_mult)
                # Write to temporary signal file
                with open(sig_tmp, 'w') as sig_tmp_filehandle:
                    subprocess.call(['grep', '^chr', os.path.join(signal_output_dir, signal_file)],
                                    stdout=sig_tmp_filehandle)
                # Run bedGraph to bigWig conversion
                bedGraph_to_bw.run(
                    Parameter(sig_tmp),
                    Parameter(chrNL_txt),
                    Parameter(os.path.join(signal_output_dir,
                                           'Signal.{}.unstranded.bw'.format(i_mult)))
                )

    # Sort transcriptome BAM to ensure order of reads to make RSEM output deterministic
    if step >= 4:
        # Set BAM file paths, mv transcriptome BAM to temporary name
        star_outfile_prefix = os.path.join(output_dir,
                                           lib_prefix + ('.' if lib_prefix[-1] != '.' else ''))
        transcriptome_bam = star_outfile_prefix + 'Aligned.toTranscriptome.out.bam'
        tr_bam = star_outfile_prefix + 'Tr.bam'
        staging_delete.append(tr_bam)
        subprocess.call(['mv', transcriptome_bam, tr_bam])

        # Template command
        merge_cmd = 'cat <({input1}) <({input2}) | {compress} > {output}'
        input1_cmd = '{samtools} view -H {bam}'
        compress_cmd = 'samtools view -@ {threads} -bS -'

        if run_is_paired_end:
            input2_cmd = ('{samtools} view -@ {threads} {bam} | ' +
                          'awk \'{{printf "%s", $0 ""; getline; print}}\' | ' +
                          'sort -S {ram} -T {tmpdir} | ' +
                          'tr \' \' \'\\n\'')
        else:
            input2_cmd = ('{samtools} view -@ {threads} {bam} | ' +
                          'sort -S {ram} -T {tmpdir}')

        subprocess.call(merge_cmd.format(
            input1=input1_cmd.format(
                samtools=config['samtools']['path'],
                bam=tr_bam
            ),
            input2=input2_cmd.format(
                samtools=config['samtools']['path'],
                threads=config['RSEM']['threads'],
                bam=tr_bam,
                ram='32G',
                tmpdir=os.path.join(output_dir, 'tmp')
            ),
            compress=compress_cmd.format(
                threads=config['RSEM']['threads']
            ),
            output=transcriptome_bam
        ), shell=True, executable='/bin/bash')

