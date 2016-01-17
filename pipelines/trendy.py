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

    # Keep list of items to delete
    staging_delete = ['tmp']

    try:
        forward_adapter = options['extra_info']['forward_adapter']
        reverse_adapter = options['extra_info']['reverse_adapter']
        sailfish_libtype = options['extra_info']['sailfish_libtype']
    except KeyError, e:
        # TODO Make better exception message
        raise KeyError('Some needed extra_info not given.')

    # Establish Software instances
    cat = Software('cat', '/bin/cat')
    cutadapt = Software('cutadapt', config['cutadapt']['path'])
    kallisto = Software('kallisto', config['kallisto']['path'])
    sailfish = Software('sailfish', config['sailfish']['path'])

    # Combine reads with extra sequencing depth
    if step <= 1 and len(reads) >= 2:
        if run_is_paired_end:
            # Aggregate read1s and read2s
            read1s, read2s = [], []
            for read in reads:
                read1, read2 = read.split(':')
                read1s.append(read1)
                read2s.append(read2)

            # Combine reads groups
            combined_reads = []
            for name, reads_group in [('read1', read1s), ('read2', read2s)]:
                combined_read_filename = os.path.join(output_dir, '{}.combined.{}.fastq.gz'.format(lib_prefix, name))
                combined_reads.append(combined_read_filename)
                cat.run(
                    Parameter(*[read for read in reads_group]),
                    Redirect(type='1>', dest=combined_read_filename)
                )

            # Update reads list
            reads = [','.join(combined_reads)]
        else:
            # Combine reads
            combined_read_filename = os.path.join(output_dir, '{}.combined.fastq.gz'.format(lib_prefix))
            cat.run(
                Parameter(*[read for read in reads]),
                Redirect(type='1>', dest=combined_read_filename)
            )

            # Update reads list
            reads = [combined_read_filename]

    # Trim adapters with cutadapt
    if step <= 2:
        reads = reads[0]
        if run_is_paired_end:
            # Get paired-end reads, construct new filenames
            read1, read2 = reads.split(':')
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

    # Run Kallisto
    if step <= 3:
        if run_is_paired_end:
            read1, read2 = reads.split(':')
            kallisto.run(
                Parameter('--index={}'.format(config['kallisto']['index-path'])),
                Parameter('--output-dir={}'.format(lib_prefix + '_kallisto_quant')),
                Parameter(read1),
                Parameter(read2)
            )
        else:
            kallisto.run(
                Parameter('--index={}'.format(config['kallisto']['index-path'])),
                Parameter('--output-dir={}'.format(lib_prefix + '_kallisto_quant')),
                Parameter(reads[0])
            )

    # Run Sailfish
    if step <= 4:
        if run_is_paired_end:
            read1, read2 = reads.split(':')
            sailfish.run(
                Parameter('--index', config['sailfish']['index-path']),
                Parameter('--libType', '\"{}\"'.format(sailfish_libtype)),
                Parameter('-1', '<(zcat {})'.format(read1)),
                Parameter('-2', '<(zcat {})'.format(read2)),
                Parameter('--output', lib_prefix + '_sailfish_quant')
            )
        else:
            sailfish.run(
                Parameter('--index', config['sailfish']['index-path']),
                Parameter('--libType', '\"{}\"'.format(sailfish_libtype)),
                Parameter('-r', '<(zcat {})'.format(reads[0])),
                Parameter('--output', lib_prefix + '_sailfish_quant')
            )

        # Delete staged items
        for item in staging_delete:
            subprocess.call(['rm', item])
