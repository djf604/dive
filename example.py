__author__ = 'Dominic Fitzgerald'
import sys
from dive.components import Software, Parameter, Redirect

sample = sys.argv[1]

# Instantiate software
# Software(software_name, software_path)
picard = Software('picard', '/path/to/java /path/to/picard.jar')

# Run software
# Put as many Parameter and Redirect as needed
# Order matters, so generally Redirect should be last

# Parameter('arguments', 'separated', 'by', 'spaces')
picard.run(
    Parameter('I=' + sample),
    Parameter('O=/path/to/output'),
    Parameter('-T', 'SplitNCigarReads'),
    Redirect(type='>', dest='out.txt')
)

# Will produce and execute:
# /path/to/java /path/to/picard.jar I=/path/to/input.bam O=/path/to/output -T SplitNCigarReads > out.txt
