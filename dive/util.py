__author__ = 'Dominic Fitzgerald'
import os


def clean_up(*args):
    for filename in args:
        os.remove(filename)
