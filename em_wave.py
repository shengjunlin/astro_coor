#! /usr/bin/env python
#-*- coding:utf-8 -*-
#Written on 8/29/2018 by Sheng-Jun Lin

from astro_coor import *
from sys import argv
from getopt import getopt

def main():

    help_doc = 'em_wave.py:\n\
Calculating frequency and wavelength of the input EM wave.\n\
Usage: em_wave.py freq/wavelength\n\
freq/wavelength need to be suffixed with their units.'

    opts, args = getopt(argv[1:], "h")

    for opt, value in opts:
        if opt == '-h':
            print(help_doc)
            exit()

    wave_unit(*args, screen=True)

if __name__ == '__main__':
    main()
