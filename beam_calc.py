#! /usr/bin/env python
#-*- coding:utf-8 -*-
#Written on 8/29/2017 by Sheng-Jun Lin

from astro_coor import *
from sys import argv
from getopt import getopt

def main():

    help_doc = 'beam_calc.py:\n\
Calculating beam sizes for Circular Fraunhofer (1st order approx.), and IRAM 30m.\n\
Usage: beam_calc.py freq/wavelength diameter\n\
E.g. ./beam_calc.py 350GHz 12m\n\
freq/wavelength/diameter need to be suffixed with their units.'

    opts, args = getopt(argv[1:], "h")

    for opt, value in opts:
        if opt == '-h':
            print(help_doc)
            exit()

    beam_size(*args)

if __name__ == '__main__':
    main()
