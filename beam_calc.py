#! /usr/bin/env python
#-*- coding:utf-8 -*-
#Written on 8/29/2017 by Sheng-Jun Lin

import numpy as np
from sys import argv
from getopt import getopt

c_SI = 299792458.

def demical(quan):

    '''
    Convert a string with demical suffix to a float in SI unit
    '''
    demdict = {'G': 1e9, 'M': 1e6, 'k': 1e3,
               'c': 1e-2, 'm': 1e-3, 'u': 1e-6, 'n': 1e-9}
    nonSI = True
    for prefix in demdict:
        if prefix == quan[-1]:
            quan = float(quan[:-1])
            factor = demdict[prefix]
            nonSI = False
            return quan * factor
    if nonSI:
        return float(quan)

def len_unit(quan):

    '''
    Convert a string to a float in SI unit
    '''
    if 'm' == quan[-1]:
        quan_SI = demical(quan[:-1])
    else:
        raise ValueError('Counld not recoginze the length units.')
    return quan_SI

def freq_unit(quan):

    '''
    Convert a string to a float in SI unit
    '''
    if 'Hz' == quan[-2:]:
        quan_SI = demical(quan[:-2])
    else:
        raise ValueError('Counld not recoginze the frequency units.')
    return quan_SI

def wave_unit(quan):

    if 'Hz' == quan[-2:]:
        quan_SI = freq_unit(quan)
        quan_type = 'freq'
    elif 'm' == quan[-1]:
        quan_SI = len_unit(quan)
        quan_type = 'wlen'
    else:
        raise ValueError('Counld not recoginze the units.')
    return (quan_SI, quan_type)

def beam_size(quan, diameter):

    quan_SI, quan_type = wave_unit(quan)
    if quan_type == 'freq':
        quan_SI = c_SI / quan_SI
    diameter_SI = len_unit(diameter)
    ratio_arcsec = quan_SI / diameter_SI / np.pi * 180 * 3600
    print('Primary Beam (FWHP) [1.02x] = {:.2f} arcsec'.format(1.02 * ratio_arcsec))
    print('The 1st Null [1.22x]        = {:.2f} arcsec'.format(1.22 * ratio_arcsec))

def main():

    help_doc = 'beam_calc.py:\n\
Usage: beam_calc.py freq/wavelength diameter\n\
freq/wavelength/diameter need to be suffixed with their units.'

    opts, args = getopt(argv[1:], "h")

    for opt, value in opts:
        if opt == '-h':
            print(help_doc)
            exit()

    beam_size(*args)

if __name__ == '__main__':
    main()
