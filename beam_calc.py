#! /usr/bin/env python
#-*- coding:utf-8 -*-
#Written on 8/29/2017 by Sheng-Jun Lin

import numpy as np
from scipy.interpolate import interp1d
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
        quan_SI = c_SI / quan_SI  # becomes wavelength
    diameter_SI = len_unit(diameter)
    ratio_arcsec = quan_SI / diameter_SI / np.pi * 180 * 3600
    print('Primary Beam (FWHP) [1.02x] = {:.2f}" '.format(1.02 * ratio_arcsec))
    print('The 1st Null [1.22x]        = {:.2f}"'.format(1.22 * ratio_arcsec))
    print('')
    if diameter_SI == 12.:
        print('ALMA Primary Beam [1.13x] = {:.2f}"'.format(1.13 * ratio_arcsec))
        print('ALMA MRS ~ 0.5 * Primary beam')
    if diameter_SI == 30.:
        # Beam sizes of the Error beams
        EBs, EBs_err = IRAM30m_EBs(wavelen_SI=quan_SI)
        # Beam eff. of the Error beams
        b_eff, P1_pr, P2_pr, P3_pr = IRAM30m_eff(wavelen_SI=quan_SI)
        print('IRAM 30m MB (HPBW) [1.166x] = {0:.2f}"; B_eff = {1:.2f}'.format(1.166 * ratio_arcsec, b_eff))
        print('IRAM 30m EB1 (HPBW)         = {0:.2f} ({1:.2f})"; P1_pr = {2:.2f}'.format(EBs[0], EBs_err[0], P1_pr))
        print('IRAM 30m EB2 (HPBW)         = {0:.2f} ({1:.2f})"; P2_pr = {2:.2f}'.format(EBs[1], EBs_err[1], P2_pr))
        print('IRAM 30m EB3 (HPBW)         = {0:.2f} ({1:.2f})"; P3_pr = {2:.2f}'.format(EBs[2], EBs_err[2], P3_pr))

def IRAM30m_EBs(wavelen_SI):

    # theata["] * freq[GHz] = k (Equ. 2 in Kramer+2013)
    freq_SI = c_SI / wavelen_SI
    freq_GHz = freq_SI / 1e9
    # k values of the 1st, 2nd, and 3rd error beams
    k_EBs = np.array([13e3, 50e3, 175e3])
    k_EBs_err = np.array([1e3, 2e3, 3e3])
    th_EBs = k_EBs / freq_GHz
    th_EBs_err = k_EBs_err / freq_GHz
    return (th_EBs, th_EBs_err)

def IRAM30m_eff(wavelen_SI):

    freq_SI = c_SI / wavelen_SI
    freq_GHz = freq_SI / 1e9
    freq_measured = np.array([86., 115., 145., 210., 230., 280., 340., 345.]) # GHz
    # B, P1', P2', P3', F effs measured at freq_measured [GHz] (Table 2 in Kramer+2013)
    # B eff: Main beam eff.
    Beff_measured = np.array([.81, .78, .74, .63, .59, .49, .35, .34])
    # P1'_eff: The 1st error beam eff.
    P1pr_measured = np.array([0., .01, .02, .04, .04, .03, .02, .02])
    # P2'_eff: The 2nd error beam eff.
    P2pr_measured = np.array([.07, .08, .09, .11, .11, .11, .11, .11])
    # P3'_eff: The 3rd error beam eff.
    P3pr_measured = np.array([.06, .06, .06, .09, .11, .15, .14, .14])
    # F_eff: Forward eff. = B_eff + P1'_eff + P2'_eff + P3'_eff + fss_eff
    Feff_measured = np.array([.95, .94, .93, .94, .92, .87, .81, .80])
    # Note: 1.00 = F_eff + rss_eff (rearward spillover&scattering eff.)
    # These effs are normailized to a solid angle of 4pi.
    # i.e. eff = eta_gain(?) * int@Ω(Pn(Ω)dΩ) / int@Ω_4pi(Pn(Ω)dΩ)
    #          = G/4pi*int@Ω_4pi(Pn(Ω)dΩ) * int@Ω(Pn(Ω)dΩ) / int@Ω_4pi(Pn(Ω)dΩ)
    #          = G / 4pi * int@Ω(Pn(Ω)dΩ) (?)

    # Linear interpolation of these effs with frequencies
    Beff = interp1d(freq_measured, Beff_measured)
    P1pr = interp1d(freq_measured, P1pr_measured)
    P2pr = interp1d(freq_measured, P2pr_measured)
    P3pr = interp1d(freq_measured, P3pr_measured)
    Feff = interp1d(freq_measured, Feff_measured)

    # Efficiencies normailized to 4pi
    effs = [Beff(freq_GHz), P1pr(freq_GHz), P2pr(freq_GHz), P3pr(freq_GHz)]
    Feff_0 = np.asscalar(Feff(freq_GHz))

    # Efficiencies normailized to 2pi, *_eff/F_eff
    # Note (Kutner&Ulich 1981, Tools of Radio Astronomy: Ch7&8):
    # The antenna temperature, T_A^* (IRAM30m data),
    # is correscted for atm. attenuation, resistive losses,
    # and rearward spillover & scattering.
    # Or call it as forward beam brightness T.
    # i.e. T_A^* = T_A * exp(tau*A) / eta_gain / eta_rss
    #      eta_rss = (int@Ω_2pi/int@Ω_4pi)(Pn(Ω)dΩ)
    #      eta_gain = G/4pi * int@Ω_4pi(Pn(Ω)dΩ)
    #      eta_rss * eta_gain = F_eff for IRAM30m
    # The radiation temperature, T_R^*, is also corrected for
    # forward spillover & scattering.
    # i.e. T_R^* = T_A^* / eta_fss = T_MB for IRAM30m ...assuming T_P's = 0
    #      eta_fss = (int@Ω_diffraction/int@Ω_2pi)(Pn(Ω)dΩ)
    #              = B_eff / F_eff for IRAM30m [Ω_4pi] ...assuming P'_effs = 0
    #              = eta_MB for IRAM30m [Ω_2pi] ...assuming P'_effs = 0
    # =>   T_A^* = T_R^* (B_eff/F_eff + Sum[ P'_eff ]/F_eff)
    #               + (T[->0] * fss_eff/F_eff)?
    # T_R is the source radiation temperature.
    #      T_R = T_R^* / eta_c
    #      eta_c = int@Ω_s(Pn(Ψ-Ω)*B(Ψ)dΨ) / int@Ω_diffraction(Pn(Ω)dΩ)
    #            ~= (θ_s)^2 / (θ_s ^2 + θ_MB ^2)
    renormalized2Feff = [np.asscalar(e) / Feff_0 for e in effs]
    return tuple(renormalized2Feff)

def main():

    help_doc = 'beam_calc.py:\n\
Calculating beam sizes for Circular Fraunhofer (1st order approx.), and IRAM 30m.\n\
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
