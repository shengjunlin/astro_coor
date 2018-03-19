#! /usr/bin/env python
#-*- coding:utf-8 -*-
# Written on 6/14/2014 by Sheng-Jun Lin
# 8/15/2017: Rewrite the method to classify types of coor. forms

import numpy as np
from functools import total_ordering
import sys

deg2rad = np.pi / 180.
rad2deg = 180. / np.pi

@total_ordering
class ra(object):

    '''
    ra class stores the RA part of coordinates.
    RA forms: h:m:s / d:arcm:arcs / degree(all_d)
    Any one of forms is acceptable for creating a ra class
    since ra.converter() automatically calculates other forms.
        e.g. ra(h=4, m=31, s=30), ra(d=67, arcm=52, arcs=30),
             or ra(all_d=67.875)
    While print a ra class, only the hms form will be returned.
    ra classes are ordered,
        i.e. the operators: '<', '>', '==', and etc. are valid.
    '''
    def __init__(self, h=np.nan, m=np.nan, s=np.nan,
                 d=np.nan, arcm=np.nan, arcs=np.nan, all_d=np.nan):
        # priority = 1
        self.h = float(h)  # hour
        self.m = float(m)  # minute
        self.s = float(s)  # second
        # priority = 2
        self.d = float(d)  # deg
        self.arcm = float(arcm)  # arcmin
        self.arcs = float(arcs)  # arcsec
        # priority = 3
        self.all_d = float(all_d)  # all in deg
        self.priority = 0 # Initialize the state of priority.
        # Different numbers are corresponding to different system (hms/dms/all deg).

    def __eq__(self, other): # defined for = operator. Eg. ra(all_d=87) == ra(all_d=50)
        self.converter() # Make sure that the "all in deg" form exists.
        other.converter()
        return self.all_d == other.all_d

    def __lt__(self, other): # defined for < operator. Eg. ra(all_d=87) < ra(all_d=50)
        self.converter()
        other.converter()
        return self.all_d < other.all_d
    # Since "@total_ordering", it is not necessary to define > operator.
    # Once = and < are defined, then the comparison relations are well-defined.

    def __repr__(self): # defined for print(). Eg. print(ra(all_d=87))
        self.converter()
        return '{0}:{1}:{2}'.format(self.h, self.m, self.s)

    def set_h(self, h):
        self.h = int(h)  # hour
        self.priority = 1

    def set_m(self, m):
        self.m = int(m)  # minute
        self.priority = 1

    def set_s(self, s):
        self.s = float(s)  # second
        self.priority = 1

    def set_d(self, d):
        self.d = int(d)  # deg
        self.priority = 2

    def set_arcm(self, arcm):
        self.arcm = int(arcm)  # arcmin
        self.priority = 2

    def set_arcs(self, arcs):
        self.arcs = float(arcs)  # arcsec
        self.priority = 2

    def set_all_d(self, all_d):
        self.all_d = float(all_d)  # all in deg
        self.priority = 3

    def check_unit(self):
        if self.priority == 0:
            if (not np.isnan(self.h)) and (not np.isnan(self.m)) \
                    and (not np.isnan(self.s)):
                self.priority = 1  # hour-min-sec system
            elif (not np.isnan(self.d)) and (not np.isnan(self.arcm)) \
                    and (not np.isnan(self.arcs)):
                self.priority = 2  # deg-arcmin-arcsec system
            elif not np.isnan(self.all_d):
                self.priority = 3  # degree system
            else:
                raise ValueError('RA doesn\'t be assigned completely.')

    def hms2all_d(self):
        if self.h == np.nan or self.h < 0:
            raise ValueError("The hms form of RA isn't assigned completely!")
        all_d = self.h * 15.
        all_d += self.m * 0.25
        all_d += self.s * 0.25 / 60.
        self.set_all_d(all_d)

    def all_d2hms(self):
        if self.all_d == np.nan or self.all_d < 0:
            raise ValueError("The degree form of RA isn't assigned completely!")
        self.set_h(int(self.all_d / 15.))
        self.set_m(int((self.all_d - self.h * 15.) / 0.25))
        self.set_s((self.all_d - self.h * 15. - self.m * 0.25) / (0.25 / 60.))

    def dms2all_d(self):
        if self.d == np.nan or self.d < 0:
            raise ValueError("The dms form of RA isn't assigned completely!")
        all_d = self.d
        all_d += self.arcm / 60.
        all_d += self.arcs / 3600.
        self.set_all_d(all_d)

    def all_d2dms(self):
        if self.all_d == np.nan or self.all_d < 0:
            raise ValueError("The degree form isn't assigned completely!")
        self.set_d(int(self.all_d))
        self.set_arcm(int((self.all_d - self.d) * 60.))
        self.set_arcs((self.all_d - self.d) * 3600. - self.arcm * 60.)

    def converter(self):
        """
        Calculate all the forms of RA (hms/dms/all deg).
        """
        if self.priority == 0:
            # None of forms exists. Need to determine
            self.check_unit()
        if self.priority == 1:
            # The h:m:s form exists
            self.hms2all_d()
            self.all_d2dms()
            # If h and m in h:m:s are not integers, recalc the h:m:s form.
            h = self.h
            m = self.m
            if h % 1 != 0. or m % 1 != 0.:
                self.all_d2hms()
            else:
                self.set_h(int(h))
                self.set_m(int(m))
        elif self.priority == 2:
            # The d:m:s form exists
            self.dms2all_d()
            self.all_d2hms()
            # If d and m in d:m:s are not integers, recalc the d:m:s form.
            d = self.d
            arcm = self.arcm
            if d % 1 != 0. or arcm % 1 != 0.:
                self.all_d2dms()
            else:
                self.set_d(int(d))
                self.set_arcm(int(arcm))
        else:
            # The degree form exists
            self.all_d2hms()
            self.all_d2dms()


@total_ordering
class dec(object):

    '''
    dec class stores the Dec part of coordinates.
    Dec forms: d:arcm:arcs / degree(all_d)
    Any one of forms is acceptable for creating a dec class
    since dec.converter() automatically calculates other forms.
        e.g. dec(sign='+', d=18, arcm=12, arcs=30),
             or dec(all_d=18.208333)
    Note that d is always nonagtive.
    While print a dec class, only the dms form will be returned.
    dec classes are ordered,
        i.e. the operators: '<', '>', '==', and etc. are valid.
    '''
    def __init__(self, sign='+', d=np.nan, arcm=np.nan, arcs=np.nan,
                 all_d=np.nan):
        self.sign = sign  # positive/negative ('+'/'-')
        # priority = -2
        self.d = float(d)  # deg
        self.arcm = float(arcm)  # arcmin
        self.arcs = float(arcs)  # arcsec
        # priority = -1
        self.all_d = float(all_d)  # all in deg
        self.priority = 0 # Initialize the state of priority.
        # Different numbers are corresponding to different system (All deg/dms).

    def __eq__(self, other):
        self.converter()
        other.converter()
        return self.all_d == other.all_d

    def __lt__(self, other):
        self.converter()
        other.converter()
        return self.all_d < other.all_d

    def __repr__(self):
        self.converter()
        return '{0}{1}:{2}:{3}'.format(self.sign, self.d, self.arcm, self.arcs)

    def set_sign(self, pn):
        self.sign = pn  # positive/negative
        self.priority = -2

    def set_d(self, d):
        self.d = int(d)  # deg
        self.priority = -2

    def set_arcm(self, arcm):
        self.arcm = int(arcm)  # arcmin
        self.priority = -2

    def set_arcs(self, arcs):
        self.arcs = float(arcs)  # arcsec
        self.priority = -2

    def set_all_d(self, all_d):
        self.all_d = float(all_d)  # all in deg
        self.priority = -3

    def check_unit(self):
        if self.priority == 0:
            if (not np.isnan(self.d)) and (not np.isnan(self.arcm)) \
                    and (not np.isnan(self.arcs)):
                self.priority = -2  # deg-arcmin-arcsec system
            elif not np.isnan(self.all_d):
                self.priority = -3  # degree sytem
            else:
                raise ValueError('Dec doesn\'t be assigned completely.')

    def dms2all_d(self):
        if self.d == np.nan:
            raise ValueError("The dms form of Dec isn't assigned completely!")
        all_d = self.d
        all_d += self.arcm / 60.
        all_d += self.arcs / 3600.
        if self.sign == '-':
            all_d *= -1
        self.set_all_d(all_d)

    def all_d2dms(self):
        if self.all_d == np.nan:
            raise ValueError(
                    "The degree form of Dec isn't assigned completely!")
        if self.all_d < 0:
            self.sign = '-'
        else:
            self.sign = '+'
        abs_all_d = abs(self.all_d)
        self.set_d(int(abs_all_d))
        self.set_arcm(int((abs_all_d - self.d) * 60.))
        self.set_arcs((abs_all_d - self.d) * 3600. - self.arcm * 60.)

    def converter(self):
        """
        Calculate all the forms of Dec (dms/all deg).
        """
        if self.priority == 0:
            # None of forms exists. Need to determine
            self.check_unit()
        if self.priority == -2:
            # The d:m:s form exists
            self.dms2all_d()
            # If d and m in d:m:s are not integers, recalc the d:m:s form.
            d = self.d
            arcm = self.arcm
            if d % 1 != 0. or arcm % 1 != 0.:
                self.all_d2dms()
            else:
                self.set_d(int(d))
                self.set_arcm(int(arcm))
        else:
            # The degree form exists
            self.all_d2dms()


class Coor(object):

    '''
    A wrapper class contains a pair of ra and dec classes.
        e.g. Coor(input_ra, input_dec)
    input_ra/dec could be:
        1. strings of RA (h:m:s/d:m:s/degree) and Dec (d:m:s/degree)
        2. floats of RA and Dec in pure degree forms
        3. the ndarray elements return from astropy.wcs.all_pix2world()
        4. ra class/dec class objects
    While print a Coor class, 'h:m:s, d:m:s' will show up.
    '''
    def __init__(self, input_ra, input_dec):
        if isinstance(input_ra, basestring):
            self.RA = resolve_RA(input_ra)
        elif isinstance(input_ra, float):
            self.RA = ra(all_d=input_ra)
        elif isinstance(input_ra, np.ndarray) \
                and len(np.shape(input_ra)) == 0:
            # For astropy.wcs
            self.RA = ra(all_d=np.asscalar(input_ra))
        else:
            self.RA = input_ra
        if isinstance(input_dec, basestring):
            self.Dec = resolve_Dec(input_dec)
        elif isinstance(input_dec, float):
            self.Dec = dec(all_d=input_dec)
        elif isinstance(input_dec, np.ndarray) \
                and len(np.shape(input_dec)) == 0:
            # For astropy.wcs
            self.Dec = dec(all_d=np.asscalar(input_dec))
        else:
            self.Dec = input_dec
        self.RA.converter()
        self.Dec.converter()

    def converter(self):
        self.RA.converter()
        self.Dec.converter()

    def __repr__(self):
        return '{0}, {1}'.format(self.RA.__repr__(), self.Dec.__repr__())


def converter(co):

    '''
    Calculate all the representation of a Coor/ra/dec class.
    '''
    co.converter()


def resolve_RA(string):

    '''
    Analysis the string of RA in either h:m:s, 'd'd:m:s, or degree forms.
        e.g. '4:31:30', 'd67:52:30', or '67.875'.
    Note that the dms form should be prefixed by 'd'.
    '''
    RA = ra()
    try:
        if 'd' in string:
            tmp = string.split('d')
            tmp = tmp[1].split(':')
            RA.d = float(tmp[0])
            RA.arcm = float(tmp[1])
            RA.arcs = float(tmp[2])
        elif ':' in string:
            tmp = string.split(':')
            RA.h = float(tmp[0])
            RA.m = float(tmp[1])
            RA.s = float(tmp[2])
        else:
            RA.all_d = float(string)
        RA.converter()
    except:
        raise ValueError('Please check your RA forms!')
    return RA


def resolve_Dec(string):

    '''
    Analysis the string of Dec in either h:m:s, d:m:s, or degree forms.
        e.g. '+18:12:30', or '+18.208333'
    '''
    Dec = dec()
    try:
        if ':' in string:
            tmp = string.split(':')
            if '-' in tmp[0]:
                Dec.sign = '-'
            if '+' in tmp[0]:
                tmp[0].replace('+', '')
            Dec.d = abs(float(tmp[0]))
            Dec.arcm = float(tmp[1])
            Dec.arcs = float(tmp[2])
        else:
            Dec.all_d = float(string)
        Dec.converter()
    except:
        raise ValueError('Please check your Dec forms!\n'
                'If there are +/- signs, '
                'the sign and number can\'t be separated by space.')
    return Dec

# The marks for DistanceMode() and DisplacementMode():
#      C          N
#     / \         ^
#  aa/   \bb      |
#   /     \   E<--
#  /   cc  \
# B - - - - A
#(no. 2)    (no. 1)

def CosineLaw_side(aa, bb, C):

    '''
    Get cc which belongs to [0, pi]
    '''
    cos_cc = np.cos(aa) * np.cos(bb) + np.sin(aa) * np.sin(bb) * np.cos(C)
    return np.arccos(cos_cc)

def CosineLaw_angle(aa, bb, cc):

    '''
    Get C which belongs to [0, pi]
    '''
    if aa == 0. or bb == 0.:
        # aa = 0 => A = 0: point C == point B
        # or bb = 0 => B = 0: point C == point A
        return np.nan
    else:
        # aa != 0 and bb != 0
        cos_C = \
            (np.cos(cc) - np.cos(aa) * np.cos(bb)) / (np.sin(aa) * np.sin(bb))
        if -1. <= cos_C <= 1.:
            # Normal case: A = (0, pi) and B = (0, pi)
            return np.arccos(cos_C)
        else:
            # (1) cc = aa + bb <=> A = 0, B = 0; C = pi
            # (2) bb = aa + cc <=> A = 0, B = pi; C = 0
            # (3) aa = bb + cc <=> A = pi, B = 0; C = 0
            raise ValueError('Please check if (A=0 and B=0), (A=0 and B=pi), '
                             'or (A=pi and B=0).')

def DistanceMode(*args):

    '''
    DistanceMode(ra_p1, dec_p1, ra_p2, dec_p2)
     or DistanceMode(Coor_p1, Coor_p2),
     where ra_p#, dec_p#, and Coor_p# are ra/dec/Coor classes,
    calculate the angular distance between p1 and p2.

    Return a list:
    [PA of p2 wrt p1[deg], PA of p1 wrt p2[deg], dist(p1, p2)[deg]]
    '''
    if len(args) == 4:
        # args = ra(), dec(), ra(), dec()
        ra_1 = args[0]
        dec_1 = args[1]
        ra_2 = args[2]
        dec_2 = args[3]
    elif len(args) == 2:
        # args = Coor(), Coor()
        ra_1 = args[0].RA
        dec_1 = args[0].Dec
        ra_2 = args[1].RA
        dec_2 = args[1].Dec
    else:
        raise ValueError('Please check input for DistanceMode!')
    # Calculate all forms of RA and Dec
    for i in [ra_1, dec_1, ra_2, dec_2]:
        i.converter()

    # Transfer the orgin of Dec to be North pole
    # and obtain aa, bb, and C in rad.
    aa = (90. - dec_2.all_d) * deg2rad
    bb = (90. - dec_1.all_d) * deg2rad
    C = (ra_2.all_d - ra_1.all_d) * deg2rad
    # Obtain A, B, cc; all in rad
    cc = CosineLaw_side(aa, bb, C)
    if cc == 0.:
        # <=> aa = bb, C = 0
        # point A == point B
        A = np.nan
        B = np.nan
    else:
        if C == 0.:
            if aa > bb:
                A = np.pi
                B = 0.
            elif bb > aa:
                A = 0.
                B = np.pi
        else:
            A = CosineLaw_angle(bb, cc, aa)
            B = CosineLaw_angle(aa, cc, bb)
            # Determine the direction of PA
            A = A if C > 0. else -A
            B = B if C < 0. else -B
    return [A * rad2deg, B * rad2deg, cc * rad2deg]

def DisplacementMode(*args):

    '''
    DisplacementMode(ra_p1, dec_p1, x["], pa[deg])
     or DisplacementMode(Coor_p1, x["], pa[deg]),
     where ra_p1, dec_p1 are ra/dec classes,
    calculate the displaced point [ra_p2, dec_p2].

    Return:
    [ra_p2, dec_p2] or Coor(ra_p2, dec_p2)
    '''
    if len(args) == 4:
        # args = ra(), dec(), x, pa
        ra_p1 = args[0]
        dec_p1 = args[1]
    elif len(args) == 3:
        # args = Coor(), x, pa
        ra_p1 = args[0].RA
        dec_p1 = args[0].Dec
    else:
        raise ValueError('Please check input for DisplacementMode!')
    x = args[-2]
    pa = args[-1]
    # Calculate all forms of RA and Dec
    for i in [ra_p1, dec_p1]:
        i.converter()

    # Obtain A, bb, cc; all in rad
    A = pa * deg2rad
    bb = (90. - dec_p1.all_d) * deg2rad
    cc = x / 3600. * deg2rad
    # Obtain aa, C; all in rad
    if cc == 0.:
        aa = bb
        C = 0.
    else:
        aa = CosineLaw_side(bb, cc, A)
        if pa % 180. == 0.:
            C = 0.
        else:
            C = CosineLaw_angle(aa, bb, cc)
            towardE = True if np.sin(A) > 0. else False
            C = C if towardE else -C
    # Calculate all forms of RA and Dec of the result
    ra_p2 = ra()
    dec_p2 = dec()
    ra_p2.all_d = ra_p1.all_d + C * rad2deg
    ra_p2.converter()
    dec_p2.all_d = 90. - aa * rad2deg
    dec_p2.converter()
    if len(args) == 4:
        return [ra_p2, dec_p2]
    if len(args) == 3:
        return Coor(ra_p2, dec_p2)
