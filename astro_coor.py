#! /usr/bin/env python
#-*- coding:utf-8 -*-
# Written on 6/14/2014 by Sheng-Jun Lin
# 08/15/2017: Rewrite the method to classify types of coor. forms
# 07/30/2018: Use property methods.

import numpy as np
from functools import total_ordering

deg2rad = np.pi / 180.
rad2deg = 180. / np.pi

def nanint(value):

    if isinstance(value, basestring):
        return int(value)
    elif np.isnan(value):
        return value
    elif isinstance(value, (int, float)):
        return value
    else:
        raise ValueError('"{0}" is not an int.'.format(value))

def nanfloat(value):

    if isinstance(value, basestring):
        return float(value)
    elif np.isnan(value):
        return value
    elif isinstance(value, (float, int)):
        return value
    else:
        raise ValueError('"{0}" is not a float.'.format(value))

@total_ordering
class ra(object):

    '''
    ra class stores the RA part of coordinates.
    RA forms:
        h:m:s (hms) / d:arcm:arcs (dms) / all_d (degree) / all_rad (radian)
    Any one of forms is acceptable for creating a ra class,
    then ra.converter() would be automatically applied to calculate the other forms.
        e.g. ra(h=4, m=31, s=30),
             ra('4:31:30'),
             ra(d=67, arcm=52, arcs=30),
             ra('d67:52:30'), # [would be interpreted as the dms form]
             ra(all_d=67.875),
             ra('67.875'), # [would be interpreted as the degree form]
             or ra(all_rad=1.18).
    While print a ra class, only the hms form will be returned.
    To retrieve certain forms of strings, there are some methods:
        ra.hms(sep=':'), ra.dms(), ra.degree(), or ra.radian(),
        where ':' is the default of sep.
    ra classes are ordered,
        i.e. the operators: '<', '>', '==', and etc. are valid.
    '''
    def __init__(self, ra_str='',
                 h=np.nan, m=np.nan, s=np.nan,
                 d=np.nan, arcm=np.nan, arcs=np.nan,
                 all_d=np.nan, all_rad=np.nan):

        if ra_str:
            self.resolve_RA(ra_str)
        else:
            # type = 1
            self.h = h  # hour
            self.m = m  # minute
            self.s = s  # second

            # type = 2
            self.d = d  # deg
            self.arcm = arcm  # arcmin
            self.arcs = arcs  # arcsec

            # type = 3
            self.all_d = all_d  # all in deg

            # type = 4
            self.all_rad = all_rad  # all in radian

        self.type = 0 # Initialize the state of type. Check with ra.converter().

    def resolve_RA(self, string):

        '''
        Analysis the string of RA in either h:m:s, 'd'd:m:s, or degree forms.
            e.g. '4:31:30', 'd67:52:30', or '67.875'.
        Note that the dms form should be prefixed by 'd',
        and the radian form is not acceptable.
        '''
        try:
            if 'd' in string:
                tmp = string.split('d')
                tmp = tmp[1].split(':')
                self.d = tmp[0]
                self.arcm = tmp[1]
                self.arcs = tmp[2]
            elif ':' in string:
                tmp = string.split(':')
                self.h = tmp[0]
                self.m = tmp[1]
                self.s = tmp[2]
            else:
                self.all_d = string
            self.converter()
        except:
            raise ValueError('Please check your RA forms "{0}"!'.format(string))

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
        return self.hms()

    @property
    def h(self):
        return self._h
    @h.setter
    def h(self, h):
        self._h = nanint(h)  # hour
        self.type = 1

    @property
    def m(self):
        return self._m
    @m.setter
    def m(self, m):
        self._m = nanint(m)  # minute
        self.type = 1

    @property
    def s(self):
        return self._s
    @s.setter
    def s(self, s):
        self._s = nanfloat(s)  # second
        self.type = 1

    def hms(self, sep=':'):
        return '{0:0>2d}{3}{1:0>2d}{3}{2:0>5.2f}'.format(self.h, self.m, self.s, sep)

    @property
    def d(self):
        return self._d
    @d.setter
    def d(self, d):
        self._d = nanint(d)  # deg
        self.type = 2

    @property
    def arcm(self):
        return self._arcm
    @arcm.setter
    def arcm(self, arcm):
        self._arcm = nanint(arcm)  # arcmin
        self.type = 2

    @property
    def arcs(self):
        return self._arcs
    @arcs.setter
    def arcs(self, arcs):
        self._arcs = nanfloat(arcs)  # arcsec
        self.type = 2

    def dms(self, sep=':'):
        return '{0:0>3d}{3}{1:0>2d}{3}{2:0>5.2f}'.format(self.d, self.arcm, self.arcs, sep)

    @property
    def all_d(self):
        return self._all_d
    @all_d.setter
    def all_d(self, all_d):
        self._all_d = nanfloat(all_d)  # all in deg
        self.type = 3

    def degree(self):
        return '{0:.5f}'.format(self.all_d)

    @property
    def all_rad(self):
        return self._all_rad
    @all_rad.setter
    def all_rad(self, all_rad):
        self._all_rad = nanfloat(all_rad)  # all in rad
        self.type = 4

    def radian(self):
        return '{0:.7f}'.format(self.all_rad)

    def check_unit(self):

        if self.type == 0:
            if (not np.isnan(self.h)) and (not np.isnan(self.m)) \
                    and (not np.isnan(self.s)):
                self.type = 1  # hour-min-sec system
            elif (not np.isnan(self.d)) and (not np.isnan(self.arcm)) \
                    and (not np.isnan(self.arcs)):
                self.type = 2  # deg-arcmin-arcsec system
            elif not np.isnan(self.all_d):
                self.type = 3  # degree system
            elif not np.isnan(self.all_rad):
                self.type = 4  # radian system
            else:
                raise ValueError("RA doesn't be assigned completely.")

    def hms2all_d(self):

        if self.h == np.nan or self.h < 0:
            raise ValueError("The hms form of RA isn't assigned completely!")
        all_d = (self.h + self.m/60. + self.s/3600.) * 15
        self.all_d = all_d

    def all_d2hms(self):

        if self.all_d == np.nan or self.all_d < 0:
            raise ValueError("The degree form of RA isn't assigned completely!")
        all_hr = self.all_d/15.
        self.h = int(all_hr)
        self.m = int( (all_hr - self.h)*60. )
        self.s = (all_hr - self.h)*3600. - self.m*60.

    def dms2all_d(self):

        if self.d == np.nan or self.d < 0:
            raise ValueError("The dms form of RA isn't assigned completely!")
        all_d = self.d
        all_d += self.arcm / 60.
        all_d += self.arcs / 3600.
        self.all_d = all_d

    def all_d2dms(self):

        if self.all_d == np.nan or self.all_d < 0:
            raise ValueError("The degree form isn't assigned completely!")
        self.d = int(self.all_d)
        self.arcm = int((self.all_d - self.d)*60.)
        self.arcs = (self.all_d - self.d)*3600. - self.arcm*60.

    def all_d2all_rad(self):

        if self.all_d == np.nan or self.all_d < 0:
            raise ValueError("The degree form isn't assigned completely!")
        self.all_rad = self.all_d * deg2rad

    def all_rad2all_d(self):

        if self.all_rad == np.nan or self.all_rad < 0:
            raise ValueError("The radian form isn't assigned completely!")
        self.all_d = self.all_rad * rad2deg

    def converter(self):

        '''
        Calculate all the forms of RA (hms/dms/all deg/all rad).
        '''
        if self.type == 0:
            # None of forms exists. Need to determine
            self.check_unit()

        if self.type == 1:
            # The h:m:s form exists
            self.hms2all_d()
            self.all_d2dms()
            self.all_d2all_rad()
        elif self.type == 2:
            # The d:m:s form exists
            self.dms2all_d()
            self.all_d2hms()
            self.all_d2all_rad()
        elif self.type == 3:
            # The degree form exists
            self.all_d2hms()
            self.all_d2dms()
            self.all_d2all_rad()
        else:
            # The radian form exists
            self.all_rad2all_d()
            self.all_d2hms()
            self.all_d2dms()


@total_ordering
class dec(object):

    '''
    dec class stores the Dec part of coordinates.
    Dec forms:
        d:arcm:arcs (dms) / all_d (degree) / all_rad (radian)
    Any one of forms is acceptable for creating a dec class,
    then dec.converter() would be automatically applied to calculate the other forms.
        e.g. dec(sign='+', d=18, arcm=12, arcs=30),
             dec('+18:12:30'),
             dec(all_d=18.208333),
             dec('18.208333'), [would be interpreted as the degree form]
             or dec(all_rad=0.317795)
    Note that d (deg in the dms form) is always non-nagtive.
    While print a dec class, only the dms form will be returned.
    To retrieve certain forms of strings, there are some methods:
        dec.dms(sep=':'), dec.degree(), or dec.radian(),
        where ':' is the default of sep.
    dec classes are ordered,
        i.e. the operators: '<', '>', '==', and etc. are valid.
    '''
    def __init__(self, dec_str='',
                 sign='+', d=np.nan, arcm=np.nan, arcs=np.nan,
                 all_d=np.nan, all_rad=np.nan):

        if dec_str:
            self.resolve_Dec(dec_str)
        else:
            # type = -2
            self.sign = sign  # positive/negative ('+'/'-')
            self.d = d  # deg
            self.arcm = arcm  # arcmin
            self.arcs = arcs  # arcsec

            # type = -3
            self.all_d = all_d  # all in deg

            # type = -4
            self.all_rad = all_rad  # all in rad

        self.type = 0 # Initialize the state of type. Check with dec.converter().

    def resolve_Dec(self, string):

        '''
        Analysis the string of Dec in either h:m:s, d:m:s, or degree forms.
            e.g. '+18:12:30', or '+18.208333'
        Note that the radian form is not acceptable.
        '''
        try:
            if ':' in string:
                tmp = string.split(':')
                # self.d and self.sign should be given separately.
                self.sign = '+'
                if '-' in tmp[0]:
                    self.sign = '-'
                    tmp[0].replace('-', '')
                # self.d is always non-negative.
                self.d = tmp[0]
                self.arcm = tmp[1]
                self.arcs = tmp[2]
            else:
                self.all_d = string
            self.converter()
        except:
            raise ValueError('Please check your Dec forms "{0}"!\n'
                    'If there are +/- signs, '
                    'the sign and number can\'t be separated by space.'.format(string))

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
        return self.dms()

    @property
    def sign(self):
        return self._sign
    @sign.setter
    def sign(self, pn):
        self._sign = pn  # positive/negative
        self.type = -2

    @property
    def d(self):
        return self._d
    @d.setter
    def d(self, d):
        self._d = nanint(d)  # deg
        # self.d is always non-negative.
        # self.d and self.sign should be given separately.
        if self._d < 0.:
            self._sign = '-'
            self._d = abs(self._d)
        self.type = -2

    @property
    def arcm(self):
        return self._arcm
    @arcm.setter
    def arcm(self, arcm):
        self._arcm = nanint(arcm)  # arcmin
        self.type = -2

    @property
    def arcs(self):
        return self._arcs
    @arcs.setter
    def arcs(self, arcs):
        self._arcs = nanfloat(arcs)  # arcsec
        self.type = -2

    def dms(self, sep=':'):
        return '{0}{1:0>2d}{4}{2:0>2d}{4}{3:0>5.2f}'.format(self.sign, self.d, self.arcm, self.arcs, sep)

    @property
    def all_d(self):
        return self._all_d
    @all_d.setter
    def all_d(self, all_d):
        self._all_d = nanfloat(all_d)  # all in deg
        self.type = -3

    def degree(self):
        return '{0:+.6f}'.format(self.all_d)

    @property
    def all_rad(self):
        return self._all_rad
    @all_rad.setter
    def all_rad(self, all_rad):
        self._all_rad = nanfloat(all_rad)  # all in rad
        self.type = -4

    def radian(self):
        return '{0:+.8}'.format(self.all_rad)

    def check_unit(self):

        if self.type == 0:
            if (not np.isnan(self.d)) and (not np.isnan(self.arcm)) \
                    and (not np.isnan(self.arcs)):
                self.type = -2  # deg-arcmin-arcsec system
            elif not np.isnan(self.all_d):
                self.type = -3  # degree sytem
            elif not np.isnan(self.all_rad):
                self.type = -4  # radian sytem
            else:
                raise ValueError("Dec doesn't be assigned completely.")

    def dms2all_d(self):

        if self.d == np.nan:
            raise ValueError("The dms form of Dec isn't assigned completely!")
        all_d = self.d
        all_d += self.arcm / 60.
        all_d += self.arcs / 3600.
        # self.d is always positive.
        # self.d and self.sign should be given separately.
        if self.sign == '-':
            all_d *= -1
        self.all_d = all_d

    def all_d2dms(self):

        if self.all_d == np.nan:
            raise ValueError(
                    "The degree form of Dec isn't assigned completely!")
        if self.all_d < 0:
            self.sign = '-'
        else:
            self.sign = '+'
        abs_all_d = abs(self.all_d)
        self.d = int(abs_all_d)
        self.arcm = int((abs_all_d - self.d)*60.)
        self.arcs = (abs_all_d - self.d)*3600. - self.arcm*60.

    def all_rad2all_d(self):

        if self.all_rad == np.nan:
            raise ValueError("The radian form isn't assigned completely!")
        self.all_d = self.all_rad * rad2deg

    def all_d2all_rad(self):

        if self.all_d == np.nan:
            raise ValueError("The degree form isn't assigned completely!")
        self.all_rad = self.all_d * deg2rad

    def converter(self):

        """
        Calculate all the forms of Dec (dms/all deg).
        """
        if self.type == 0:
            # None of forms exists. Need to determine
            self.check_unit()

        if self.type == -2:
            # The d:m:s form exists
            self.dms2all_d()
            self.all_d2all_rad()
        elif self.type == -3:
            # The degree form exists
            self.all_d2dms()
            self.all_d2all_rad()
        else:
            # The radian form exists
            self.all_rad2all_d()
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
            self.RA = ra(ra_str=input_ra)
        elif isinstance(input_ra, float):
            self.RA = ra(all_d=input_ra)
        elif isinstance(input_ra, np.ndarray) \
                and len(np.shape(input_ra)) == 0:
            # For astropy.wcs
            self.RA = ra(all_d=np.asscalar(input_ra))
        else:
            self.RA = input_ra

        if isinstance(input_dec, basestring):
            self.Dec = dec(dec_str=input_dec)
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
    Note that the dms form should be prefixed by 'd',
    and the radian form is not acceptable.
    '''
    RA = ra(ra_str=string)
    return RA


def resolve_Dec(string):

    '''
    Analysis the string of Dec in either h:m:s, d:m:s, or degree forms.
        e.g. '+18:12:30', or '+18.208333'
    Note that the radian form is not acceptable.
    '''
    Dec = dec(dec_str=string)
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

