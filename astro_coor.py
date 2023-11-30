#! /usr/bin/env python
#-*- coding:utf-8 -*-
# #########################################################
# Author : Sheng-Jun Lin
# Email : sj.lin@gapp.nthu.edu.tw
# Description : Calculating the coordinate system and beam sizes
#    2014-JUL-14: Init. The astro_coor module has the following functions.
#        CosineLaw_side(), CosineLaw_angle(): Spherical cosine law.
#        DistanceMode(): Compute the great-circle distance.
#        DisplacementMode(): Compute a new location by an angular shift.
#                            This is a feature that astropy.SkyCoord (v2.x) doesn't include.
#    2017-AUG-15: Write the method to classify types of coord. forms.
#    2017-AUG-29: Add some functions about units, and calculating beam sizes (for calc_beam.py).
#    2018-JUL-30: Use property methods.
#    2018-AUG-29: Modify wave_unit a little bit (for em_wave.py).
#    *2019-FEB*   astropy v3.1 adds a new method: SkyCoord.directional_offset_by()
#                 that provides the same feature as DisplacementMode().
#    2020-MAR-20: Become compatible with Python3.
#                 Add xy_comp for DistanceMode()
#
# #########################################################
import sys
import numpy as np
from functools import total_ordering
from scipy.interpolate import interp1d
from astropy.coordinates import SkyCoord
from astropy import units as u


if (sys.version_info.major>=3):
    # Python3
    xrange = range
    izip = zip
else:
    # Python2
    from itertools import izip

c_SI = 299792458.

deg2rad = np.pi / 180.
rad2deg = 180. / np.pi

def nanint(value):

    if isinstance(value, str):
        return int(value)
    elif np.isnan(value):
        return value
    elif isinstance(value, (int, float)):
        return value
    else:
        raise ValueError('"{0}" is not an int.'.format(value))

def nanfloat(value):

    if isinstance(value, str):
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
        e.g. [hms]: ra(h=4, m=31, s=30.0),
                    ra('4:31:30.0'), ra('4h31m30.0s'), ra('4 31 30.0'),
             [dms]: ra(d=67, arcm=52, arcs=30.0),
                    ra('d67:52:30.0'), ra('d67 52 30.0'),
             [deg]: ra(all_d=67.875),
                    ra('67.875'),
             [rad]: ra(all_rad=1.18).
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
            e.g. [hms]: '4:31:30.0', '4h31m30.0s', '4H31M30.0S', '4 31 30.0',
                 [dms]: 'd67:52:30.0', 'd67 52 30.0',
                 [deg]: '67.875'.
        Note that the dms form should be prefixed by 'd',
        and the radian form is not acceptable.
        '''
        string0 = string
        try:
            string = string.lower().replace('h',' ').replace('m',' ').replace('s',' ')
            string = string.replace(':',' ')
            ls_str = string.split()
            if 'd' in string:
                tmp = string.split('d')
                tmp = tmp[1].split()
                self.d = tmp[0]
                self.arcm = tmp[1]
                self.arcs = tmp[2]
            elif len(ls_str) == 3:
                self.h = ls_str[0]
                self.m = ls_str[1]
                self.s = ls_str[2]
            else:
                self.all_d = string
            self.converter()
        except:
            raise ValueError('Please check your RA forms "{0}"!'.format(string0))

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
        return '{0:0>2d}{3}{1:0>2d}{3}{2:0>7.4f}'.format(self.h, self.m, self.s, sep)

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
        e.g. [dms]: dec(sign='+', d=18, arcm=12, arcs=30.0),
                    dec('+18:12:30.0'), dec('+18d12m30.0s'),
                    dec('+18d12\'30.0\"'), dec('+18 12 30.0'),
             [deg]: dec(all_d=18.208333),
                    dec('18.208333'),
             [rad]: dec(all_rad=0.317795).
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
            e.g. [dms]: '+18:12:30.0', '+18d12m30.0s', '+18D12M30.0S', '+18 12 30.0'
                 [deg]: '+18.208333'
        Note that the radian form is not acceptable.
        '''
        string0 = string
        try:
            string = string.lower().replace('d',' ').replace('m',' ').replace('s',' ')
            string = string.replace("'",' ').replace('"',' ')
            string = string.replace(':',' ')
            ls_str = string.split()
            if len(ls_str) == 3:
                # self.d and self.sign should be given separately.
                self.sign = '+'
                if '-' in ls_str[0]:
                    self.sign = '-'
                    ls_str[0].replace('-', '')
                # self.d is always non-negative.
                self.d = ls_str[0]
                self.arcm = ls_str[1]
                self.arcs = ls_str[2]
            else:
                self.all_d = string
            self.converter()
        except:
            raise ValueError('Please check your Dec forms "{0}"!\n'
                    'If there are +/- signs, '
                    'the sign and number can\'t be separated by space.'.format(string0))

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
        return '{0}{1:0>2d}{4}{2:0>2d}{4}{3:0>6.3f}'.format(self.sign, self.d, self.arcm, self.arcs, sep)

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
    def __init__(self, input_ra, input_dec, frame='fk5'):

        if isinstance(input_ra, str):
            self.RA = ra(ra_str=input_ra)
        elif isinstance(input_ra, float):
            self.RA = ra(all_d=input_ra)
        elif isinstance(input_ra, np.ndarray) \
                and len(np.shape(input_ra)) == 0:
            # For astropy.wcs
            self.RA = ra(all_d=np.asscalar(input_ra))
        else:
            self.RA = input_ra

        if isinstance(input_dec, str):
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
        self.frame = frame

    def converter(self):

        self.RA.converter()
        self.Dec.converter()

    def __repr__(self):

        return '{0}, {1}'.format(self.RA.__repr__(), self.Dec.__repr__())

    def hms(self, sep=':'):
        return self.RA.hms(sep)

    def dms(self, sep=':'):
        return self.Dec.dms(sep)

    def to_SkyCoord(self, frame='', **kwargs):
        if not frame:
            frame = self.frame
        return SkyCoord(self.hms(), self.dms(), unit=(u.hourangle, u.deg), **kwargs)


class Loc(Coor):

    """Location of a target.
    """

    def __init__(self, RA, Dec, name='', dist_pc=0., frame='fk5'):
        """
        Args:
          RA: str. Any format.
          Dec: str. Any format.
          name: str. Target name.
          dist_pc: float. Distance in parsec.
        """
        super(Loc, self).__init__(RA, Dec)
        self.name = name
        self.dist_pc = dist_pc
        # The rad/deg forms of RA and Dec
        self.RA_rad = self.RA.all_rad
        self.Dec_rad = self.Dec.all_rad
        self.RA_deg = self.RA.all_d
        self.Dec_deg = self.Dec.all_d
        self.frame = frame

    @property
    def pix_size_amin(self):
        """Set the pixel size in arcmin or retrieve the value.
        E.g.,
        >>> A = Loc("00:00:00.0", "+00:00:00.0")
        >>> A.pix_size_amin(0.5)  # Set the pixel size as 0.5 arcmin
        >>> print(A.pix_size_amin)  # Obtain the pixel size in arcmin
        >>> print(A.pix_size_asec)  # The pixel size in arcsec is also computed
        """
        return self._pix_size_amin

    @pix_size_amin.setter
    def pix_size_amin(self, pix_size_amin):
        self._pix_size_amin = pix_size_amin
        self._pix_size_asec = pix_size_amin * 60.

    @property
    def pix_size_asec(self):
        """Set the pixel size in arcsec or retrieve the value.
        E.g.,
        >>> A = Loc("00:00:00.0", "+00:00:00.0")
        >>> A.pix_size_asec(30.)  # Set the pixel size as 30 arcsec
        >>> print(A.pix_size_asec)  # Obtain the pixel size in arcsec
        >>> print(A.pix_size_amin)  # The pixel size in arcmin is also computed
        """
        return self._pix_size_asec

    @pix_size_asec.setter
    def pix_size_asec(self, pix_size_asec):
        self._pix_size_asec = pix_size_asec
        self._pix_size_amin = pix_size_asec / 60.

    @property
    def pix_size_au(self):
        """Once the pixel size is set by the methods, pix_size_asec
        or pix_size_amin, we can retrieve the pixel size in au.
        E.g.,
        >>> A = Loc("00:00:00.0", "+00:00:00.0")
        >>> A.pix_size_asec(30.)  # Set the pixel size as 30 arcsec
        >>> print(A.pix_size_au)  # Obtain the pixel size in au
        """
        self._pix_size_au = self.dist_pc * self.pix_size_asec
        return self._pix_size_au

    @property
    def pix_size_pc(self):
        """Once the pixel size is set by the methods, pix_size_asec
        or pix_size_amin, we can retrieve the pixel size in parsec.
        E.g.,
        >>> A = Loc("00:00:00.0", "+00:00:00.0")
        >>> A.pix_size_asec(30.)  # Set the pixel size as 30 arcsec
        >>> print(A.pix_size_pc)  # Obtain the pixel size in pc
        """
        self._pix_size_pc = self.pix_size_au * au_SI / pc_SI
        return self._pix_size_pc

    @property
    def pix_size_m(self):
        """Once the pixel size is set by the methods, pix_size_asec
        or pix_size_amin, we can retrieve the pixel size in meter.
        E.g.,
        >>> A = Loc("00:00:00.0", "+00:00:00.0")
        >>> A.pix_size_asec(30.)  # Set the pixel size as 30 arcsec
        >>> print(A.pix_size_m)  # Obtain the pixel size in meter
        """
        self._pix_size_m = self.pix_size_au * au_SI
        return self._pix_size_m

    @property
    def pix_size_cm(self):
        """Once the pixel size is set by the methods, pix_size_asec
        or pix_size_amin, we can retrieve the pixel size in centimeter.
        E.g.,
        >>> A = Loc("00:00:00.0", "+00:00:00.0")
        >>> A.pix_size_asec(30.)  # Set the pixel size as 30 arcsec
        >>> print(A.pix_size_cm)  # Obtain the pixel size in centimeter
        """
        self._pix_size_cm = self.pix_size_au * au_SI * 100.
        return self._pix_size_cm

    @property
    def asec_in_au(self):
        """Once self.dist_pc is set,
        we can retrieve the 1-arcsec projected scale in au.
        E.g.,
        >>> A = Loc("00:00:00.0", "+00:00:00.0", dist_pc=140.)
        >>> print(A.asec_in_au)  # Projected length of 1 arcsec in au.
        """
        return self.dist_pc

    @property
    def amin_in_au(self):
        """Once self.dist_pc is set,
        we can retrieve the 1-arcmin projected scale in au.
        E.g.,
        >>> A = Loc("00:00:00.0", "+00:00:00.0", dist_pc=140.)
        >>> print(A.amin_in_au)  # Projected length of 1 arcmin in au.
        """
        return self.dist_pc * 60.

    @property
    def au_in_asec(self):
        """Once self.dist_pc is set,
        we can retrieve the 1-au projected length in the units of arcsec.
        E.g.,
        >>> A = Loc("00:00:00.0", "+00:00:00.0", dist_pc=140.)
        >>> print(A.au_in_asec)  # Projected length of 1 au in arcsec.
        """
        return 1. / self.asec_in_au

    @property
    def au_in_amin(self):
        """Once self.dist_pc is set,
        we can retrieve the 1-au projected length in the units of arcmin.
        E.g.,
        >>> A = Loc("00:00:00.0", "+00:00:00.0", dist_pc=140.)
        >>> print(A.au_in_amin)  # Projected length of 1 au in arcmin.
        """
        return 1. / self.amin_in_au

    def to_SkyCoord(self, frame=''):
        """Return the Loc object to be an astropy.coordinates.SkyCoord object.

        Arg:
          frame: str. {'fk5', 'fk4'} Default value = ''.
        """
        if not frame:
            frame = self.frame
        return SkyCoord(self.hms(), self.dms(), frame=frame, unit=(u.hourangle, u.deg))


def SkyCoord_to_Coor(skycoord_obj):
    """Create a Coor class instance by the input astropy.SkyCoord object
    """
    return Coor(str(skycoord_obj.ra.to_string(unit=u.hourangle)),
                str(skycoord_obj.dec.to_string()),
                frame=skycoord_obj.frame.name)


def SkyCoord_to_Loc(skycoord_obj, name='', dist_pc=0.):
    """Create a Loc class instance by the input astropy.SkyCoord object

    2018-DEC: astropy.SkyCoord provides a new method, directional_offset_by(),
              in version 3.1.2, but astropy 2.0.16 is the last version that
              supports python 2.
    """
    return Loc(str(skycoord_obj.ra.to_string(unit=u.hourangle)),
               str(skycoord_obj.dec.to_string()),
               name, dist_pc, frame=skycoord_obj.frame.name)


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
# 1. C at the pole.
# 2. line(CA), line(CB) are longitude lines.
#
#      C = dRA                N
#     / \                     ^
#  aa/   \bb = pi/2 - Dec     |
#   /     \               E<--
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


def DistanceMode(*args, **kwargs):

    '''
    DistanceMode(ra_p1, dec_p1, ra_p2, dec_p2, xy_comp=False)
     or DistanceMode(Coor_p1, Coor_p2, xy_comp=False),
     where ra_p#, dec_p#, and Coor_p# are ra/dec/Coor classes,
    calculate the angular distance between p1 and p2.

    xy_comp: bool. Output dx["] and dy["]. (Default value=False)

    Return a tuple:
    if xy_comp=False:
      (dist(p1, p2)[deg], PA of p2 wrt p1[deg], PA of p1 wrt p2[deg])
    else:
      (dist[deg], PAof2wrt1[deg], PAof1wrt2[deg],
       dx of dist(p1,p2) wrt p1["],
       dy of dist(p1,p2) wrt p1["])
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

    xy_comp = False
    if kwargs:
        xy_comp = kwargs['xy_comp']

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
            else:  # bb == aa
                A = np.pi/2
                B = np.pi/2
        else:
            A = CosineLaw_angle(bb, cc, aa)
            B = CosineLaw_angle(aa, cc, bb)
            # Determine the direction of PA
            A = A if C > 0. else -A # A ∈ (-pi/2, +pi/2]
            B = B if C < 0. else -B # B ∈ (-pi/2, +pi/2]

    dist_deg = cc * rad2deg
    PAof2wrt1_deg = A * rad2deg
    PAof1wrt2_deg = B * rad2deg

    if not xy_comp:
        # Output:
        # (dist(p1, p2)[deg], PA of p2 wrt p1[deg], PA of p1 wrt p2[deg])
        return (dist_deg, PAof2wrt1_deg, PAof1wrt2_deg)
    else:
        # dRA (make dec_2 become dec_1)
        dx_deg, PA21_deg, _ = DistanceMode(ra_1, dec_1, ra_2, dec_1, xy_comp=False)
        dx_deg = -dx_deg if PA21_deg < 0. else dx_deg
        dx_as = dx_deg * 3600.
        # dDec
        dy_deg = dec_2.all_d - dec_1.all_d
        dy_as = dy_deg * 3600.

        # Output:
        # (dist[deg], PA2wrt1[deg], PA1wrt2[deg],
        #  dRA of dist(p2,p1) wrt p1["], dDec of dist(p1,p2) wrt p1["])
        return (dist_deg, PAof2wrt1_deg, PAof1wrt2_deg, dx_as, dy_as)


def DisplacementMode(*args, **kwargs):

    '''
    DisplacementMode(ra_p1, dec_p1, shift_as["], pa_deg[deg])
     or DisplacementMode(ra_p1, dec_p1, dRA["], dDec["], xy_comp=True)
     or DisplacementMode(Coor_p1, shift_as["], pa_deg[deg]),
     or DisplacementMode(Coor_p1, dRA["], dDec["], xy_comp=True),
     where ra_p1, dec_p1 are ra/dec classes,
    calculate the displaced point [ra_p2, dec_p2].

    xy_comp: bool. Input dx["] and dy["]. (Default value=False)

    Return:
    (ra_p2, dec_p2) or Coor(ra_p2, dec_p2)
    '''
    if len(args) == 4:
        # args = ra(), dec(), shift_as["], pa_deg[deg]
        # args = ra(), dec(), dRA["], dDec["]; kwargs = {'xy_comp': True}
        ra_p1 = args[0]
        dec_p1 = args[1]
    elif len(args) == 3:
        # args = Coor(), shift_as["], pa_deg[deg]
        # args = Coor(), dRA["], dDec["]; kwargs = {'xy_comp': True}
        ra_p1 = args[0].RA
        dec_p1 = args[0].Dec
    else:
        raise ValueError('Please check input for DisplacementMode!')

    xy_comp = False
    if kwargs:
        xy_comp = kwargs['xy_comp']

    # Calculate all forms of RA and Dec
    for i in [ra_p1, dec_p1]:
        i.converter()

    # Calculate shift_as["] and pa_deg[deg]
    if not xy_comp:
        shift_as = args[-2]
        pa_deg = args[-1]
    else:
        # Convert from dRA["] and dDec["]
        dx_as = args[-2]
        dy_as = args[-1]
        #dx = dRA * np.cos(dec_p1)
        pa_deg, shift_as = offsets_to_shift_PA(dx_as, dy_as)

    # Obtain A, bb, cc; all in rad
    A = pa_deg * deg2rad
    bb = (90. - dec_p1.all_d) * deg2rad
    cc = shift_as / 3600. * deg2rad
    # Obtain aa, C; all in rad
    if cc == 0.:
        aa = bb
        C = 0.
    else:
        aa = CosineLaw_side(bb, cc, A)
        if pa_deg % 180. == 0.:
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
        return (ra_p2, dec_p2)
    if len(args) == 3:
        return Coor(ra_p2, dec_p2)


def offsets_to_shift_PA(dx_as, dy_as):
    '''Convert offsets (dx_as, dy_as) to displacement (as) and P.A. (deg)

      Args:
        dx_as: float. Horizontial offset in arcsec.
        dy_as: float. Vertial offset in arcsec.
      Returns:
        (displacement_as, PA_deg)
    '''
    shift_as = np.sqrt(dx_as**2 + dy_as**2)
    pa_deg = np.arctan2(dx_as, dy_as) * rad2deg
    return (pa_deg, shift_as)


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


def wave_unit(quan, screen=False):

    if 'Hz' == quan[-2:]:
        quan_SI = freq_unit(quan)
        quan_type = 'freq'
        if screen:
            print('{0:f}GHz; {1:f}mm'.format(quan_SI/1e9, c_SI/quan_SI*1e3))
    elif 'm' == quan[-1]:
        quan_SI = len_unit(quan)
        quan_type = 'wlen'
        if screen:
            print('{0:f}GHz; {1:f}mm'.format(c_SI/quan_SI/1e9, quan_SI*1e3))
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
    if diameter_SI == 6.:
        print('SMA !')
    elif diameter_SI == 12. or diameter_SI == 7.:
        print('ALMA Primary Beam [1.13x] = {:.2f}"'.format(1.13 * ratio_arcsec))
        print('ALMA MRS ~ 0.5 * Primary beam')
    elif diameter_SI == 30.:
        # Beam sizes of the Error beams
        EBs, EBs_err = IRAM30m_EBs(wavelen_SI=quan_SI)
        # Beam eff. of the Error beams
        b_eff, P1_pr, P2_pr, P3_pr, Feff = IRAM30m_eff(wavelen_SI=quan_SI)
        print('IRAM 30m MB (HPBW) [1.166x] = {0:.2f}"; B_eff = {1:.2f}'.format(1.166 * ratio_arcsec, b_eff))
        print('IRAM 30m EB1 (HPBW)         = {0:.2f} ({1:.2f})"; P1_pr = {2:.2f}'.format(EBs[0], EBs_err[0], P1_pr))
        print('IRAM 30m EB2 (HPBW)         = {0:.2f} ({1:.2f})"; P2_pr = {2:.2f}'.format(EBs[1], EBs_err[1], P2_pr))
        print('IRAM 30m EB3 (HPBW)         = {0:.2f} ({1:.2f})"; P3_pr = {2:.2f}'.format(EBs[2], EBs_err[2], P3_pr))
        print('IRAM 30m Forward eff. = {0:.2f}'.format(Feff))
    elif diameter_SI == 100.:
        fwhm_arcsec, diffr_coef = GBT_fwhm(wavelen_SI=quan_SI)
        print('GBT 100m resolution (FWHM) [{0}x] = {1:.2f}"'.format(diffr_coef, fwhm_arcsec))


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
    Feff_freq = np.asscalar(Feff(freq_GHz))

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
    renormalized2Feff = [np.asscalar(e) / Feff_freq for e in effs]
    return tuple(renormalized2Feff + [Feff_freq])

def GBT_fwhm(wavelen_SI):

    freq_SI = c_SI / wavelen_SI
    freq_GHz = freq_SI / 1e9
    if freq_GHz > 1.:
        Te_Db = 14. # The egde taper is Te = 14+-2 Db for the Gregorian feed.
        # A formula: FWHM_arcsec = 747.6"~763.8" / freq_GHz (PG for the GBT, 2017)
    else:
        Te_Db = 18. # The edge taper is Te = 18+-2 Db for the prime focus receivers.
        # A formula: FWHM_arcsec = 763.8"~797.4" / freq_GHz (PG for the GBT, 2017)
    diffr_coef = (1.02 + 0.0135 * Te_Db)
    FWHM_rad = diffr_coef * wavelen_SI/100.
    FWHM_arcsec = FWHM_rad / np.pi * 180. * 3600.
    return (FWHM_arcsec, diffr_coef)
