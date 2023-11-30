# #########################################################
# Author : Sheng-Jun Lin
# Email : shengjunlin@asiaa.sinica.edu.tw
# Description : astro_coor
# #########################################################

"""
# Description :

A Python module to calculate all the RA/Dec forms, the angular distance,
the angular displaccement of the coordinates. There are also additional
functions handling the beam size calculation.

This module was written because the feature of DisplacementMode()
was not implemented by astropy 2.x. While, astropy 3.1 starts to provide such
a feature with SkyCoord.directional_offset_by(), and then the astro_coor.Coor
class can be almost replaced with astropy.coordinates.SkyCoord now.
However, the command-line interface is still useful.
The commands are the following.

**RADec_convert.py**

Features:
1. Convert given RA and Dec into different forms.
2. Calculate the distance and P.A. between two given points.
3. Calculate the shifted location from a given point with the displacement.

See more by 'python RADec_convert.py -h'

**beam_calc.py**

A python script calculating the beam size from given freq/wavelength and
dish diameter.
See more by 'python beam_calc.py -h'

**em_wave.py**

A python script calculating the freq from a given wavelength and vice versa.
See more by 'python em_wave.py -h'


# History:

2014-JUL-14: Init. The astro_coor module has the following functions.
    CosineLaw_side(), CosineLaw_angle(): Spherical cosine law.
    DistanceMode(): Compute the great-circle distance.
    DisplacementMode(): Compute a new location by an angular shift.
                        This is a feature that astropy.SkyCoord (v2.x) doesn't include.
2017-AUG-15: Write the method to classify types of coord. forms.
2017-AUG-29: Add some functions about units, and calculating beam sizes (for calc_beam.py).
2018-JUL-30: Use property methods.
2018-AUG-29: Modify wave_units a little bit (for em_wave.py).
*2019-FEB*   astropy v3.1 adds a new method: SkyCoord.directional_offset_by()
             that provides the same feature as DisplacementMode().
2020-MAR-20: Become compatible with Python3.
             Add xy_comp for DistanceMode()


# Usage:

The astro_coor module provides the command-line interface.
Please do the following to check their usage,
% export PATH=[PATH of the directory of this file]:${PATH}
% python [RADec_convert.py|beam_calc.py|em_wave.py] -h

The astro_coor module also several useful class and functions.
To import this module for using the Coor() class and the other functions , please do
>>> import sys
>>> sys.path.append([PATH of the directory of this file])
>>> from astro_coor import *

"""


# import modules
################

from .astro_coor import *

###################################################################################################
