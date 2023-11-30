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

This module was written because the feature of `astro_coor.DisplacementMode()`
was not implemented by `astropy 2`. However, `astropy 3.1` starts to provide such
a feature with `astropy.coordinates.SkyCoord.directional_offset_by()`,
and then the `astro_coor.Coor` class can almost be replaced by
`astropy.coordinates.SkyCoord`. Despite the great compatibility with the other
 `astropy` feature of using `SkyCoord`, the command-line interface provided by
 `astro_coor` remains very useful.
The CLI commands are the following.

1. `RADec_convert.py`

* Convert given RA and Dec into different forms.
* Calculate the distance and P.A. between two given points.
* Calculate the shifted location from a given point with the displacement.

2. `beam_calc.py` A python script calculating the beam size from given
freq/wavelength and dish diameter.

3. `em_wave.py` A python script calculating the freq from a given wavelength
and vice versa.


## Usage:

The `astro_coor` module provides the command-line interface.
Please do the following to check their usage.
```
$ export PATH=[the path to this 'astro_coor/' directory]:${PATH}
$ python RADec_convert.py -h
$ python beam_calc.py -h
$ python em_wave.py -h
```

The `astro_coor` module also provides several useful class and functions.
To import this module for using the `Coor()` class and the other functions,
please do
```
>>> import sys
>>> sys.path.append([the path to the parent directory of astro_coor/])
>>> from astro_coor import *
```

## History:

2014-JUL-14: Init. The `astro_coor` module has the following functions.
* `CosineLaw_side()`, `CosineLaw_angle()`: Spherical cosine law.
* `DistanceMode()`: Compute the great-circle distance.
* `DisplacementMode()`: Compute a new location by an angular shift. This is a feature that astropy (v2.x) doesn't include.

2017-AUG-15: Write the method to classify types of coord. forms.

2017-AUG-29: Add some functions about units, and calculating beam sizes (for `calc_beam.py`).

2018-JUL-30: Use property methods.

2018-AUG-29: Modify `wave_unit` a little bit (for `em_wave.py`).

*2019-FEB*:   astropy v3.1 adds a new method: `SkyCoord.directional_offset_by()` that provides the same feature as `DisplacementMode()`.

2020-MAR-20: Become compatible with Python3. Add `xy_comp` for `DistanceMode()`.


"""


# import modules
################

from .astro_coor import *

###################################################################################################
