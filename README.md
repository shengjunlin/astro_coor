# astro_coor
A Python module to calculate all the RA/Dec forms, the angular distance,
the angular displaccement of the coordinates. There are also additional
functions handling the beam size calculation.

This module was written because the feature of `astro_coor.DisplacementMode()`
was not implemented by astropy 2.x. While astropy 3.1 starts to provide such
a feature with `astropy.coordinates.SkyCoord.directional_offset_by()`,
and then the `astro_coor.Coor` class can be almost replaced with
` astropy.coordinates.SkyCoord`.
However, the command-line interface is still useful.
The commands are the following.

1. `RADec_convert.py`

Features:
* Convert given RA and Dec into different forms.
* Calculate the distance and P.A. between two given points.
* Calculate the shifted location from a given point with the displacement.

2. `beam_calc.py` A python script calculating the beam size from given freq/wavelength and diameter.

3. `em_wave.py` A python script calculating the freq from a given wavelength and vice versa.


## Usage:

The `astro_coor` module provides the command-line interface.
Please do the following to check their usage.
```
$ export PATH=[the path to this 'astro_coor/' directory]:${PATH}
$
$ python RADec_convert.py -h
RADec_convert.py:
A Python script to calculate all the RA/Dec forms, the angular distance,
                   or the angular displaccement of the coordinates.

Usage:
1) Input from a file:           python RADec_convert.py [OPTION]...  [INPUT FILE]
2) Input from the command line: python RADec_convert.py -c [OPTION].. [COORDINATES]...
  -d             Change to [Angular Distance Mode].
  -x             Change to [Angular Displacement Mode] (given disp. & PA).
  -y             Change to [Angular Displacement Mode] (given dx & dy).
  -c             Type in comand line directly and don't assign a input file.
  -n             Don't generate the output file.
  -o output_file Specify the output filename. Otherwise it will be named as RADec.output
         If the file already exists, the result will be appended to the file.
  -s             Don't print the result on screen.
  -h             Print this help information.

Input file format:
 [Converter Mode] (the default process)
     For RA, there are 3 forms to choose.
       1) hour(int) : min(int)   : sec(float)
       2) deg(int)  : arcmin(int): arcsec(float)
       3) degree (float)
     For Dec, there are 2 forms to choose.
       1) +/-deg(int) : arcmin(int): arcsec(float), where '+' could be ignored.
       2) degree (float)
     Please write coordinates into each line,
     where the 1st component is RA, and the 2nd one is Dec.
     For the 2nd form of RA, you need to prefix a keyword 'd' without space.
     For integer vars hour/min/deg/arcmin, you could input float vars,
     then the decimal part will be combined to the next subunit.
     And the line which starts from '#' will be ignored.
     Example:
              10:23:47.68  0:38:41.3
              155.948674   00.644805
              d155:56:55.2 -0:38:41.3
 [Angular Distance Mode] (option: -d)
     The Format: RA(1st)  Dec(1st)  RA(2nd)  Dec(2nd)
     For each line, please write the coordinates of 2 points, in whatever forms below.
     Example:
              10:23:47.68  0:38:41.3  10:01:55  0:30:40
 [Angular Displacement Mode]
   (1) Given displacement & P.A. (option: -x)
     The Format: RA  Dec  displacement(the default is in arcsec)  P.A.(the default is in deg)
     For each line, please write the coordinates of a point, in whatever forms below.
     Then write the displacement, and the position angle.
     One can use 'arcsec', 'as', 'arcmin', 'am', 'deg', 'd' to change the unit.
     Example:
              10:23:47.68  0:38:41.3  328.251247arcmin -91.371639
   (2) Given dx & dy (option: -y)
     The Format: RA  Dec  dx(the default is in arcsec)  dy(the default is in arcsec)
     One can use 'arcsec', 'as', 'arcmin', 'am', 'deg', 'd' to change the unit.
     Example:
              10:23:47.68  0:38:41.3  10arcsec 10arcsec

Output file format:
    The 1st line is the time of executing this command.
    [Converter Mode] List all forms of input coordinates.
    [Angular Distance Mode]
        List 2 coordinates in hms and dms form,
        the angular distance in arcsec, arcmin, and degree, and the position angle of them.
    [Angular Displacement Mode]
        List all forms of input and output coordinates,
        and the angular distance in arcsec, arcmin, and degree.
$
$ python beam_calc.py -h
beam_calc.py:
Calculating beam sizes for Circular Fraunhofer (1st order approx.), and IRAM 30m.
Usage: beam_calc.py freq/wavelength diameter
E.g. ./beam_calc.py 350GHz 12m
freq/wavelength/diameter need to be suffixed with their units.
$
$ python em_wave.py -h
em_wave.py:
Calculating frequency and wavelength of the input EM wave.
Usage: em_wave.py freq/wavelength
E.g. ./em_wave.py 350GHz
freq/wavelength need to be suffixed with their units.
```

The `astro_coor` module also several useful class and functions.
To import this module for using the `Coor()` class and the other functions, please do
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


