#! /usr/bin/env python
#-*- coding:utf-8 -*-
# Written on 6/14/2014 by Sheng-Jun Lin
# 8/10/2014: Add Angular Distance Mode
# 8/15/2017: Move classes to be the module astro_coor

from astro_coor import *
from sys import argv
from getopt import getopt
from subprocess import call, PIPE, Popen

def resolve_RA_line(string, line):

    try:
        RA = ra(ra_str=string)
        return RA
    except ValueError as err_msg:
        raise ValueError('RADec_convert: [line: {0}] {1}'.format(line, str(err_msg)))


def resolve_Dec_line(string, line):

    try:
        Dec = dec(dec_str=string)
        return Dec
    except ValueError as err_msg:
        raise ValueError('RADec_convert: [line: {0}] {1}'.format(line, str(err_msg)))

def main():

    help_doc = 'RADec_convert.py: \n\
A Python script to calculate all the RA/Dec forms, the angular distance,\n\
                   or the angular displaccement of the coordinates.\n\n\
Usage: \n\
1) Input from a file:           python RADec_convert.py [OPTION]...  [INPUT FILE]\n\
2) Input from the command line: python RADec_convert.py -c [OPTION].. [COORDINATES]...\n\
  -d             Change to [Angular Distance Mode].\n\
  -x             Change to [Angular Displacement Mode] (given disp. & PA).\n\
  -y             Change to [Angular Displacement Mode] (given dx & dy).\n\
  -c             Type in comand line directly and don\'t assign a input file.\n\
  -n             Don\'t generate the output file.\n\
  -o output_file Specify the output filename. Otherwise it will be named as RADec.output\n\
         If the file already exists, the result will be appended to the file.\n\
  -s             Don\'t print the result on screen.\n\
  -h             Print this help information.\n\n\
Input file format:\n\
 [Converter Mode] (the default process)\n\
     For RA, there are 3 forms to choose.\n\
       1) hour(int) : min(int)   : sec(float)\n\
       2) deg(int)  : arcmin(int): arcsec(float)\n\
       3) degree (float)\n\
     For Dec, there are 2 forms to choose.\n\
       1) +/-deg(int) : arcmin(int): arcsec(float), where \'+\' could be ignored.\n\
       2) degree (float)\n\
     Please write coordinates into each line,\n\
     where the 1st component is RA, and the 2nd one is Dec.\n\
     For the 2nd form of RA, you need to prefix a keyword \'d\' without space.\n\
     For integer vars hour/min/deg/arcmin, you could input float vars,\n\
     then the decimal part will be combined to the next subunit.\n\
     And the line which starts from \'#\' will be ignored.\n\
     Example:\n\
              10:23:47.68  0:38:41.3\n\
              155.948674   00.644805\n\
              d155:56:55.2 -0:38:41.3\n\
 [Angular Distance Mode] (option: -d)\n\
     The Format: RA(1st)  Dec(1st)  RA(2nd)  Dec(2nd)\n\
     For each line, please write the coordinates of 2 points, in whatever forms below.\n\
     Example:\n\
              10:23:47.68  0:38:41.3  10:01:55  0:30:40\n\
 [Angular Displacement Mode]\n\
   (1) Given displacement & P.A. (option: -x)\n\
     The Format: RA  Dec  displacement(the default is in arcsec)  P.A.(the default is in deg)\n\
     For each line, please write the coordinates of a point, in whatever forms below.\n\
     Then write the displacement, and the position angle.\n\
     One can use \'arcsec\', \'as\', \'arcmin\', \'am\', \'deg\', \'d\' to change the unit.\n\
     Example:\n\
              10:23:47.68  0:38:41.3  328.251247arcmin -91.371639\n\
   (2) Given dx & dy (option: -y)\n\
     The Format: RA  Dec  dx(the default is in arcsec)  dy(the default is in arcsec)\n\
     One can use \'arcsec\', \'as\', \'arcmin\', \'am\', \'deg\', \'d\' to change the unit.\n\
     Example:\n\
              10:23:47.68  0:38:41.3  10arcsec 10arcsec\n\n\
Output file format:\n\
    The 1st line is the time of executing this command.\n\
    [Converter Mode] List all forms of input coordinates.\n\
    [Angular Distance Mode]\n\
        List 2 coordinates in hms and dms form,\n\
        the angular distance in arcsec, arcmin, and degree, and the position angle of them.\n\
    [Angular Displacement Mode]\n\
        List all forms of input and output coordinates,\n\
        and the angular distance in arcsec, arcmin, and degree.'

    opts, args = getopt(argv[1:], "ho:nsdxyc")

    exec_cmdline = False
    on_screen = True
    gen_file = True
    distanceMode = False
    displacementMode = False
    outputfile = 'RADec.output'

    for opt, value in opts:
        if opt == '-h':
            print(help_doc)
            exit()
        if opt == '-o':
            outputfile = value
        if opt == '-s':
            on_screen = False
        if opt == '-n':
            gen_file = False
        if opt == '-d':
            distanceMode = True
        if opt == '-x':
            displacementMode = True
            xy_comp = False
        if opt == '-y':
            displacementMode = True
            xy_comp = True
        if opt == '-c':
            exec_cmdline = True

    if gen_file:
        cmd = Popen(['ls'], stdout=PIPE)
        shell_result = cmd.communicate()[0]
        if outputfile not in shell_result:
            call(['touch', outputfile])
            txt = open(outputfile, 'a')
            call(['date'], stdout=txt, shell=True)
        else:
            txt = open(outputfile, 'a')
            call(['date'], stdout=txt, shell=True)
    if exec_cmdline:
        items = args
        database = [1]
    else:
        inputfile = open(args[0], 'r')
        database = inputfile.readlines()
        inputfile.close()

    n = 1
    result = ''
    for i in xrange(len(database)):
        # Reading the input coordinates
        if not exec_cmdline:
            items = database[i].split()
        if len(items) == 0 or '#' in items[0]:
            continue
        data_ra1 = items[0]
        data_dec1 = items[1]

        # Resolving the input coordinates
        RA1 = resolve_RA_line(data_ra1, i + 1)
        Dec1 = resolve_Dec_line(data_dec1, i + 1)

        # Calculating all the representation of input coordinates
        RA1.converter()
        Dec1.converter()

        if distanceMode:
            # For calculating angular distance
            data_ra2 = items[2]
            data_dec2 = items[3]
            RA2 = resolve_RA_line(data_ra2, i + 1)
            Dec2 = resolve_Dec_line(data_dec2, i + 1)
            RA2.converter()
            Dec2.converter()
            d_deg, pa21, pa12, dx_as, dy_as = DistanceMode(RA1, Dec1, RA2, Dec2, xy_comp=True)

            # Write into the outputfile [anglar distance]
            result += ('#{0:-<4d}--{1:-<12s}--{2:-<14s} '
                       'dist = {3:>11.3f}" = {4:>9.3f}\' = {5:>7.3f} deg\n'.format(
                       n, 'R.A.', 'Dec.', d_deg * 3600., d_deg * 60., d_deg))
            n += 1
            result += ('1)hms {0:>2d}:{1:0>2d}:{2:0>7.4f}  '
                       '{3}{4:0>2d}:{5:0>2d}:{6:0>6.3f}  '
                       'dx  = {7:>11.3f}" = {8:>9.3f}\' = {9:>7.3f} deg\n'.format(
                       RA1.h, RA1.m, RA1.s,
                       Dec1.sign, Dec1.d, Dec1.arcm, Dec1.arcs,
                       dx_as, dx_as / 60., dx_as / 3600.))
            result += ('2)hms {0:>2d}:{1:0>2d}:{2:0>7.4f}  '
                       '{3}{4:0>2d}:{5:0>2d}:{6:0>6.3f}  '
                       'dy  = {7:>11.3f}" = {8:>9.3f}\' = {9:>7.3f} deg\n'.format(
                       RA2.h, RA2.m, RA2.s,
                       Dec2.sign, Dec2.d, Dec2.arcm, Dec2.arcs,
                       dy_as, dy_as / 60., dy_as / 3600.))
            result += (' '*36 + 'PA   = {0:+8.3f} (2wrt1) / {1:+8.3f} (1wrt2) deg\n'.format(pa21, pa12))
        elif displacementMode:
            # For calculating angular displacement
            ## given dist, pa
            if not xy_comp:
                dist = items[2]
                pa = items[3]
                if dist[-1] == 'd':
                    dist = float(dist[:-1]) * 3600
                elif dist[-3:] == 'deg':
                    dist = float(dist[:-3]) * 3600
                elif dist[-2:] == 'am':
                    dist = float(dist[:-2]) * 60.
                elif dist[-6:] == 'arcmin':
                    dist = float(dist[:-6]) * 60.
                elif dist[-2:] == 'as':
                    dist = float(dist[:-2])
                elif dist[-6:] == 'arcsec':
                    dist = float(dist[:-6])
                else:
                    dist = float(dist)
                if pa[-2:] == 'as':
                    pa = float(pa[:-2]) / 3600.
                elif pa[-6:] == 'arcsec':
                    pa = float(pa[:-6]) / 3600.
                elif pa[-2:] == 'am':
                    pa = float(pa[:-2]) / 60.
                elif pa[-6:] == 'arcmin':
                    pa = float(pa[:-6]) / 60.
                elif pa[-1] == 'd':
                    pa = float(pa[:-1])
                elif pa[-3:] == 'deg':
                    pa = float(pa[:-3])
                else:
                    pa = float(pa)
                RA2, Dec2 = DisplacementMode(RA1, Dec1, dist, pa)
            else:
            ## given dx, dy
                dx = items[2]
                dy = items[3]
                if dx[-2:] == 'as':
                    dx = float(dx[:-2])
                elif dx[-6:] == 'arcsec':
                    dx = float(dx[:-6])
                elif dx[-2:] == 'am':
                    dx = float(dx[:-2]) * 60.
                elif dx[-6:] == 'arcmin':
                    dx = float(dx[:-6]) * 60.
                elif dx[-1] == 'd':
                    dx = float(dx[:-1]) * 3600.
                elif dx[-3:] == 'deg':
                    dx = float(dx[:-3]) * 3600.
                else:
                    dx = float(dx)
                if dy[-2:] == 'as':
                    dy = float(dy[:-2])
                elif dy[-6:] == 'arcsec':
                    dy = float(dy[:-6])
                elif dy[-2:] == 'am':
                    dy = float(dy[:-2]) * 60.
                elif dy[-6:] == 'arcmin':
                    dy = float(dy[:-6]) * 60.
                elif dy[-1] == 'd':
                    dy = float(dy[:-1]) * 3600.
                elif dy[-3:] == 'deg':
                    dy = float(dy[:-3]) * 3600.
                else:
                    dy = float(dy)
                RA2, Dec2 = DisplacementMode(RA1, Dec1, dx, dy, xy_comp=True)

            dist_deg, pa21, pa12, dx_as, dy_as = DistanceMode(RA1, Dec1, RA2, Dec2, xy_comp=True)

            # Write into the outputfile [anglar displacement]
            result += ('#{0:-<4d}--{1:-<12s}--{2:-<14s} '
                       'dist = {3:>11.3f}" = {4:>9.3f}\' = {5:>7.3f} deg\n'.format(
                       n, 'R.A.', 'Dec.', dist_deg * 3600., dist_deg * 60., dist_deg))
            n += 1
            result += ('1)hms  {0:>2d}:{1:0>2d}:{2:0>7.4f}  '
                       '{3}{4:0>2d}:{5:0>2d}:{6:0>6.3f} '
                       'dx  = {7:>11.3f}" = {8:>9.3f}\' = {9:>7.3f} deg\n'.format(
                       RA1.h, RA1.m, RA1.s,
                       Dec1.sign, Dec1.d, Dec1.arcm, Dec1.arcs,
                       dx_as, dx_as / 60., dx_as / 3600.))
            result += ('  dms {0:>3d}:{1:0>2d}:{2:0>7.4f}                '
                       'dy  = {3:>11.3f}" = {4:>9.3f}\' = {5:>7.3f} deg\n'.format(
                       RA1.d, RA1.arcm, RA1.arcs,
                       dy_as, dy_as / 60., dy_as / 3600.))
            result += ('  deg {0:>14.10f}  {1:+.10f} '
                       'PA   = {2:+8.3f} (2wrt1) / {3:+8.3f} (1wrt2) deg\n'.format(
                      RA1.all_d, Dec1.all_d, pa21, pa12))
            result += ('2)hms  {0:>2d}:{1:0>2d}:{2:0>7.4f}  '
                       '{3}{4:0>2d}:{5:0>2d}:{6:0>6.3f}\n'.format(
                       RA2.h, RA2.m, RA2.s,
                       Dec2.sign, Dec2.d, Dec2.arcm, Dec2.arcs))
            result += '  dms {0:>3d}:{1:0>2d}:{2:0>7.4f}\n'.format(
                       RA2.d, RA2.arcm, RA2.arcs)
            result += '  deg {0:>14.10f}  {1:+.10f}\n'.format(
                       RA2.all_d, Dec2.all_d)

        else:
            # Write into the outputfile [coordinates converter]
            result += '#{0:-<3d}--{1:-<12s}----{2:-<1s}\n'.format(
                      n, 'R.A.', 'Dec.')
            n += 1
            result += ('hms  {0:>2d}:{1:0>2d}:{2:0>7.4f}  '
                       '{3}{4:0>2d}:{5:0>2d}:{6:0>7.4f}\n'.format(
                       RA1.h, RA1.m, RA1.s,
                       Dec1.sign, Dec1.d, Dec1.arcm, Dec1.arcs))
            result += 'dms {0:>3d}:{1:0>2d}:{2:0>7.4f}\n'.format(
                       RA1.d, RA1.arcm, RA1.arcs)
            result += 'deg {0:>14.10f}  {1:+.10f}\n'.format(
                       RA1.all_d, Dec1.all_d)
    result += '=' * 35 + '\n'
    if gen_file:
        txt.write(result + '\n')
        txt.close()
    if on_screen:
        print(result)


if __name__ == '__main__':
    main()
