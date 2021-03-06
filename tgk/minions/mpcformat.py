# Licensed under a MIT style license - see LICENSE
"""mpcsubmit - Format the observation into MPC format."""

# tofix: (1) only support comets with periodic numbers

# INPUT

obj = '0041P'   # object, in MPC format  

obs = [['2017-03-14T10:14:31', '10h17m24.592s', '+44d01m16.50s', 14.26, 'r'], \
['2017-03-14T10:21:25', '10h17m24.658s', '+44d01m18.09s', 14.24, 'r'], \
['2017-03-14T10:25:31', '10h17m24.751s', '+44d01m19.83s', 14.14, 'r'], \
['2017-03-14T10:29:54', '10h17m24.818s', '+44d01m21.56s', 14.10, 'r']]
                # an array of observations in the format of [ [time, RA, Dec, mag, band], [time, RA, Dec, mag, band], ... ]
                # e.g. time -- str('2000-01-01T00:00:00'); RA -- str('00h00m00.000s'); Dec -- str('+00d00m00.000s'); mag -- float(20.00), band -- str('r') -- must be single letter

site = 'Q63'    # need to be read from the header?
                # Sutherland - K94
                # Siding Spring - Q63
                # Haleakala - T04
                # McDonald - V37
                # Cerro Tololo - W85
                # Tenerife - Z21
                
obsver = "\
OBS D. Bodewits, M. S. P. Kelley, M. M. Knight, T. L. Farnham, \n\
OBS S. Protopapa, J.-Y. Li, T. Lister, Q.-Z. Ye, H. Hsieh, Z.-Y. Lin, \n\
OBS E. Jehin, C. Snodgrass, A. Fitzsimmons"
                # observers, in MPC format

cat = 'UCAC-4'
                # catalog used for astrometric/photometric reduction
                
# MAIN

# header

report = "COD %s\n" % site
report += "%s\n" % obsver
report += "TEL 1.0-m reflector + CCD\n"
report += "ACK LCO41P\n"
report += "NET %s\n" % cat

# main body

for obsi in obs:
    tmp = float(float(obsi[0][17:19])/86400. + float(obsi[0][14:16])/1440. + float(obsi[0][11:13])/24.)
    report += "%s         C%s %s %s.%s %s %s %s %s %s %s          %s %s      %s\n" % \
    (obj, str(obsi[0][0:4]), str(obsi[0][5:7]), str(obsi[0][8:10]), str(round(tmp, 5))[2:], \
    str(obsi[1][0:2]), str(obsi[1][3:5]), str(obsi[1][6:11]), \
    str(obsi[2][0:3]), str(obsi[2][4:6]), str(obsi[2][7:11]), \
    str(round(obsi[3], 1)), str(obsi[4]), site)
report += '---'     # make sure the last line won't get truncated
                
# save the file

o = open("MPCReport.txt", "w")
o.write(report)
o.close()
