"""
Genaration of stretching markers

N.B.: make sure you print without scaling!

Inkscape is recommended for printing
https://inkscape.org/en/

Author: Ruslan Guseinov
"""

## Options #################################################

# List of stretching device diameters to output, mm
diameters = [200, 400, 500]

# Uniform stretching factor of membrane
stretch_factor = 3

# Output paper format (A4/A3/none)
paper = 'A4'

############################################################


import numpy as np

# Generate star

pt = np.empty([16, 2])
for i in range(16):
    alpha = np.pi / 16 + 2 * np.pi * i / 16
    pt[i] = np.array([np.cos(alpha), np.sin(alpha)])

# Output to SVG file

if paper == 'A4':
    frameX = 210
    frameY = 297
elif paper == 'A3':
    frameX = 297
    frameY = 420
else:
    frameX = frameY = max(diameters) / stretch_factor + 10
cnt = [frameX / 2, frameY / 2]

ofile = open('stretch_markers.svg', 'w')

ofile.write('<?xml version="1.0" encoding="utf-8"?>')
ofile.write('<svg ')
ofile.write('xmlns:svg="http://www.w3.org/2000/svg" ')
ofile.write('xmlns="http://www.w3.org/2000/svg" ')
ofile.write('version="1.1" ')
ofile.write('width="' + str(frameX) + 'mm" ')
ofile.write('height="' + str(frameY) + 'mm" ')
ofile.write('viewBox="0 0 ' + str(frameX) + ' ' + str(frameY) + '">\n')
ofile.write('<style type="text/css">\n')
ofile.write('	.st0{fill:#FFFFFF;stroke:#000000;stroke-miterlimit:10;}\n')
ofile.write('</style>\n')

ofile.write('<g transform="translate(' + ",".join(str(a) for a in cnt) + ')">\n')

for diam in diameters[::-1]:
    radius = diam * 0.5 / stretch_factor
    ptscaled = pt * radius
    pts = " ".join([",".join([str(a) for a in b]) for b in ptscaled])
    ofile.write('<polygon class="st0" stroke-width="0.3" points="' + pts + '"/>"\n')
    ofile.write('<text x="0" y="-' + str(radius) + '" font-size="2" text-anchor="middle" fill="black">' + str(diam) + 'mm</text>\n')
    for p in ptscaled:
        ofile.write('<circle cx="' + str(p[0]) + '" cy="' + str(p[1]) + '" r="0.5" />\n')

ofile.write('</g>\n')

ofile.write("</svg>\n")
ofile.close()
