import subprocess
import os
import itertools

openscad = 'C:/Program Files/OpenSCAD' # path to openscad

# Parameter grid
spring_th = [0.3, 0.4, 0.5, 0.6]
d = [5, 6]

# Store params and append include line
with open('params.scad', 'r') as fparams:
    lines = fparams.readlines()
with open('params.scad', 'a') as fparams:
    fparams.write("include <inputs.scad>;\n")

# Write inputs and run openscad

if not os.path.exists('output'):
    os.makedirs('output')

i = 1
total = len(spring_th) * len(d)
for e in itertools.product(spring_th, d):
    with open('inputs.scad', 'w') as finputs:
        finputs.write('spring_th = ' + str(e[0]) + ';\n')
        finputs.write('d = ' + str(e[1]) + ';\n')
    subprocess.call(openscad + '/openscad.exe bracket.scad -o output/bracket_' + str(e[0]) + '_' + str(e[1]) + '.stl')
    print(str(i) + '\t/\t' + str(total) + ' done')
    i = i + 1

# Restore <params.scad> without the include line and remove <inputs.scad>
with open('params.scad', 'w') as fparams:
    fparams.writelines(lines)
os.remove('inputs.scad')
