#!/usb/bin/env python
import math
import os

wedge_patches = ["bottomEmptyFaces","topEmptyFaces"]
axis_patches  = ["axis"] 


ifile = open('constant/polyMesh/points','r')
lines = ifile.readlines()
ifile.close()

for k, line in enumerate(lines):
    if line[0] == "(":
        start = int(k) + 1
        break
    
number = int(lines[k-1])
end = start + number 


# xy-plane
theta = 2.5/360*2*math.pi
for k, line  in enumerate(lines[start: end]):
    pieces = line[1:-2].split()
    x = float(pieces[0])
    y = float(pieces[1])    
    z = float(pieces[2])
    dist = y*math.asin(theta)
    if z < 0:
        new_z = -dist
    elif z >= 0 :
        new_z = dist
    else:
        raise ValueError()
    lines[k+start] = "(" + pieces[0] + " " + pieces[1] + " " + str(new_z) + ")\n"
    
ofile = open('constant/polyMesh/points','w+')
ofile.writelines(lines)
ofile.close()

# Collapes
os.system('collapseEdges -overwrite')


# Manipulate boundary-file
ifile = open('constant/polyMesh/boundary','r')
lines = ifile.readlines()
ifile.close()

for wedge_patch in wedge_patches:
    for k, line in enumerate(lines):
        if line.strip() == wedge_patch:
            lines[k+2] = "\t\ttype\t\t\twedge;\n"
            lines[k+3] = "\t\tinGroups\t\t1(wedge);\n"    
            break            
        
for axis_patch in axis_patches:
    for k, line in enumerate(lines):
        if line.strip() == axis_patch:
            lines[k+2] = "\t\ttype\t\t\tempty;\n"
            lines[k+3] = "\t\tinGroups\t\t1(empty);\n"    
            break             
        
ofile = open('constant/polyMesh/boundary','w+')
ofile.writelines(lines)
ofile.close()

