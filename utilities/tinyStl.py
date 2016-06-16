#!/usr/bin/python

from __future__ import division

from os.path import basename, splitext

##############################################################################

def write_xy_line(xs,ys, zs, path, 
                  label=None, closed=False, float_separator="."):

    # If no label given, set filename with no extension as label
    if not label:
        label = splitext(basename(path))[0]
    write(xy_line(xs,ys, zs, label, closed, float_separator), path)                          

def xy_line(xs,ys, zs, label="nameless", closed=False, float_separator=".", extra=0):
    
    # Allow not definining other z
    if type(zs) == float or type(zs) == int:
        z1 = zs
        z0 = 0    
    else:
        z0 = zs[0]
        z1 = zs[1]
    
    lines_str = "solid {label}".format(label=label)
    for k in range(len(xs)-1):
        lines_str += _xy_box(xs[k],xs[k+1], ys[k], ys[k+1], z0, z1)
    if closed:
        lines_str += _xy_box(xs[-1],xs[0], ys[-1], ys[0], z0, z1)
    lines_str += "\n" + "endsolid {label}".format(label=label) +"\n"
    lines_str = lines_str.replace(".",float_separator)
    return lines_str

def write(stl, path):
    ofile = open(path, 'w+')
    ofile.write(stl)
    ofile.close()

def _xy_box(x0,x1, y0, y1, z0, z1):
#    print x0,x1, y0, y1, z0, z1
#    print type(x0),type(x1), type(y0), type(y1), type(z0), type(z1)
    return """
  facet normal 0 0 0
    outer loop
      vertex {x0:e} {y0:e} {z0:e}
      vertex {x1:e} {y1:e} {z0:e}
      vertex {x1:e} {y1:e} {z1:e}
    endloop
  endfacet
  facet normal 0 0 0
    outer loop
      vertex {x0:e} {y0:e} {z0:e}
      vertex {x1:e} {y1:e} {z1:e}
      vertex {x0:e} {y0:e} {z1:e}
    endloop
  endfacet""".format(x0=x0, x1=x1, 
           y0=y0, y1=y1, z0=z0, z1=z1)

##############################################################################

if __name__ == "__main__":
    print "START"
    print """
    Unit test. See source for detail.
    
    Writes two ascii stl files with a box in xy-plane with
    small extrusion to -z-direction.
    File ./box.stl contains a box labeled \"box\" and
    file ./named_box.stl contains a similar box with sides labeled as 
    \"bottom\", \"right\", \"top\", and \"left\".    
    
    A third stl file \"./ellipse.stl\" contains an ellipse to 
    illustrate rounds shapes. An arbitrary 2d form in xy-plane
    with an arbitrary resolutions can be drawn. Just provide 
    the lists of coordinates.
    
    Written by Antti Mikkonen, a.mikkonen@iki.fi, 2015
""" 
    # Box   
    xs     = [0, 2, 2, 0] # List of x-coordinates.
    ys     = [0, 0, 1, 1] # List of y-coordinates.   
    dz     = -1e-1        # Z-direction extrusion.  
    
    write_xy_line(xs,ys, dz, path="box.stl", label="box", closed=True)
    
    # Labeled box
    bottom = xy_line([0,2], [0,0], dz, "bottom")
    right  = xy_line([2,2], [0,1], dz, "right")
    top    = xy_line([2,0], [1,1], dz, "top")    
    left   = xy_line([0,0], [1,0], dz, "left")
    
    write(bottom+right+top+left, "labeled_box.stl")
        
    # Ellipse
    try:
        import numpy as np
    except:
        print "*"*50
        print "No numpy found, skipping the ellipse example."
        print "Stl generation is unaffected. Numpy is only needeed"
        print "for the ellipase example."
        print "*"*50        
    else:
        theta = np.linspace(0,2*np.pi, 100)
        a = 2; b = 1
        xs = a * np.cos(theta)
        ys = b * np.sin(theta)
        write_xy_line(xs,ys, dz, "ellipse.stl", "ellipse", True)
        
    print "END"
