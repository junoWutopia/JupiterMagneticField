# JupiterMagneticField
JRM09 model: magnetic field data from juno after completing 9 orbits


Test.py contains 2 functions, map2d(R=0.85, MaxDeg=13) and vecfld(extshell_rad=1.0, layer=3, MaxDeg=10).

map2d generates 3 .png files as the cylindrical map of world(jrmXXtest.png), the azimuthal map of the north hemisphere(jrmXXtestNpolar) and the south hemisphere(jrmXXtestSpolar). The parameter R represents the radius of the shell(map surface).

vecfld shows 3d vector field (module: mayavi.mlab) of the layer with different radii. The parameter layer represents the number of the layers; extshell_rad means the outermost layer radius.
