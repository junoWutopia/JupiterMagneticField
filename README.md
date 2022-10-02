# JupiterMagneticField
JRM09 model: magnetic field data from juno after completing 9 orbits

JRM33 model: magnetic field data from juno after completing 33 orbits

Test.py contains 2 functions, map2d(R=0.85, MaxDeg=13) and vecfld(extshell_rad=1.0, layer=3, MaxDeg=10).

map2d generates 3 .png files as the cylindrical map of world(jrmXXtest.png), the azimuthal map of the north hemisphere(jrmXXtestNpolar) and the south hemisphere(jrmXXtestSpolar). The parameter R represents the radius of the shell(map surface).

vecfld shows 3d vector field (module: mayavi.mlab) of the layer with different radii. The parameter layer represents the number of the layers; extshell_rad means the outermost layer radius.

# References
Connerney, J. E. P., Kotsiaros, S., Oliversen, R. J., Espley, J. R., Joergensen, J. L., Joergensen, P. S., et al. (2018). A new model of Jupiter's magnetic field from Juno's first nine orbits. Geophysical Research Letters, 45, 2590– 2596. https://doi.org/10.1002/2018GL077312

Connerney, J. E. P., Timmins, S., Oliversen, R. J., Espley, J. R., Joergensen, J. L., Kotsiaros, S., et al. (2022). A new model of Jupiter's magnetic field at the completion of Juno's Prime Mission. Journal of Geophysical Research: Planets, 127, e2021JE007055. A New Model of Jupiter's Magnetic Field at the Completion of Juno's Prime Mission - Connerney - 2022 - Journal of Geophysical Research: Planets - Wiley Online Library

C.K. Goertz, D.E. Jones, B.A. Randall, E.J. Smith, M.F. Thomsen, J. Geophys. Res. 81, (1976)
