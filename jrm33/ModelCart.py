import numpy as np
from . import Model

def ModelCart(x,y,z,Deg=13):
	'''
	JRM33 Magnetic field model (see Connerney et al 2022 below). The 
	model uses right-handed System III coordinates (I think). 
	
	Inputs
	======
	x : float
		x-coordinate in Rj, R-H System III.
	y : float
		y-coordinate in Rj, R-H System III.
	z : float
		z-coordinate in Rj, R-H System III.
	Deg : int
		Maximum degree of the model, valid 1 - 30 (default = 13).
				
	Returns
	=======
	Bx : float
		x component of magnetic field, nT.
	By : float
		y component of magnetic field, nT.
	Bz : float
		z component of magnetic field, nT.
		
	If using the JRM33 model, please cite the following paper:
	
	Connerney, J. E. P., Timmins, S., Oliversen, R. J., Espley, J. R., 
	Joergensen, J. L., Kotsiaros, S., et al. (2022). A new model of 
	Jupiter's magnetic field at the completion of Juno's Prime Mission. 
	Journal of Geophysical Research: Planets, 127, e2021JE007055. 
	https://doi.org/10.1029/2021JE007055
	
	'''	
	
	#convert to spherical polar coords
	r = np.sqrt(x**2 + y**2 + z**2)
	theta = np.arccos(z/r)
	phi = (np.arctan2(y,x) + (2*np.pi)) % (2*np.pi)
	
	#call the model
	Br,Bt,Bp = Model(r,theta,phi,Deg)

	#convert to Cartesian (hopefully correctly...)
	cost = np.cos(theta)
	sint = np.sin(theta)
	cosp = np.cos(phi)
	sinp = np.sin(phi)
	Bx = Br*sint*cosp + Bt*cost*cosp - Bp*sinp
	By = Br*sint*sinp + Bt*cost*sinp + Bp*cosp
	Bz = Br*cost - Bt*sint
	
	return Bx,By,Bz
	

