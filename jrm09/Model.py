import numpy as np
from . import _SphHarm

def Model(r,theta,phi,MaxDeg=10):
	'''
	JRM09 Magnetic field model (see Connerney et al 2018 below). The 
	model uses right-handed System III coordinates (I think). 
	
	Inputs
	======
	r : float
		Radial distance in Rj.
	theta : float
		Colatitude in radians.
	phi : float
		East longitude in radians.
	MaxDeg : int
		Maximum degree of the model, valid 1 - 10 (default = 10).
		
	Returns
	=======
	Br : float
		Radial field
	Bt : float
		Meridional field
	Bp : float
		Azimuthal field
		
	If using the JRM09 model, please cite the following paper:
	
	Connerney, J. E. P., Kotsiaros, S., Oliversen, R. J., Espley, J. R., 
	Joergensen, J. L., Joergensen, P. S., et al. (2018). A new model of 
	Jupiter's magnetic field from Juno's first nine orbits. Geophysical 
	Research Letters, 45, 2590– 2596. https://doi.org/10.1002/2018GL077312
	
	'''
	
	#make all the inputs arrays
	_r = np.array([r]).flatten()
	_t = np.array([theta]).flatten()
	_p = np.array([phi]).flatten()
	
	return _SphHarm._SphHarm(_r,_t,_p,MaxDeg)
