import numpy as np
from _SphHarm import _SphHarm,_SphHarmArr

def Model(r,theta,phi,Deg=13):
	'''
	JRM33 Magnetic field model (see Connerney et al 2022 below). The 
	model uses right-handed System III coordinates (I think). 
	
	Inputs
	======
	r : float
		Radial distance in Rj.
	theta : float
		Colatitude in radians.
	phi : float
		East longitude in radians.
	Deg : int
		Maximum degree of the model, valid 1 - 30 (default = 13).
		
	Returns
	=======
	Br : float
		Radial field
	Bt : float
		Meridional field
	Bp : float
		Azimuthal field
		
	If using the JRM33 model, please cite the following paper:
	
	Connerney, J. E. P., Timmins, S., Oliversen, R. J., Espley, J. R., 
	Joergensen, J. L., Kotsiaros, S., et al. (2022). A new model of 
	Jupiter's magnetic field at the completion of Juno's Prime Mission. 
	Journal of Geophysical Research: Planets, 127, e2021JE007055. 
	https://doi.org/10.1029/2021JE007055
	
	'''
	
	#get the length of the input
	n = np.size(r)
	
	#whether it is an array
	isarr = hasattr(r,'__iter__')
	
	if n == 1 and isarr:
		return _SphHarm(r[0],theta[0],phi[0],Deg)
	elif n == 1:
		return _SphHarm(r,theta,phi,Deg)
	else:
		#get shape
		sh = np.shape(r)
		_r = r.flatten()
		_t = theta.flatten()
		_p = phi.flatten()
		
		Br,Bt,Bp = _SphHarmArr(_r,_t,_p,Deg)
		return Br,Bt,Bp
	

