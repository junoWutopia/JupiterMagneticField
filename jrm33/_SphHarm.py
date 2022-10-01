import numpy as np
from . import _CoeffGrids
from ._Legendre import _Legendre,_LegendreArr
from . import _Schmidt
from . import _ScaledCoeffs


def _SphHarm(r,theta,phi,Deg=13):
	'''
	This function calculates the JRM33 model field using spherical
	harmonics. 
	
	Inputs
	======
	r : float
		Radial coordinate in Rj (69,911 km).
	theta : float
		Colatitude in radians (RH SIII).
	phi : float
		Azimuth in radians (RH SIII).
	Deg : int
		Maximum degree of the model to use - can be 1 to 30 (default=13)
		for JRM33, where larger degrees take longer to process, but are
		better.
	
	Returns
	=======
	Br : float
		Magnetic field strength in radial direction (nT).
	Bt : float
		Magnetic field strength in latitudinal direction (nT).
	Bp : float
		Magnetic field strength in azimuthal direction (nT).
	
	'''
	

	#calculate cosmphi and sinmphi
	cosmphi = np.zeros((Deg+1),dtype='float64') + 1.0
	sinmphi = np.zeros((Deg+1),dtype='float64')
	for m in range(1,Deg+1):
		mphi = m*phi
		cosmphi[m] = np.cos(mphi)
		sinmphi[m] = np.sin(mphi)
	
	#get the grids of scaled g and h parameters
	g,h = _ScaledCoeffs()
	
	#output arrays
	Br = 0.0
	Bt = 0.0
	Bp = 0.0
	
	#start calculating p costheta and its derivative
	Pnm,dPnm = _Legendre(theta,Deg)

	#temporary arrays for the inner sums
	sumr = 0.0
	sumt = 0.0
	sump = 0.0
	
	#now sum everything up
	r1 = 1/r
	sintheta1 = 1.0/np.sin(theta)
	if np.isfinite(sintheta1) == False:
		sintheta1 = 0.0
	C = r1*r1
	for n in range(1,Deg+1):
		#define the constant (a/r)**(n+2)
		C = C*r1

		#do the innder sum for this value of n first
		#I suspect this could be vectorised, but would probably be 
		#pretty unreadable if I did
		sumr = 0.0
		sumt = 0.0
		sump = 0.0
		for m in range(0,n+1):
			gchs = g[n,m]*cosmphi[m] + h[n,m]*sinmphi[m]
			sumr += Pnm[n,m]*gchs
			sumt += dPnm[n,m]*gchs
			sump += m*Pnm[n,m]*(h[n,m]*cosmphi[m] - g[n,m]*sinmphi[m])

		#add it to the output arrays
		Br += C*(n+1)*sumr
		Bt += -C*sumt
		Bp += -C*sump
	Bp = sintheta1*Bp

	return Br,Bt,Bp


def _SphHarmArr(r,theta,phi,Deg=10):
	'''
	This function calculates the JRM09 model field using spherical
	harmonics. 
	
	Inputs
	======
	r : float
		Radial coordinate in Rj (69,911 km).
	theta : float
		Colatitude in radians (RH SIII).
	phi : float
		Azimuth in radians (RH SIII).
	Deg : int
		Maximum degree of the model to use - can be 1 to 10 (default=10)
		for JRM09, where larger degrees take longer to process, but are
		better.
	
	Returns
	=======
	Br : float
		Magnetic field strength in radial direction (nT).
	Bt : float
		Magnetic field strength in latitudinal direction (nT).
	Bp : float
		Magnetic field strength in azimuthal direction (nT).
	
	'''
	

	#calculate cosmphi and sinmphi
	nr = np.size(r)
	cosmphi = np.zeros((Deg+1,nr),dtype='float64') + 1.0
	sinmphi = np.zeros((Deg+1,nr),dtype='float64')
	for m in range(1,Deg+1):
		mphi = m*phi
		cosmphi[m] = np.cos(mphi)
		sinmphi[m] = np.sin(mphi)

	#get the grids of g and h parameters
	g,h = _ScaledCoeffs._ScaledCoeffs()
	
	#output arrays
	Br = np.zeros(r.size,dtype='float64')
	Bt = np.zeros(r.size,dtype='float64')
	Bp = np.zeros(r.size,dtype='float64')
	

	#start calculating p costheta and its derivative
	Pnm,dPnm = _LegendreArr(theta,Deg)


	#temporary arrays for the inner sums
	sumr = np.zeros((r.size,),dtype='float64')
	sumt = np.zeros((r.size,),dtype='float64')
	sump = np.zeros((r.size,),dtype='float64')
	
	#now sum everything up
	r1 = 1/r
	sintheta1 = 1.0/np.sin(theta)
	badst = np.where(np.isfinite(sintheta1) == False)[0]
	if badst.size > 0:
		sintheta1[badst] = 0.0
	C = (r1)**2
	for n in range(1,Deg+1):
		#define the constant (a/r)**(n+2)
		C = C*r1

		#do the innder sum for this value of n first
		#I suspect this could be vectorised, but would probably be 
		#pretty unreadable if I did
		sumr[:] = 0
		sumt[:] = 0
		sump[:] = 0
		for m in range(0,n+1):
			gchs = g[n,m]*cosmphi[m] + h[n,m]*sinmphi[m]
			sumr += Pnm[n,m]*gchs
			sumt += dPnm[n,m]*gchs
			sump += m*Pnm[n,m]*(h[n,m]*cosmphi[m] - g[n,m]*sinmphi[m])

		#add it to the output arrays
		Br += C*(n+1)*sumr
		Bt += -C*sumt
		Bp += -C*sump
	Bp = sintheta1*Bp

	return Br,Bt,Bp
