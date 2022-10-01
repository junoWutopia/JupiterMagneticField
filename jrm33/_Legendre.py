import numpy as np


def _Legendre(theta,Deg=13):
	'''
	This lovely little routine is used to calculate the P(cos(theta))
	Legendre polynomials for the model. It starts off with P_00, P_10
	and P_11 - then recursively calculates the higher degree 
	coefficients.
	
	Inputs
	======
	theta : float
		Colatitude coordinate in radians.
	Deg:
		Maximum degree to calculate the model for 1 - 30 (default = 30).
		
	Returns
	=======
	Pnm : float
		Polynomial coefficients for each theta and each degree, shape
		(MaxDeg+1,MaxDeg+1,n), where n is the length of theta.
	dPnm : float
		Derivatives of Pnm, with the same shape.
	
	'''
	
	Pnm = np.zeros((Deg+1,Deg+1),dtype='float64') + np.nan
	dPnm = np.zeros((Deg+1,Deg+1),dtype='float64') + np.nan
	sintheta = np.sin(theta)
	costheta = np.cos(theta)
	Pnm[0,0] = 1.0
	Pnm[1,0] = costheta
	Pnm[1,1] = sintheta
	dPnm[0,0] = 0.0
	dPnm[1,0] = -sintheta
	dPnm[1,1] = costheta
	for n in range(2,Deg+1):
		n21 = 2*n - 1
		n1 = n - 1
		for m in range(0,n+1):
			m1 = m - 1
			if m >= n1:
				Pnm[n,m] = n21*sintheta*Pnm[n1,m1]
				dPnm[n,m] = n21*(costheta*Pnm[n1,m1] + sintheta*dPnm[n1,m1])
			else:
				nm1 = n + m - 1
				d1nm = 1.0/(n-m)
				n2 = n - 2
				Pnm[n,m] = d1nm*(costheta*n21*Pnm[n1,m] - nm1*Pnm[n2,m])
				dPnm[n,m] = d1nm*(n21*(costheta*dPnm[n1,m] - sintheta*Pnm[n1,m]) - nm1*dPnm[n2,m])
				

	return Pnm,dPnm


def _LegendreArr(theta,Deg=10):
	'''
	This lovely little routine is used to calculate the P(cos(theta))
	Legendre polynomials for the model. It starts off with P_00, P_10
	and P_11 - then recursively calculates the higher degree 
	coefficients.
	
	Inputs
	======
	theta : float
		Colatitude coordinate in radians.
	Deg:
		Maximum degree to calculate the model for 1 - 10 (default = 10).
		
	Returns
	=======
	Pnm : float
		Polynomial coefficients for each theta and each degree, shape
		(MaxDeg+1,MaxDeg+1,n), where n is the length of theta.
	dPnm : float
		Derivatives of Pnm, with the same shape.
	
	'''
	
	Pnm = np.zeros((Deg+1,Deg+1,theta.size),dtype='float64') + np.nan
	dPnm = np.zeros((Deg+1,Deg+1,theta.size),dtype='float64') + np.nan
	sintheta = np.sin(theta)
	costheta = np.cos(theta)
	Pnm[0,0,:] = 1.0
	Pnm[1,0,:] = costheta
	Pnm[1,1,:] = sintheta
	dPnm[0,0,:] = 0.0
	dPnm[1,0,:] = -sintheta
	dPnm[1,1,:] = costheta
	for n in range(2,Deg+1):
		n21 = 2*n - 1
		n1 = n - 1
		for m in range(0,n+1):
			m1 = m - 1
			if m >= n1:
				Pnm[n,m] = n21*sintheta*Pnm[n1,m1]
				dPnm[n,m] = n21*(costheta*Pnm[n1,m1] + sintheta*dPnm[n1,m1])
			else:
				nm1 = n + m - 1
				d1nm = 1.0/(n-m)
				n2 = n - 2
				Pnm[n,m] = d1nm*(costheta*n21*Pnm[n1,m] - nm1*Pnm[n2,m])
				dPnm[n,m] = d1nm*(n21*(costheta*dPnm[n1,m] - sintheta*Pnm[n1,m]) - nm1*dPnm[n2,m])
				

	return Pnm,dPnm
