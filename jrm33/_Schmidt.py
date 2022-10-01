import numpy as np
import Globals

def _Schmidt():
	'''
	Calculate the Schmidt normalization parameters for each n,m 
	combination.
	

	Returns
	=======
	Snm : float
		2D array of Schmidt normalization factors, shape 
		(MaxDeg+1,MaxDeg+1)
	
	
	'''
	if not Globals.Snm is None:
		return Globals.Snm
	
	
	#output matrix
	MaxDeg = Globals.MaxDeg
	Snm = np.zeros((MaxDeg+1,MaxDeg+1),dtype='float64') + np.nan
	
	#list a bunch of factorials from 0 to (n+m)!
	facts = np.append(1,np.cumprod(np.float64(np.arange(2*MaxDeg)+1)))

	#fill the output array
	for n in range(0,MaxDeg+1):
		for m in range(0,n+1):
			if m == 0:
				delta = 1
			else:
				delta = 2
			Snm[n,m] = np.sqrt(delta*((facts[n-m]/facts[n+m])))
			
	#save in globals
	Globals.Snm = Snm
	return Snm
