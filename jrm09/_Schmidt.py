import numpy as np
from . import Globals

def _Schmidt(MaxDeg=10):
	'''
	Calculate the Schmidt normalization parameters for each n,m 
	combination.
	
	Inputs
	======
	MaxDeg : int
		Maximum degree of the model
		
	Returns
	=======
	Snm : float
		2D array of Schmidt normalization factors, shape 
		(MaxDeg+1,MaxDeg+1)
	
	
	'''
	if MaxDeg in list(Globals.Snm.keys()):
		return Globals.Snm[MaxDeg]
	
	#output matrix
	Snm = np.zeros((MaxDeg+1,MaxDeg+1),dtype='float64') + np.nan
	
	#list a bunch of factorials from 0 to (n+m)!
	facts = np.append(1,np.cumprod(np.arange(2*MaxDeg)+1,dtype=np.longlong))
	#print(facts)
	#fill the output array
	for n in range(0,MaxDeg+1):
		for m in range(0,n+1):
			if m == 0:
				delta = 1
			else:
				delta = 2
			#print(delta*((facts[n-m]/facts[n+m])))
			Snm[n,m] = np.sqrt(delta*((facts[n-m]/facts[n+m])))
			
	#save in globals
	Globals.Snm[MaxDeg] = Snm
	return Snm
