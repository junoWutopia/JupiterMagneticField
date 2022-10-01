import numpy as np
import Globals
import _ReadCoeffs

def _CoeffGrids():
	'''
	This function returns the grids of "g" and "h" coefficients for the 
	model. If they do not already exist in memory, they are read in from 
	file and stored in jrm09.Globals
	
	Returns
	=======
	g : float
		2D array of "g" coefficients, shape (MaxDeg+1,MaxDeg+1), where
		MaxDeg = 10 for JRM09.
	h : float
		2D array of "h" coefficients, shape (MaxDeg+1,MaxDeg+1), where
		MaxDeg = 10 for JRM09.
	
	'''
	#check if we have loaded them yet
	if not Globals.g is None and not Globals.h is None:
		return Globals.g,Globals.h
	
	
	#read the file with the coefficients in first
	if Globals.coeffs is None:
		Globals.coeffs = _ReadCoeffs._ReadCoeffs()
	coeffs = Globals.coeffs
	
	#create the grids (n,n) shape
	MaxDeg = np.max(coeffs.n)
	g = np.zeros((MaxDeg+1,MaxDeg+1),dtype='float64')
	h = np.zeros((MaxDeg+1,MaxDeg+1),dtype='float64')
		
	#fill it in
	g[coeffs.n,coeffs.m] = coeffs.g
	h[coeffs.n,coeffs.m] = coeffs.h

	#add to globals
	Globals.g = g
	Globals.h = h
	
	return g,h
	
	
