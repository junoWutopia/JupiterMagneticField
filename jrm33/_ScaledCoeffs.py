import numpy as np
import Globals
import _CoeffGrids
import _Schmidt

def _ScaledCoeffs():
	'''
	Return a scaled set of g and h coefficients. These are scaled by
	begin multiplied by the Schmidt normalization coefficients.
	
	Returns
	=======
	gs : float
		g*Snm
	hs : float
		h*Snm
	
	'''

	#return the globals variables if they have already been defined
	if not Globals.gs is None and not Globals.hs is None:
		return Globals.gs,Globals.hs
		
	#otherwise calculate them
	g,h = _CoeffGrids._CoeffGrids()
	Snm = _Schmidt._Schmidt()
	
	Globals.gs = g*Snm
	Globals.hs = h*Snm
	
	return Globals.gs,Globals.hs
