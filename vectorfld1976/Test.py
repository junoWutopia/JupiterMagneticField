import numpy as np
import time

def vecfld1976(extshell_rad=1.5,layer=5,MaxDeg=10):
	try:
		from mpl_toolkits.axes_grid1 import make_axes_locatable
		from mayavi import mlab
	except:
		raise SystemError('This function requires "matplotlib" to be instaled')
	
	radlay=layer+1
	r = np.linspace(0.0,extshell_rad,radlay)
	lat = np.linspace(-np.pi/2,np.pi/2,181)
	lon = np.linspace(0,2*np.pi,361)
	latg,rg,long = np.meshgrid(lat,r,lon)

	rc = 0.5*(r[1:] + r[:-1])
	latc = 0.5*(lat[1:] + lat[:-1])
	lonc = 0.5*(lon[1:] + lon[:-1])
	longc,latgc = np.meshgrid(lonc,latc)
	longcr = longc
	latgcr = (np.pi/2 - latgc)

	latgc,rgc,longc = np.meshgrid(latc,rc,lonc)
	
	alpha, b_0, C, D, M = 0.7, 0.09, 10, 1, 4.2
	for i in range(0,radlay-1):
		X,Y,Z=rgc[i]*np.sin(np.pi/2-latgc[i])*np.cos(longc[i]),rgc[i]*np.sin(np.pi/2-latgc[i])*np.sin(longc[i]),rgc[i]*np.cos(np.pi/2-latgc[i])
		
		s=(X**2+Y**2)**0.5
		p=np.zeros_like(X)	#p=np.arctan(Y/X)
		for i in range(X.shape[0]):
			for j in range(X.shape[1]):
				if (X[i][j]>0 and Y[i][j]>0):
					p[i][j]=np.arctan(Y[i][j]/X[i][j])
				elif (X[i][j]<0 and Y[i][j]>0):
					p[i][j]=np.arctan(Y[i][j]/X[i][j])+np.pi
				elif (X[i][j]<0 and Y[i][j]<0):
					p[i][j]=np.arctan(Y[i][j]/X[i][j])+np.pi
				elif (X[i][j]>0 and Y[i][j]<0):
					p[i][j]=np.arctan(Y[i][j]/X[i][j])+2*np.pi
		z=Z
		r=(X**2+Y**2+Z**2)**0.5
		Bs=3*M*z*s/(r**5)+ b_0*np.tanh(z/D)/(D*s*(r**alpha))- alpha*b_0*z*(np.log(np.cosh(z/D))+C)/(s*(r**(alpha+2)))
		Bp=0-6.12*0.001*(np.e**(s/500))*s*Bs
		Bz=M*(2*(z**2)-s**2)/(r**5)
		mlab.quiver3d(X,Y,Z,Bs*np.cos(p)-Bp*np.sin(p),Bs*np.sin(p)+Bp*np.cos(p),Bz)
	mlab.show()

if __name__ == '__main__':
	vecfld1976()
	

