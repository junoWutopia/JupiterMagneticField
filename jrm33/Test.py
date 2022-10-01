import numpy as np
import Model
import time
import pstats
import cProfile

def map2d(R=0.85,MaxDeg=13):
	'''
	This is a simple function to test the model by recreating a plot in 
	Connerney et al 2018 (figure 4, sort of).
	
	Inputs
	======
	R : float
		The radial distance to evaluate the model at.
	MaxDeg : int
		Maximum model degree to calculate.
	
	'''
	try:
		import matplotlib.pyplot as plt
		import matplotlib.colors as colors
		from mpl_toolkits.axes_grid1 import make_axes_locatable
	except:
		raise SystemError('This function requires "matplotlib" to be instaled')

	#get the coordinates to calculate the model at
	lat = np.linspace(-90,90,181)
	lon = np.linspace(0.0,360.0,361)
	latc = 0.5*(lat[1:] + lat[:-1])
	lonc = 0.5*(lon[1:] + lon[:-1])
	
	long,latg = np.meshgrid(lon,lat)
	longc,latgc = np.meshgrid(lonc,latc)
	longcr = longc*np.pi/180.0
	latgcr = (90.0 - latgc)*np.pi/180.0
	r = np.zeros(longcr.shape) + R
	

	#calculate the model
	Br,Bt,Bp = Model.Model(r,latgcr,longcr,MaxDeg)
	
	#B = np.sqrt(Br**2 + Bt**2 + Bp**2)
	
	#convert to Gauss
	Bg = Br.reshape(longcr.shape)*1e-5


	#get the scale
	scale = [-60.0,60.0]
	
	#set norm
	norm = colors.Normalize(vmin=scale[0],vmax=scale[1])	
		
	maps = [1,1,0,0]
	fig1=plt.figure()
	ax1 = plt.subplot2grid((maps[1],maps[0]),(maps[3],maps[2]))
	ax1.set_aspect(1.0)
	ax1.set_xlabel('SIII East Longitude ($^\circ$)')
	ax1.set_ylabel('SIII Latitude ($^\circ$)')

	sm = ax1.pcolormesh(long,latg,Bg,cmap='RdYlBu_r',norm=norm)
	ct = ax1.contour(longc,latgc,Bg,colors='grey',levels=np.linspace(-50,50,11))
	ax1.clabel(ct, inline=True, fontsize=8,fmt='%2d')

	divider = make_axes_locatable(ax1)
	cax1 = divider.append_axes("right", size="5%", pad=0.05)

	cbar1 = plt.colorbar(sm,cax=cax1) 
	cbar1.set_label('$B_r$ (Gauss) at $r$ = {:4.2f}'.format(R) + ' R$_{j}$')

	fig1.savefig('jrm33test.png')
	
	
	Ns = np.linspace(90, 0, 91) 
	Ntheta = np.linspace(0, 2*np.pi, 361) 
	Nsc = 0.5*(Ns[1:] + Ns[:-1])
	Nthetac = 0.5*(Ntheta[1:] + Ntheta[:-1])
	NTheta, NS = np.meshgrid(Ntheta, Ns)
	NThetac, NSc = np.meshgrid(Nthetac, Nsc)

	fig2 = plt.figure()
	ax2 = plt.subplot2grid((maps[1],maps[0]),(maps[3],maps[2]))
	ax2.set_aspect(1.0)
	ax2.axes.xaxis.set_ticklabels([])
	ax2.axes.yaxis.set_ticklabels([])

	BgNpole=Bg[90:]
	BgNpole=BgNpole[::-1]
	sm = ax2.pcolormesh((90-NS)*np.cos(NTheta-np.pi/2),(90-NS)*np.sin(NTheta-np.pi/2),BgNpole,cmap='RdYlBu_r',norm=norm)
	ct = ax2.contour((90-NSc)*np.cos(NThetac-np.pi/2),(90-NSc)*np.sin(NThetac-np.pi/2),BgNpole,colors='grey',levels=np.linspace(-50,50,11))
	ax2.clabel(ct, inline=True, fontsize=8,fmt='%2d')

	for i in range(4):
		ax2.add_patch(plt.Circle((0, 0), 90-30*i,fill=0,color='k',linewidth=0.6))
		ax2.annotate(str(30*i)+'°N',(0,90-30*i),verticalalignment='bottom',horizontalalignment='right')
	for j in range(4):
		ax2.plot([0,90*np.cos(np.pi/2*j)],[0,90*np.sin(np.pi/2*j)],color='k',linewidth=0.6)
		ax2.annotate(str(90*j)+'°E',(90*np.cos(np.pi/2*(j-1)),90*np.sin(np.pi/2*(j-1))),verticalalignment='top',horizontalalignment='right')

	divider = make_axes_locatable(ax2)
	cax2 = divider.append_axes("right", size="5%", pad=0.05)
	cbar2 = plt.colorbar(sm,cax=cax2) 
	cbar2.set_label('$B_r$ (Gauss) at $r$ = {:4.2f}'.format(R) + ' R$_{j}$')
	
	fig2.savefig('jrm33testNpolar.png')
	

	Ss = np.linspace(-90, 0, 91) 
	Stheta = np.linspace(0, 2*np.pi, 361) 
	Ssc = 0.5*(Ss[1:] + Ss[:-1])
	Sthetac = 0.5*(Stheta[1:] + Stheta[:-1])
	STheta, SS = np.meshgrid(Stheta, Ss)
	SThetac, SSc = np.meshgrid(Sthetac, Ssc)

	fig3 = plt.figure()
	ax3 = plt.subplot2grid((maps[1],maps[0]),(maps[3],maps[2]))
	ax3.set_aspect(1.0)
	ax3.axes.xaxis.set_ticklabels([])
	ax3.axes.yaxis.set_ticklabels([])


	BgSpole=Bg[:90]
	sm = ax3.pcolormesh((90+SS)*np.cos(-STheta+np.pi/2),(90+SS)*np.sin(-STheta+np.pi/2),BgSpole,cmap='RdYlBu_r',norm=norm)
	ct = ax3.contour((90+SSc)*np.cos(-SThetac+np.pi/2),(90+SSc)*np.sin(-SThetac+np.pi/2),BgSpole,colors='grey',levels=np.linspace(-50,50,11))
	ax3.clabel(ct, inline=True, fontsize=8,fmt='%2d')
	
	for i in range(4):
		ax3.add_patch(plt.Circle((0, 0), 90-30*i,fill=0,color='k',linewidth=0.6))
		ax3.annotate(str(90-30*i)+'°S',(0,30*i),verticalalignment='bottom',horizontalalignment='right')
	for j in range(4):
		ax3.plot([0,90*np.cos(np.pi/2*j)],[0,90*np.sin(np.pi/2*j)],color='k',linewidth=0.6)
		ax3.annotate(str(90*j)+'°E',(90*np.cos(np.pi/2*(j+1)),90*np.sin(np.pi/2*(j+1))),verticalalignment='top',horizontalalignment='right')

	divider = make_axes_locatable(ax3)
	cax3 = divider.append_axes("right", size="5%", pad=0.05)
	cbar3 = plt.colorbar(sm,cax=cax3) 
	cbar3.set_label('$B_r$ (Gauss) at $r$ = {:4.2f}'.format(R) + ' R$_{j}$')
	
	fig3.savefig('jrm33testSpolar.png')

	return ax1,ax2,ax3

def vecfld(extshell_rad=1.0,layer=3,MaxDeg=10):
	try:
		import matplotlib.pyplot as plt
		import matplotlib.colors as colors
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

	fig = plt.figure()
	maps = [1,1,0,0]
	ax = plt.subplot2grid((maps[1],maps[0]),(maps[3],maps[2]))
	ax.set_aspect(1.0)
	ax.axes.xaxis.set_ticklabels([])
	ax.axes.yaxis.set_ticklabels([])
	
	for i in range(0,radlay-1):
		#print(r[i+1])
		Br,Bt,Bp = Model.Model(np.zeros(longcr.shape) + r[i+1],latgcr,longcr,MaxDeg)
		Br = Br.reshape(longcr.shape)*1e-5
		Bt = Bt.reshape(longcr.shape)*1e-5
		Bp = Bp.reshape(longcr.shape)*1e-5
		X,Y,Z=rgc[i]*np.sin(np.pi/2-latgc[i])*np.cos(longc[i]),rgc[i]*np.sin(np.pi/2-latgc[i])*np.sin(longc[i]),rgc[i]*np.cos(np.pi/2-latgc[i])
		
		mlab.quiver3d(X,Y,Z,
			(Br*np.sin(np.pi/2-latgc[i])*np.cos(longc[i])+Bt*np.cos(np.pi/2-latgc[i])*np.cos(longc[i])-Bp*np.sin(longc[i])),
			(Br*np.sin(np.pi/2-latgc[i])*np.sin(longc[i])+Bt*np.cos(np.pi/2-latgc[i])*np.sin(longc[i])+Bp*np.cos(longc[i])),
			(Br*np.cos(np.pi/2-latgc[i])-Bt*np.sin(latgc[i]))
			)
	
	mlab.show()

if __name__ == '__main__':
	map2d()
	vecfld()