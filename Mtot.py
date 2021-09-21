#EVOLUTION OF THE TOAL MASS, GIVEN THE EVAPORATION RADIUS

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import Path

AU=1.495978707e11 #m
pc=3.0856775813e16 #m
Myr=60*60*24*365*1e6 #s
Msun=1.98847e30 #kg
J_to_erg=1e-7 #1J=1e7 erg
G=6.67408e-11 #m^3/kg/s^2
 
 #write time for file names
def files(tin, tfin, step):  
	time_=np.arange(tin, 28.75, 0.25)  #start, stop, step	
	time_str=[]
	for i in range(len(time_)):
		t_=format(time_[i], '.2f')
		fname1=folder+'single.40_'+t_ #check if the file exists
		myfile=Path(fname1)
		if myfile.is_file(): time_str.append(t_)
	return time_str


def escaped(t, Rev):  #finds n° of ejected stars and their total mass, at time t
	nesc = 0
	mesc=[] 

	fname1=folder+'single.40_'+str(t)
	fname2=folder+'binary.40_'+str(t)

	#########  singles  #############à
	f = open(fname1, "r")	
	idx,x,y,z= np.genfromtxt(f,dtype="float", comments="#", usecols=(0,2,3,4), unpack=True)
	f.close()
	
	r=np.sqrt(x**2+y**2+z**2)
	
	######### binaries ##############
	myfile=Path(fname2)  #read binary file if it exists
	if myfile.is_file():
		f = open(fname2, "r")
		idx1,idx2,xcm,ycm,zcm = np.genfromtxt(f,dtype="float", comments="#", usecols=(0,1,5,6,7), unpack=True)
		f.close()
		
		rcm=np.sqrt(xcm**2+ycm**2+zcm**2)
				
		idx=np.append(idx,(idx1, idx2)) #add binary quantities, if binary file is present
		r=np.append(r,(rcm,rcm)) 
	else: print('no binary file')

	#check that the total number of particles is = ntot
	#if len(idx)!= ntot : print('WARNING: total n° of particles < ', ntot)
	
	#check which particles escaped and add their mass: 
	m=np.zeros(len(idx))
	for i in range(len(idx)):    #loop over all stars
		
		m[i] = idmass_dict.get(idx[i])   # array with all the masses
		if r[i] > Rev : 	
			nesc+=1
			mesc.append(m[i])
	mesc_tot=sum(mesc)
	
	if t=='28.50': 

		density=False
				
		bins=np.logspace(-1,np.log10(150),50)
		plt.hist(m, bins=bins, density=density, histtype='step',  color='tab:blue', lw=2)
		plt.hist(mesc, bins=bins, density=density, histtype='step',  color='tab:red', lw=2)
		plt.xscale('log')
		plt.yscale('log')
		plt.xlabel('M [M$_\odot$]', fontsize=14)
		plt.ylabel('Number', fontsize=14)	
		plt.legend(['IMF', 'ejected stars'], fontsize=12)
		plt.title('mass function', weight='heavy', fontsize=14)
		plt.tight_layout()
		plt.show()
	
	return float(t), nesc, mesc_tot  

def nstars(Rev):
	time_str = files(0., 28.75, 0.25) #generate file names
		
	t=np.zeros(len(time_str)) 
	nesc=np.zeros(len(time_str)) #n° of escaped stars
	mesc=np.zeros(len(time_str)) #total escaped mass 

	for i in range(len(time_str)):  #find Nstars for each file (at each time)
		t[i], nesc[i], mesc[i]=  escaped(time_str[i], Rev )
		
	t*=tscale

	Nstars = (ntot-nesc)   #n° of remaining stars
	Mstars = (mtot-mesc)  		#total remaining mass

	print('\n\nevaporation radius: ', Rev, ' pc')	
	print('remaining stars:', Nstars[-1], '	', Nstars[-1]/ntot*100, ' %')
	print('remaining mass: ', Mstars[-1], '	', Mstars[-1]/mtot*100, ' %'   )		

	return t, Nstars, Mstars

	

folder='N5000_/'
ntot=5000

####### initial masses and indeces: 
fname=folder+'single.40_0.00' #initial file with all stars
f = open(fname, "r")
index, mass = np.genfromtxt(f,dtype="float", comments="#", usecols=(0,1), unpack=True)
f.close()
index=index.astype(int)
idmass_dict = dict(zip(index, mass)) #dictionary with initial ids, masses
mtot=sum(mass)
print('total mass of the cluster: ', mtot)
print('average mass: ', mtot/len(mass))

###### read Lagr.7 and global.30: 
fname=folder+'lagr.7'
f = open(fname, "r")	
t_l,r1,r2,r3 , RC_l= np.genfromtxt(f,dtype="float", comments="#", usecols=(0,15,16,17,19), unpack=True) #t[NB], 90,95,99%L.R.[NB], Rcore[NB]
f.close()

fname=folder+'global.30'
f = open(fname, "r")	
t_g,RC_g, ccm = np.genfromtxt(f,dtype="float", comments="#", usecols=(1,8,15), unpack=True) #t[Myr], Rcore[pc], Rcm[NB]
f.close()	

###### find tscale and lscale to convert to Physical units	(physical units = scale * N.B.units)
tscale=t_g[2]/t_l[2]
lscale=RC_g[2]/RC_l[2]
print('lscale = ', lscale, '\ntscale = ', tscale)

###### check CM position:
t_g=np.delete(t_g,0)
ccm=np.delete(ccm,0)
ccm=ccm*lscale #to [pc]
print('max cluster CM distance from the origin: ', np.amax(ccm))

####### Evaporation radii to test:
print('R_evaporation: \n', r1[0], ' Mtot: ',r1[1]*lscale,'\n', r2[0], ' Mtot: ',r2[1]*lscale, '\n', r3[0], ' Mtot: ',r3[1]*lscale)


############################################################################
#					PLOT
#############################################################################

t, Nstars, Mstars = nstars(r3[1]*lscale) 


#PLOT CLUSTER MASS
fig=plt.figure(figsize=(8,6))
plt.plot(t, Mstars, color='blue')
plt.xlabel('t [Myr]', fontsize=14)
plt.ylabel('M [M$_\odot$]', fontsize=14)
plt.title('Mass of the cluster', weight='heavy')
fig.tight_layout()
plt.show()



#PLOT Nstars and Mstars together
'''
fig, ax1 = plt.subplots()

color = 'tab:blue'
ax1.set_xlabel('t [Myr]')
ax1.set_ylabel('$N_{stars}$', color=color)
ax1.plot(time, Nstars, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:red'
ax2.set_ylabel('$M_{tot}$', color=color)  # we already handled the x-label with ax1
ax2.plot(time, Mstars, color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()	
'''	
	







