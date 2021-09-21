### find ejected stars with 2 methods: 
#1) count all stars that cross R_evaporation as ejected
#2) conut stars outside R_evaporation at each time

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
	print('\n\n')
	time_=np.arange(tin, 28.75, 0.25)  #start, stop, step	
	time_str=[]
	for i in range(len(time_)):
		t_=format(time_[i], '.2f')
		fname1=folder+'single.40_'+t_ #check if the file exists
		myfile=Path(fname1)
		if myfile.is_file(): time_str.append(t_)
	return time_str


def escaped(t, Rescape, idxesc1):  #finds n° of ejected stars at time t

	nesc1, nesc2 = 0,0

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
	
	#check which particles escaped: 
	for i in range(len(idx)):   #loop over all stars

		if idx[i] in idxesc1 : ##METHOD 1: add previously escaped stars
			nesc1+=1
				
		elif r[i] > Rescape :  #add new stars	
			nesc1+=1
			idxesc1.append(idx[i])
			
		if r[i] > Rescape : ##METHOD 2: 	
			nesc2+=1
	
	return float(t), nesc1, nesc2, idxesc1

def nstars(color,Rev):
	time_str = files(0., 28.75, 0.25) #generate file names
		
	t=np.zeros(len(time_str)) 
	nesc1=np.zeros(len(time_str))  
	nesc2=np.zeros(len(time_str))  
	idxesc1=[]

	for i in range(len(time_str)):  #find Nstars for each file (at each time)
		t[i], nesc1[i], nesc2[i], idxesc1=  escaped(time_str[i], Rev , idxesc1)
		
	t*=tscale

	Nstars1 = (ntot-nesc1)/ntot
	Nstars2 = (ntot-nesc2)/ntot
	print('evaporation radius: ', Rev, ' pc')	
	print('fraction of stars:\nmethod 1 (all ejected): ', Nstars1[-1], '\nmethod 2 (inside the sphere at t): ', Nstars2[-1]   )		

	plt.plot(t, Nstars1, ls='--', color=color) #criterium 1: eliminate stars that cross R_evaporation
	plt.plot(t, Nstars2, color=color) #criterium 2: only stars outside R_evaporation at the moment
	return  Nstars2[-10]
	
	
	

folder='N5000_/'
ntot=5000

####### initial masses and indeces: 
'''
fname=folder+'single.40_0.00' #initial file with all stars
f = open(fname, "r")
index, mass = np.genfromtxt(f,dtype="float", comments="#", usecols=(0,1), unpack=True)
f.close()
index=index.astype(int)
idmass_dict = dict(zip(index, mass)) #dictionary with initial ids, masses
mtot=sum(mass)
print('total mass of the cluster: ', mtot)
print('average mass: ', mtot/len(mass))
'''

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

'''
t1,x1=lagr7(19) # t[NB]   19)Rcore[N.B.]
t2,x2=global30(8) #t[Myr], 8)Rcore[pc] 
tscale=t2[2]/t1[2]
lscale=x2[2]/x1[2]
print('lscale = ', lscale, '\ntscale = ', tscale)

###### check CM position:
t,x=global30(15) #N.B.
t=np.delete(t,0)
x=np.delete(x,0)
x=x*lscale #pc
print('max cluster CM distance from the origin: ', np.amax(x))

find lagrangian radii
tl,x = lagr7(15)   #15,16,17 
r1 = x[1]*lscale
print(x[0],' mass Lagr. radius')

tl,x = lagr7(16)   #15,16,17 
r2 = x[1]*lscale
print(x[0],' mass Lagr. radius')

tl,x = lagr7(17)   #15,16,17 
r3 = x[1]*lscale
print(x[0],' mass Lagr. radius')
'''

############################################################################
#					N ejected stars
#############################################################################

#7.5, 8, 9

fig=plt.figure(figsize=(8,6))

y1=nstars('orange',8. ) #r1[1]*lscale)   ### call all finctions from here
y2=nstars('red', 9.5 ) #r2[1]*lscale )
y3=nstars('blue', r3[1]*lscale) #r3[1]*lscale )

plt.text(18, y1-0.05, '90 % ', color='orange', fontsize=12)
plt.text(18, y2+0.02, '95 % ', color='red', fontsize=12)
plt.text(18, y3+0.02, '99 % ', color='blue', fontsize=12)
plt.xlabel('t [Myr]', fontsize=14)
plt.ylabel('N/N$_{tot}$', fontsize=14)
plt.title('Number of stars', weight='heavy')
plt.tight_layout()
plt.show()


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
	
	





#EVAPORATION RATE
dNesc90=np.zeros(len(Nesc90)-1)
tt=np.zeros(len(Nesc90)-1)
for i in range(len(Nesc90)-1):
	tt[i]=(t[i+1]+t[i])/2
	dNesc90[i]=(Nesc90[i+1]-Nesc90[i])#/(t[i+1]-t[i])

plt.plot(tt, dNesc90, color='green')
plt.show()





#CLUSTER MASS
t_,Mcore = global30(11) 

plt.plot(time, Mstars90)
#plt.plot(t_, Mcore)
plt.show()

'''




'''
plt.plot(time, mesc)
plt.xlabel('t [Myr]', fontsize=14)
plt.ylabel('$M_{esc}\, [M_\odot]$', fontsize=14)
plt.show()

plt.plot(time, Eint)
plt.xlabel('t [Myr]', fontsize=14)
plt.ylabel('$E_{int} [erg]$', fontsize=14)
plt.title('binary internal energy')
plt.show()
'''




'''
#chekk energy variation
t,de =global30(4) 
x=np.zeros(len(de)-1) 	#deltaE/E
for i in range(len(de)-1): 
	x[i]=(de[i+1]-de[i])/de[i]
t=t[:-1]
plt.plot(t,x)
#plt.yscale('log')
plt.ylabel('$\Delta E / E$')
plt.xlabel('t [Myr]')
plt.show()
'''



'''
t,m =status36(5)    
plt.plot(t,m)
plt.show()
'''
