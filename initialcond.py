#GENERATE CLUSTER INITIAL CONDITIONS (MONTE CARLO)
import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt

N=int(5e3)

#/////////////////////  MASSES: Salpeter /////////////////////////////////////////
#parameters for the Salpeter mass function:
m_min, m_max = 0.1, 150. #masses in M_sun
a=2.3
p=-1.3  #1-a

#normalize to 1
I = lambda m : (m**p-m_min**p)/p  #cumulative distribution P(m)   (not notmalized to 1)
if I(m_max)!=1 : C = I(m_max) 	
else: C=1 	

#generate the random sample
#npr.seed(123) 		#set seed for reproducibility
x = npr.rand(N) 	# uniform distr [0,1]
m = ( C*p*x + m_min**p )**(1/p)

Mtot = sum(m)
print("total mass of the cluster = ", Mtot)
print("average mass", Mtot/N)


#////////////////////////   VELOCITIES: Maxwellian   //////////////////////////////// 
c, s = 0., 0.5  #parameters for the MaAxwellian: mean, sigma (km/s)

#[BoxMuller] generate 3 gaussian variables
P = np.random.rand(N)	#cumulative distribution of r, theta: uniform [0,1]
r=np.sqrt(-2*np.log(1-P)*s**2)	#find r, theta with inverse random sampling

P = np.random.rand(N)  
theta=2*np.pi*P

P = np.random.rand(N)	 #cumulative distribution of r, theta: uniform [0,1]
r2=np.sqrt(-2*np.log(1-P)*s**2) #find r, theta with inverse random sampling

P = np.random.rand(N)  
theta2=2*np.pi*P

vx=r*np.cos(theta)	#gaussian distributed variables
vy=r*np.sin(theta) 
vz=r2*np.cos(theta2)

v = np.sqrt(vx**2 + vy**2 + vz**2)  #Maxwellian distributed variable (module of the vector)

#geneate theta, phi UNIFORM (to get the direction of the vector)
P = np.random.rand(N)
theta = np.arccos(1-2*P)

P = np.random.rand(N)
phi = 2*np.pi*P

#convert to cartesian coordinates
vx = v*np.sin(theta)*np.cos(phi)	
vy = v*np.sin(theta)*np.sin(phi)
vz = v*np.cos(theta)

#///////////////////////   POSITIONS: Plummer sphere density    /////////////////////////////// 
a = 1. # [pc]	system's lenght scale (parameter for the Plummer)
P = np.random.rand(N)	
r = np.sqrt( a**2/(P**(-2/3)-1) ) 	

P = np.random.rand(N)	#angular coord: UNIFORM
theta = np.arccos(1-2*P)

P = np.random.rand(N)
phi = 2*np.pi*P

x = r*np.sin(theta)*np.cos(phi) 	#conversion to cartesian coord
y = r*np.sin(theta)*np.sin(phi)
z = r*np.cos(theta) 


#///////////////// VIRIALIZATION //////////////////////
#conversion factors:    
pc_m = 3.086e16		# 1 pc = 3.086e16 m 
pc_km = 3.086e13	# 1 pc = 3.086e13 km 
Msun_kg = 1.989e30	# M_sun = 1.989e30 kg 

G = 6.67408e-11 	 #m^3 kg^-1 s^-2
G *= (Msun_kg/(pc_m**3)) #pc^3 Msun^-1 s^-2  checked value: OK

#kinetic energy
v /= pc_km 
K = 0.5*sum(m*v**2)
print("kinetic energy", K)

#gravitational potential
U=0.
for i in range(len(m)):
	s=0.
	for j in range (i+1, len(m) ):
		s += (m[j] / np.sqrt((x[i]-x[j])**2 + (y[i]-y[j])**2 + (z[i]-z[j])**2 ) )
	U += m[i]*s
U *= G	
print("potential energy", U )


Qvir = 2*K/U	#virial ratio
print("Qvir", Qvir)

vx_vir = vx / np.sqrt(Qvir) #IMPOSE VIRIALIZATION
vy_vir = vy / np.sqrt(Qvir)
vz_vir = vz / np.sqrt(Qvir)

#check if Q_vir=1
v_vir = np.sqrt(vx_vir**2 + vy_vir**2 + vz_vir**2)
v_vir /= pc_km
K_vir = 0.5*sum(m*v_vir**2)
Qvir_ = 2*K_vir/U
print("kinetic energy after virialization", K_vir)
print("Qvir after virialization", Qvir_)

#WRITE TO FILE
# col.0: masses/msun, col.1: x/pc, col.2: y/pc, col.3: z/pc, col. 4:vx/kms, col. 5:vy/kms, col.6: vz/kms
f = open("initialcond.txt", "w")
for i in range (N):
	f.write(str(m[i])+' '+str(x[i])+' '+str(y[i])+' '+str(z[i])+' '+str(vx_vir[i])+' '+str(vy_vir[i])+' '+str(vz_vir[i])+'\n')

f.close()


ff = open("initialcond_prameters", "w")
ff.write('#	'+'Nbodies		'+'Mtot[M_sun]		'+'Qvir'+'\n' )
ff.write('	'+str(N)+' 	'+ str(Mtot)+' 	'+ str(Qvir)  )
ff.close()




