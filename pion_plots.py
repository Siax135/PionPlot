import matplotlib.pyplot as plt  # Used for managing and creating plots 
import numpy as np  # Used for its math library

MPI=140 # Mass of the Pion in MeV
FMTOMEV=197 # Factor for conversion from fermi to MeV

def normsq(n1,n2,n3):
	"""
	normsq(n1,n2,n3)

	n1,n2,n3 are integers that represent quanta of momentum in the x, y, and z directions respectfully.

	returns the magnitude squared on the vector {n1,n2,n3}.
	"""
	return (n1**2)+(n2**2)+(n3**2)

def pmagsq(n1,n2,n3,L):
	"""
	pmagsq(n1,n2,n3,L)

	n1,n2,n3 are integers that represent quanta of momentum in the x, y, and z directions respectfully.
	L is the length of one side of the volume the pions are in. (Assumes they are in a cubic volume)

	returns the square of the momentum for a single pion in MeV^2
	"""
	return ((2*np.pi/L)**2)*normsq(n1,n2,n3)*(FMTOMEV**2)

def energy(n1,n2,n3,L):
	"""
	energy(n1,n2,n3,L)

	n1,n2,n3 are integers that represent quanta of momentum in the x, y, and z directions respectfully.
	L is the length of one side of the volume the pions are in. (Assumes they are in a cubic volume)

	returns the energy for a single pion in cube of size L in the lab frame.	
	"""
	return np.sqrt(MPI**2 + pmagsq(n1,n2,n3,L))

def cmenergy(n1,n2,n3,d1,d2,d3,L):
	"""
	cmenergy(n1,n2,n3,d1,d2,d3,L)

	n1,n2,n3 are integers that represent quanta of momentum in the x, y, and z directions respectfully.
	d1,d2,d3 are integers that represent boosts in the x, y, and z directions respectfully.
	L is the length of one side of the volume the pions are in. (Assumes they are in a cubic volume)

	returns the center of mass energy for two pions in cube of size L.
	"""
	return np.sqrt(((energy(n1,n2,n3,L)+energy(d1-n1,d2-n2,d3-n3,L))**2) - pmagsq(d1,d2,d3,L))

Lrange = np.arange(6., 30., 0.2) # x-axis range with steps of 0.2
ecmValues=[] # empty list, used to make sure energy levels are double plotted
dValue=input("Input magnitude for d: ") # ask for input (integer) from the user

# add threshold to plot
plt.plot(Lrange,cmenergy(0,0,0,0,0,0,Lrange),'r--',label="threshold")

# nested for loops that run over possible boosts and momentums for two pions in a cubic volume of size L
for d1 in range(0,3):
	for d2 in range(0,3):
		for d3 in range(0,3):
			if normsq(d1,d2,d3) == dValue: 
				for n1 in range(1,4):
					for n2 in range(0,3):
						for n3 in range(0,3):
							ecm6=cmenergy(n1,n2,n3,d1,d2,d3,6)
							if ecm6 not in ecmValues:
								ecmValues.append(ecm6) # add energy values to list if they aren't already in it. This prevents double plotting
								plt.plot(Lrange,cmenergy(n1,n2,n3,d1,d2,d3,Lrange),label="|n|^2="+str(normsq(n1,n2,n3)))


# set axis labels and title
plt.xlabel("L (fm)")
plt.ylabel("Energy (MeV)")
plt.title("Energy Spectrum |d|^2="+str(dValue))

# uncomment line below to show plot legend
#plt.legend()

# show plot
plt.show()