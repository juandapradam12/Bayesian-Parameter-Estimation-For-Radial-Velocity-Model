import numpy as np 
import matplotlib.pyplot as plt 

Datos = np.loadtxt('Datos.dat')

plt.plot(Datos[:,0], Datos[:,1])
plt.savefig("Datos.png")
#plt.show()
plt.clf()
