import numpy as np 
import matplotlib.pyplot as plt 

Datos = np.loadtxt('Datos.dat')

radio = Datos[:,0]
v_rot = Datos[:,1]

plt.plot(radio, v_rot)
plt.xlabel("Radio")
plt.ylabel("Velocidad de Rotacion")
plt.savefig("Datos.png")
#plt.show()
plt.clf()
