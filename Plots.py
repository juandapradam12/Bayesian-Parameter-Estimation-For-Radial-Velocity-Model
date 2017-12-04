import numpy as np 
import matplotlib.pyplot as plt 

#### DATOS

Datos = np.loadtxt('Datos.dat')

R_exp = Datos[:,0]
V_Rad_exp = Datos[:,1]

#### AJUSTE 

Ajuste = np.loadtxt('Ajuste.dat')

#### GRAFICA

plt.plot(R_exp, V_Rad_exp)
plt.plot(R_exp, Ajuste)
plt.xlabel("Radio")
plt.ylabel("Velocidad de Rotacion")
plt.savefig("Datos.png")
plt.savefig('DatosyModelo.png')
plt.show()
plt.clf()








