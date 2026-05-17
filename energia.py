import numpy as np
import matplotlib.pyplot as plt

energy = []
x = []
x = np.linspace(0, 10, 100)
e = open('nc_energy_data.txt', 'r')
energy_data = e.readlines()[0:100]
for i in range(100):
    energy.append(float(energy_data[i]))


derivada = np.gradient(energy,x) 
# plt.title("Diferença entre a posição dos centros de Energia e Massa:")
plt.xlabel("t")
plt.ylabel("energy")
plt.grid(True, alpha=0.6)
# plt.xlim(-100, 0)
# plt.ylim(25.0, 37.0)
plt.plot(x, energy)
plt.show()
e.close()
