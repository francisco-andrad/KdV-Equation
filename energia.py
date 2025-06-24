import numpy as np
import matplotlib.pyplot as plt

energy = []
x = []
x = np.linspace(0, 10, 100)
e = open('energy_data.txt', 'r')
energy_data = e.readlines()[0:100]
for i in range(100):
    energy.append(float(energy_data[i]))


plt.title("Variação da energia:")
plt.xlabel("tempo")
plt.ylabel("energia")
# plt.xlim(-100, 0)
# plt.ylim(40.2, 40.6)
plt.plot(x, energy)
plt.show()
e.close()
