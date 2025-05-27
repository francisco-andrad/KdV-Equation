import numpy as np
import matplotlib.pyplot as plt

mass = []
x = []
x = np.linspace(0, 10, 200)
m = open('mass_data.txt', 'r')
mass_data = m.readlines()[0:200]
for i in range(200):
    mass.append(float(mass_data[i]))


plt.title("Variação da massa:")
plt.xlabel("tempo")
plt.ylabel("massa")
# plt.xlim(-100, 0)
plt.ylim(42.60, 42.80)
plt.plot(x, mass)
plt.show()
m.close()
