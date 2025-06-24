import numpy as np
import matplotlib.pyplot as plt


x_axis = np.linspace(-200, 200, 20001)
massa_eixo = np.linspace(0, 100, 1000)
u_axis = []
mass = []
m = open('mass_data.txt', 'r')
mass_data = m.readlines()[0:1000]
# for k in range(1000):
#     mass.append(float(mass_data[k]))
a = 0
b = 20001

aux = []
aux_mass = []

for i in range(200):
    plt.title("Solução - KdV defocusing")
    plt.xlabel("x")
    plt.ylabel("-u")
    plt.xlim(-200, 200)
    plt.ylim(-5, 15)
    f = open('kdv_data.txt', 'r')
    aux = f.readlines()[a:b]
    for j in range(20001):
        u_axis.append(float(aux[j]))

    # for k in range(1000):
    #     aux_mass.append(float(mass_data[i]))
    #     k = k+1

    # plt.plot(massa_eixc, aux_mass)
    plt.plot(x_axis, u_axis)
    plt.savefig(f'images_kdv/{i:00003}', dpi=100, facecolor='white')
    plt.close()
    u_axis.clear()
    aux_mass.clear()
    a = a + 20002
    b = b + 20002
    f.close()


# for k in range(1000):
#     mass.append(float(mass_data[k]))
#
# print(len(mass))
# plt.title("Simulação na semirreta -  colisão com o limite do intervalo ")
# plt.xlabel("x")
# plt.ylabel("massa")
# # plt.xlim(-100, 0)
# # plt.ylim(-5, 15)
# plt.plot(massa_eixo, mass)
# plt.show()
m.close()
