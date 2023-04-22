import matplotlib.pyplot as plt
import numpy as np

# Lire le fichier ligne par ligne
with open('sol_analyse.csv', 'r') as f:
    lines = f.readlines()

# Données
N = 20
L = 1

x = np.linspace(0, L, N)
y = np.linspace(0, L, N)

# Extraire les données
temps = []
Cx_ana = []
Cx_num = []
Cy_ana = []
Cy_num = []

for i in range(0, len(lines), 4*20+1):
    time = lines[i].rstrip()
    temps.append(float(time))
    for k in range(1,21):
        Cx_ana.append(float(lines[i+k].strip().split(" ")[0]))
    for k in range(1,21):
        Cx_num.append(float(lines[i+k].strip().split(" ")[0]))
    for k in range(1,21):
        Cy_ana.append(float(lines[i+k].strip().split(" ")[0]))
    for k in range(1,21):
        Cy_num.append(float(lines[i+k].strip().split(" ")[0]))

# Tracer les courbes
plt.figure(1)
k = 7
print(Cx_ana[int(k*20):int((k+1)*20)])
plt.plot(x, Cx_ana[int(k*20):int((k+1)*20)], label="Cx_ana")
plt.plot(x, Cx_num[int(k*20):int((k+1)*20)], label="Cx_num")
plt.xlabel("x")
plt.ylabel("Cx")
plt.legend()
plt.grid()

plt.figure(2)
plt.plot(y, Cy_ana[0:20], label="Cy_ana")
plt.plot(y, Cy_num[0:20], label="Cy_num", linestyle='dotted')
plt.xlabel("y")
plt.ylabel("Cy")
plt.legend()
plt.grid()

plt.figure(3)
for k in range(len(temps)):
    plt.plot(x, Cx_ana[int(k*20):int((k+1)*20)])
plt.xlabel("x")
plt.ylabel("Cx")
plt.grid()

plt.show()