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
C_moy = []

for i in range(0, len(lines), 4*N+2):
    time = lines[i].rstrip()
    temps.append(float(time))
    for k in range(1,21):
        Cx_ana.append(float(lines[i+k].strip().split(" ")[0]))
        #print("k1",k)
    for k in range(1,21):
        Cx_num.append(float(lines[i+k+N].strip().split(" ")[0]))
        #print("k2",k)
    for k in range(1,21):
        Cy_ana.append(float(lines[i+k+2*N].strip().split(" ")[0]))
        #print("k3",k)
    for k in range(1,21):
        Cy_num.append(float(lines[i+k+3*N].strip().split(" ")[0]))
        #print("k4",k)
    #print(i)
    moy = lines[i+4*N+1].rstrip()
    C_moy.append(float(moy))
    print("moy",moy)

### Tracer les courbes ###

# Validation du programme
"""
plt.figure(1)
k = 2
print(Cx_num[int(k*20):int((k+1)*20)])
plt.plot(y, Cy_ana[int(k*20):int((k+1)*20)], label="Cy_ana")
plt.plot(y, Cy_num[int(k*20):int((k+1)*20)],'o',label="Cy_num")
plt.xlabel("y")
plt.ylabel("C_y")
plt.xlim(0,1)
plt.ylim(0,1)
plt.legend()
plt.grid()

plt.figure(2)
plt.plot(x, Cx_ana[int(k*20):int((k+1)*20)], label="Cx_ana")
plt.plot(x, Cx_num[int(k*20):int((k+1)*20)],'o',label="Cx_num")
plt.xlabel("x")
plt.ylabel("C_x")
plt.xlim(0,1)
plt.ylim(0,1)
plt.legend()
plt.grid()

plt.figure(3)
plt.plot(y, Cy_ana[0:20], label="Cy_ana")
plt.plot(y, Cy_num[0:20], label="Cy_num", marker='o')
plt.xlabel("y")
plt.ylabel("Cy")
plt.legend()
plt.grid()
"""
"""
plt.figure(4)
for k in range(len(temps)):
    #print(Cx_num[int(k*20):int((k+1)*20)])
    #print(Cx_ana[int(k*20):int((k+1)*20)])
    #print(temps[k])
    plt.plot(y, Cy_ana[int(k*20):int((k+1)*20)])
    plt.plot(y, Cy_num[int(k*20):int((k+1)*20)],'o')
plt.xlabel("y")
plt.ylabel("C_y")
plt.grid()
"""

# Concentration moyenne

plt.figure(5)
plt.plot(temps, C_moy)
plt.xlabel("temps [s]")
plt.ylabel("C_moy")
plt.xlim(0,temps[-1])
plt.ylim(0,C_moy[-1]+0.05)
plt.grid()
# Dessiner une ligne horizontale et verticale là ou la courbe atteint 0.475 
lim = C_moy[-1]*0.95 # valeur seuil
plt.axhline(y=lim, color='r', linestyle='--')
plt.axvline(x=79.98, color='r', linestyle='--')
plt.text(0.01, 0.48, 0.474, color='r')
plt.text(82, 0.01, '79.98', color='r')

plt.show()