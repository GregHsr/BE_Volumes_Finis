import matplotlib.pyplot as plt
import numpy as np

###------------------------------------------------------------------------###
###----------------------------- Initialisation ---------------------------###
###------------------------------------------------------------------------###

## Lire le fichier ligne par ligne
with open('sol_analyse.csv', 'r') as f:        # fichier d'analyse produit par le code Fortran
    lines = f.readlines()

## Données
N = 20                                         # A MODIFIER EN FONCTION DU NOMBRE DE MAILLES -> Nx = Ny = N
L = 1

###------------------------------------------------------------------------###
###-------------- Lecture des données (ne pas modifier) -------------------###
###------------------------------------------------------------------------###
x = np.linspace(0, L, N)
y = np.linspace(0, L, N)

## Extraire les données
temps = []
Cx_ana = []
Cx_num = []
Cy_ana = []
Cy_num = []
C_moy = []

print(len(lines))

for i in range(0, len(lines), 4*N+2):
    time = lines[i].rstrip()
    temps.append(float(time))
    for k in range(1,N+1):
        Cx_ana.append(float(lines[i+k].strip().split(" ")[0]))
        #print("k1",k)
    for k in range(1,N+1):
        Cx_num.append(float(lines[i+k+N].strip().split(" ")[0]))
        #print("k2",k)
    for k in range(1,N+1):
        Cy_ana.append(float(lines[i+k+2*N].strip().split(" ")[0]))
        #print("k3",k)
    for k in range(1,N+1):
        Cy_num.append(float(lines[i+k+3*N].strip().split(" ")[0]))
        #print("k4",k)
    #print(i)
    moy = lines[i+4*N+1].rstrip()
    C_moy.append(float(moy))
    #print("moy",moy)

###------------------------------------------------------------------------###
### ------------------Tracer les courbes (analyse) ------------------------###
###------------------------------------------------------------------------###

### Validation du programme

## Concentration en fonction de y
# k -> numéro de l'itération à tracer
# Cette partie trace C(y,k*dt) en x=L/2

plt.figure(1)
k = 2
print(Cx_num[int(k*N):int((k+1)*N)])
plt.plot(y, Cy_ana[int(k*N):int((k+1)*N)], label="Cy_ana")
plt.plot(y, Cy_num[int(k*N):int((k+1)*N)],'o',label="Cy_num")
plt.xlabel("y")
plt.ylabel("C_y")
plt.xlim(0,1)
plt.ylim(0,1)
plt.legend()
plt.grid()

## Concentration en fonction de y
# k -> numéro de l'itération à tracer
# Cette partie trace C(x,k*dt) en y=L/2

plt.figure(2)
plt.plot(x, Cx_ana[int(k*N):int((k+1)*N)], label="Cx_ana")
plt.plot(x, Cx_num[int(k*N):int((k+1)*N)],'o',label="Cx_num")
plt.xlabel("x")
plt.ylabel("C_x")
plt.xlim(0,1)
plt.ylim(0,1)
plt.legend()
plt.grid()


## Concentration en fonction de y
# On trace toutes les itérations (peu lisible)

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


## Concentration moyenne
# On trace C_moy en fonction du temps  (des itérations)
# Variables à modifier en fonction du cas :
#   - val_x : valeur ou C_moy atteint 0.45
#   - ligne 115 et 116 : valeur seuil et valeur x correspondante (pour affichage de la valeur sur la courbe)

val_x = 11.3
print(val_x)
plt.figure(5)
plt.plot(temps, C_moy)
plt.xlabel("temps [s]")
plt.ylabel("C_moy (Pe=50)")
plt.xlim(0,temps[-1])
plt.ylim(0,C_moy[-1]+0.05)
plt.grid()
# Dessiner une ligne horizontale et verticale là ou la courbe atteint 0.475 
lim = C_moy[-1]*0.90 # valeur seuil
print(lim)
plt.axhline(y=lim, color='r', linestyle='--')
plt.axvline(x=val_x, color='r', linestyle='--')
plt.text(0.01, lim + 0.01, '0.45', color='r')
plt.text(val_x, 0.01, '11.3', color='r')

## Temps adimensionnels (Complété à la main)

Pe = [0.5,5,50] 
D = [0.2,0.02,0.002]
alpha = [0.005,0.05,0.5]

td = [0.91,8.7,57.3]
tc = [18.1,17.3,11.3]

tad = [td[k]/(L*L/D[k]) for k in range(len(td))]
tac = [tc[k]/(L/alpha[k]) for k in range(len(tc))]

plt.figure(6)
plt.plot(Pe, tad, label="temps caractéristique de diffusion")
plt.plot(Pe, tac, label="temps caractéristique de convection")
plt.xlabel("Pe")
plt.ylabel("temps adimensionnels")
plt.legend()
plt.grid()
plt.loglog()


### Print et enregistrement manuel des valeurs à un t fixe (20s) pour la convergence en maillage ###
# Complété à la main

print(Cy_num[-N:])

# Pour N = 5
Cy_6 = [0.910114288, 0.651918232, 0.419149011, 0.323979318, 0.184287861, 0.0460528918]

# Pour N = 10
Cy_10 = [0.957710087, 0.848909199, 0.691412449, 0.52886045, 0.426742643, 0.372594804, 0.30068332, 0.193292707, 0.0935369283, 0.0260876156]

# Pour N = 20
Cy_20 = [0.981295407, 0.940559626, 0.890447617, 0.826979101, 0.750044286, 0.664939284, 0.581839681, 0.511778235, 0.460939765, 0.427513063, 0.403013617, 0.376344204, 0.339088738, 0.289105982, 0.230689853, 0.171529293, 0.118474558, 0.0749225691, 0.0406267494, 0.0127814431]

# Pour N = 30
Cy_30 = [0.987869203, 0.962610006, 0.934431612, 0.901657641, 0.863057494, 0.818086028, 0.767138004, 0.711736023, 0.654517949, 0.598893523, 0.548354983, 0.505619526, 0.471922547, 0.446728259, 0.42791599, 0.412318081, 0.396390319, 0.376997679, 0.352074444, 0.320978165, 0.284519792, 0.244638905, 0.203819439, 0.164453268, 0.1283613, 0.0965714902, 0.0693244711, 0.0462105796, 0.0263480712, 0.00854795706]

plt.figure(7)
plt.plot(np.linspace(0, 1, 6), Cy_6, label="Cy_6")
plt.plot(np.linspace(0, 1, 10), Cy_10, label="Cy_10")
plt.plot(np.linspace(0, 1, 20), Cy_20, label="Cy_20")
plt.plot(np.linspace(0, 1, 30), Cy_30, label="Cy_30")
plt.xlabel("y")
plt.ylabel("C(L/2,y)")
plt.legend()    
plt.grid()

plt.show()
