#!/usr/bin/env python
# -*- coding: utf-8 -*-
import scipy.integrate as integrate
import scipy.optimize as optimize
import numpy as np
import matplotlib.pyplot as plt

# Aqui definimos el factor de Hubble inverso y normalizado como funcion de z
def E(z,OmDE):
    return 1/np.sqrt((1-OmDE)*(1+z)**3+OmDE)

# Aqui calculamos la distancia luminica
def dl(z,OmDE,h):
    inte=integrate.quad(E,0,z,args=(OmDE))
    # Velocidad del sonido en km/s
    c = 299792.458
    # Factor de Hubble
    Ho = 100*h
    return c*(1+z)/Ho * inte[0]


# Funcion para calcular el logaritmo del likelihood
def loglike(params):
    OmDE = params[0]
    h = params[1]
    delta = []
    for i in zandmu:
# Ahora quiero calcular la diferencia entre el valor reportado y el calculado
        delta.append(5*np.log10(dl(i[0],OmDE,h))+25-i[1])
    delta = np.array(delta)

    chisquare=np.dot(delta,np.dot(np.linalg.inv(covariance),delta))
    return -chisquare/2

# Aqui calculamos la cadena de Markov
def markovchain(pasos, covparams, pasoinicial):
    chain=[pasoinicial]
    likechain=[loglike(chain[0])]
    chol = np.linalg.cholesky(covparams)
    aceptados = 0
    rechazados = 0

    for i in range(pasos):
        #newpoint = [np.random.normal(chain[i][0],ancho_OmDE),np.random.normal(chain[i][1],ancho_h)]
        # Intentando usar la matriz de covarianza para encontrar nuevos valores
        rand = np.random.normal(0.,1.,2)
        newpoint = chain[i] + np.dot(chol,rand)
        liketry = loglike(newpoint)
        if np.isnan(liketry) :
            print 'Paso algo raro'
            liketry = -1E50
        elif liketry > likechain[i]:
            accept = 1
        else:
            accept = np.exp(liketry - likechain[i])

        if accept >= np.random.uniform(0.,1.):
            chain.append(newpoint)
            likechain.append(liketry)
            aceptados += 1
        else:
            chain.append(chain[i])
            likechain.append(likechain[i])
            rechazados += 1

    print "Razon de aceptacion", aceptados,"/",pasos,   "=",float(aceptados)/float(pasos)

    return chain, likechain


################################3333
#
# Aquí empieza el programa
#
######################################3



# Inicializando los datos
global zandmu
global covariance
zandmu = np.loadtxt('data/SCPUnion2.1_mu_vs_z.txt', skiprows=5,usecols=(1,2))
covariance = np.loadtxt('data/SCPUnion2.1_covmat_sys.txt')

# Aqui quiero minimizar el (menos) logaritmo del likeligood
#minimo = optimize.maximize(loglike,[0.7,0.7],bounds=((0,1),(0.5,4)))
#print "Valores mas probables de Omega_De , h: "
#print minimo.x


startchain=[0.7,0.7]
ancho_OmDE = 0.01
ancho_h = 0.01
paramcovariance = [[ancho_OmDE**2, 0.],[0., ancho_h**2]]

pasos = 100
burnin = 10
pixels = 30

# Aquí calcula una primer cadena para calcular la matriz de covarianza

chain ,likechain = markovchain(burnin, paramcovariance, startchain)

# Calculando la matriz de covarianza inicial
omega =[]
ache = []
omegaporache = []
for line in chain:
        omega.append(line[0])
        ache.append(line[1])
        omegaporache.append(line[1]*line[0])
omega_prom = np.mean(omega)
ache_prom = np.mean(ache)
print "Valores esperados Omega_DE = ", omega_prom, ", h = ", ache_prom
sigma_omega = np.sqrt(np.mean(np.square(omega))-omega_prom**2)
sigma_ache = np.sqrt(np.mean(np.square(ache)) - ache_prom**2)
omegaporache_prom = np.mean(omegaporache)
print "Incertidumbres sigmaOmega_DE = ", sigma_omega, ", sigmah = ", sigma_ache
paramcovariance = [[sigma_omega**2, omegaporache_prom-omega_prom*ache_prom],[omegaporache_prom-omega_prom*ache_prom, sigma_ache**2]]


# Ahora sí, está es la buena
startchain = [omega_prom, ache_prom]
chain ,likechain = markovchain(pasos, paramcovariance, startchain)

omega =[]
ache = []
for line in chain:
        omega.append(line[0])
        ache.append(line[1])
#### What are we really using the counter i for? Only cause it counts the line number? 
#### Does the line counter start in 0 at the same time i=1?
# Parece que esto solo es necesario si se usa **2
#omega = np.array(omega)
#ache = np.array(ache)
omega_prom = np.mean(omega)
ache_prom = np.mean(ache)
print "Valores esperados Omega_DE = ", omega_prom, ", h = ", ache_prom
sigma_omega = np.sqrt(np.mean(np.square(omega))-omega_prom**2)
sigma_ache = np.sqrt(np.mean(np.square(ache)) - ache_prom**2)
omegaporache_prom = np.mean(omegaporache)
print "Incertidumbres sigmaOmega_DE = ", sigma_omega, ", sigmah = ", sigma_ache

#plt.scatter(omega, ache)
plt.plot(omega, ache,'-o')
plt.figure()
plt.plot(likechain,'ro')
plt.figure()
plt.hist(omega,bins=30)
plt.figure()
plt.hist(ache,bins=30)
plt.figure()

# Esto es un histograma 2-D
min_om = np.min(omega)
max_om = np.max(omega)
min_ache = np.min(ache)
max_ache = np.max(ache)
axis_om = np.linspace(min_om, max_om, pixels)
axis_ache = np.linspace(min_ache, max_ache, pixels)
# Primero creamos una matriz cuadrada
graph = np.zeros([pixels,pixels])
range_om = (max_om - min_om)/pixels
range_ache = (max_ache - min_ache)/pixels
axis_om, axis_ache = np.meshgrid(axis_om,axis_ache)
vector_om = min_ache
i = 0
for line in chain:
      if line[0] == max_om:
          om_index = pixels - 1
          ache_index=int((line[1] - min_ache)/range_ache)
          if line[1] == max_ache:
              ache_index = pixels -1
      elif line[1] == max_ache:
          ache_index = pixels - 1
          om_index=int((line[0] - min_om)/range_om)
          if line[0] == max_om:
              om_index = pixels -1
      else:
          om_index=int((line[0] - min_om)/range_om)
          ache_index=int((line[1] - min_ache)/range_ache)
      try:
          graph[om_index, ache_index] += 1
      except:
          print max_om,max_ache,om_index,ache_index,line
#plt.pcolormesh(axis_om,axis_ache,graph,shading='gouraud')
plt.pcolormesh(axis_om,axis_ache,graph)
plt.figure()

# Graficando los contornos
suma = 0.
N_max=np.max(graph)
prob=0.
while prob < 0.68:
    suma = 0
    suma_hist=0
    for line in graph:
        for pix in line:
            if pix>=N_max:
                suma += pix
    prob = suma/float(pasos)
    N_max -= 1
contour = np.zeros([pixels,pixels])
for i in range(pixels):
    for j in range(pixels):
        if graph[i,j] > N_max:
            contour[i,j] = 1



while prob < 0.95:
    suma = 0.
    for line in graph:
        for pix in line:
            if pix>=N_max:
                suma += pix
    prob = suma/float(pasos)
    N_max -= 1
for i in range(pixels):
    for j in range(pixels):
        if graph[i,j] > N_max:
            contour[i,j] += 1

plt.pcolormesh(axis_om,axis_ache,contour)
plt.figure()
plt.contour(contour,[0,1,2])
plt.show()

