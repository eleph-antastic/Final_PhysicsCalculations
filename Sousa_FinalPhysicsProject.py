#!/usr/bin/python
# coding: utf-8

#Julia Sousa
#2345424
#jsousa@chapman.edu
#CPSC236-01
#Final: Assigned to Physics
#Stability Derivatives/Mach Number Plotting (Constant Gamma = 1.4, Varying Mach Number)

#Imports
import numpy as np
import matplotlib.pyplot as plt

#Known Data
#Given
sigma_0d = .70
sigma_0c = .912
sigma_0n = .955
eff_comb = .96
stdfuelLHV = 43500000 #J/kg
T_03val = 2400 #Kelvin

#Heat Properties
y = 1.4
cp = 1005 #J/kg
R = 287.15 #J/kgK


#Functions
def pStagRatio(machNum , gamma):
    return ( 1 + (gamma-1)/2 * machNum**2 ) ** (gamma/ (gamma-1) )

def tStagRatio(machNum , gamma):
    return ( 1 + (gamma-1)/2 * machNum**2 )

def tStagRatioInv(tRatio , gamma):
    return ( (tRatio-1)*2/(gamma-1) )**0.5

def pStagRatioInv(pRatio , gamma):
    return (((pRatio)**((gamma-1)/gamma)- 1) * (2/(gamma-1)))**0.5

def fuelRatio(T_03 , T_0a , eff_comb , lhv , c_p):
    return (T_03/T_0a -1) / (eff_comb * lhv/(c_p*T_0a) - T_03/T_0a)

#Brayton Cycle
# @A: Atmosphere with H = 10600 m (from tables)
M_a = np.arange(2.0 , 5.2 , 0.2)
T_a = 219.3652 * np.ones(M_a.size) #K
P_a = 24220 * np.ones(M_a.size) #Pa
a_a = (y*R*T_a)**0.5

# @0A: Atmosphere slowed down to Stagnation State
T_0a = np.ones(M_a.size)
P_0a = np.ones(M_a.size)

for i in range(M_a.size):
    T_0a[i] = tStagRatio(M_a[i] , y) * T_a[i]
    P_0a[i] = pStagRatio(M_a[i] , y) * P_a[i]
    
    
# @02: End of Diffuser with Pressure Loss
T_02 = T_0a.copy()
P_02 = sigma_0d * P_0a

# @ 03: End of Combustor Stage with Pressure Loss
T_03 = T_03val * np.ones(M_a.size)
P_03 = sigma_0c * P_02

# @ 05: Beginning of Nozzle Stage with Pressure Loss
T_05 = T_03.copy()
P_05 = sigma_0n * P_03


# @ 5: End of Nozzle Stage with Pressure Loss (assuming we go back to P_a)
M_5 = np.ones(M_a.size)
T_5 = np.ones(M_a.size)
P_5 = np.ones(M_a.size)

for i in range(M_a.size):
    M_5[i] = pStagRatioInv(P_05[i]/P_a[i] , y)
    
    T_5[i] = (1/tStagRatio(M_5[i] , y)) * T_05[i]
    P_5[i] = (1/pStagRatio(M_5[i] , y)) * P_05[i]

a_5 = (y*R*T_5)**0.5
    
#Processing
# Fuel Ratios
fuelToAirRatios = np.ones(M_a.size)
for i in range(M_a.size): fuelToAirRatios[i] = fuelRatio(T_03[i] , T_0a[i] , eff_comb , stdfuelLHV , cp)
    
# Propulsive Efficiencies
effPropulsion = np.ones(M_a.size)
for i in range(M_a.size): effPropulsion[i] = 2*(M_a[i]*a_a[i]) / (M_5[i]*a_5[i] + M_a[i]*a_a[i])
    
# Thermal Efficiencies
effThermal = np.ones(M_a.size)
for i in range(M_a.size): effThermal[i] = ( (1+fuelToAirRatios[i])*(M_5[i]*a_5[i])**2/2 - (M_a[i]*a_a[i])**2/2 )/(fuelToAirRatios[i] * stdfuelLHV)

# Overall Efficiencies
effOverall = effPropulsion*effThermal

#Thrust Specific Fuel Consumption
tsfc = (fuelToAirRatios/((fuelToAirRatios+1)*(M_5*a_5)-(M_a*a_a)))*3600 #kg/N/hr 

#Specific Thrust
f = fuelToAirRatios
spThrust = (1+f)*M_5*a_5 - M_a*a_a 

#Plotting and Printing
fig1 , ax1  = plt.subplots()
ax1.plot(M_a , fuelToAirRatios , label = "Fuel to Air Ratio")
ax1.plot(M_a , effPropulsion , label = "Propulsive Efficiency")
ax1.plot(M_a , effThermal , label = "Thermal Efficiency")
ax1.plot(M_a , effOverall , label = "Overall Efficiency")
ax1.plot(M_a , tsfc , label = "TSPC")
ax1.set_xlabel("Mach Number") 
ax1.legend()
fig1.show()

fig2 , ax2  = plt.subplots()
ax2.plot(M_a , spThrust , label = "Specific Thrust")
ax2.set_xlabel("Mach Number")
ax2.set_ylabel("N*sec/kg")
ax2.legend()
fig2.show()

print("Fuel Ratio @ M = 5 : " + str(fuelToAirRatios[-1]))
print("Propulsive Efficiency @ M = 5 : " + str(effPropulsion[-1]))
print("Thermal Efficiency @ M = 5 : " + str(effThermal[-1]))
print("Overall Efficiency @ M = 5 : " + str(effOverall[-1]))
print("Thrust Specific Fuel Comsumption @ M = 5 : " + str(tsfc[-1]) + " kg/N/hr")
print("Specific Thrust @ M = 5 : " + str(spThrust[-1]) + " N*sec/kg")

