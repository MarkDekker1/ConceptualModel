#Preambule
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from numpy.linalg import inv
matplotlib.style.use('ggplot')
from Class_Eisenman import *


#Model parameters
tmax=100
Dt=0.01

#Physical parameters in seconds
T0=10.+273.
S0=35.
t_r=25.*24.*3600.
alpha_t=10.**(-4.)
alpha_s=7.6*10.**(-4.)
H=4500
t_d=180.*365.*24.*3600.
qv=1.72*10.**(-4.)
I=4.*10.**(-6.)/(365.*24.*3600.)
k0=0.05/(365.*24.*3600.)
M=0.1/(365.*24.*3600.)
alpha=0.0003
Fs=100.
alpha_e=0.9
alpha_p=0.9
G=0.00001
StefBoltz=5.6*10.**(-8.)
Emissivity=0.1
V=100000.
SedRate=0.1/(365.*24.*3600.)
beta=k0
Period=20000.
phase=0.

#Initial conditions
T_eq_ocean_0=25.+273.
T_pole_ocean_0=5.+273.
T_eq_atm_0=20.+273.
T_pole_atm_0=-10.+273.
A_0=2.
S_eq_ocean_0=35.
S_pole_ocean_0=25.
P_0=0.01
Ca_0=0.01
Co_0=0.0001
CO3_0=0.001

#Dictionary
Parameter_initial = {
    'tmax':tmax,
    'Dt':Dt,
    'T0':T0,
    'S0':S0,
    't_r':t_r,
    'alpha_t':alpha_t,
    'alpha_s':alpha_s,
    'H':H,
    't_d':t_d,
    'q/v':qv,
    'I':I,
    'k0':k0,
    'M':M,
    'alpha':alpha,
    'alpha_e':alpha_e,
    'alpha_p':alpha_p,
    'Fs':Fs,
    'V':V,
    'StefBoltz':StefBoltz,
    'Emissivity':Emissivity,
    'G':G,
    'T_eq_ocean_0':T_eq_ocean_0,
    'T_pole_ocean_0':T_pole_ocean_0,
    'T_eq_atm_0':T_eq_atm_0,
    'T_pole_atm_0':T_pole_atm_0,
    'S_eq_ocean_0':T_eq_ocean_0,
    'S_pole_ocean_0':S_pole_ocean_0,
    'A_0':A_0,
    'P_0':P_0,
    'Ca_0':Ca_0,
    'Co_0':Co_0,
    'CO3_0':CO3_0,
    'SedRate':SedRate,
    'beta':Co_0,
    'Period':Period,
    'phase':phase
    }
    
#Initiate model
m1 = Eisenman(Parameter_initial,method='Upwind',initialvalue=0)
m1.integrateModel()