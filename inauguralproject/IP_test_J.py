#importing packages
from types import SimpleNamespace

import numpy as np
from scipy import optimize
from scipy import minimize_scalar

import pandas as pd
import matplotlib.pyplot as plt


#Defining the class
class household : 
    def __init__(self):
        

# given data with baseline parameters

x = np.linspace(0,24,49)

Lm = x
Lf = x 
Hm = x
Hf = x


# Baseline Parameters 
par= SimpleNamespace()
par.rho = 2.0
par.nu = 0.001 
par.epsilon = 1.0
par.omega = 0.5
par.alpha = 0.5
par.sigma = 1
par.omega_m = 1
par.omega_f = 1

alpha_dict = [0.25 , 0.50 , 0.75]
sigma_dict = [0.5 , 1.0 , 1.5]


if par.sigma == 0 :
    H= min(Hm,Hf)
if par.sigma == 1 : 
    H= Hm**(1-par.sigma) * Hf**par.sigma 
else : 
    H= ((1-par.alpha)*Hm**((par.sigma-1)/par.sigma) + par.alpha* Hf**((par.sigma-1)/par.sigma))**(par.sigma/(par.sigma-1))


C = par.omega_m * Lm + par.omega_f * Lf
Q = C**par.omega * H**(1-par.omega)
Tm = Lm + Hm 
Tf = Lf + Hf 
Lm = Hm = Lf = Hf >= 0 
Tm = Tf  <= 24 


# Defining the model
epsilon = par.epsilon
rho = par.rho
nu = par.nu
def u_func(Q,Tm,Tf,epsilon,rho,nu) : 
    return (Q^(1-rho))/(1-rho)-nu * (Tm^(1+ 1/(epsilon))/(1+ 1/epsilon) + Tf^(1+ 1/epsilon)/(1+ 1/epsilon))



