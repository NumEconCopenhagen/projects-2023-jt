#importing packages
from types import SimpleNamespace

import numpy as np
from scipy import optimize

import pandas as pd
import matplotlib.pyplot as plt

#Defining the class
class householdclass:
    def __init__(self):
        
        # creating namespaces
        par = self.par = SimpleNamespace()
        sol = self.sol = SimpleNamespace()

        # preferences
        par.rho = 2.0
        par.nu = 0.001
        par.epsilon = 1.0
        par.omega = 0.5

        # household production
        par.alpha = 0.5
        par.sigma = 1.0

        # wages
        par.wM = 1.0
        par.wF = 1.0
        par.wF_vec = np.linspace(0.8,1.2,5)

        #targets
        par.beta0_target = 0.4
        par.beta1_target = -0.1


        #solution
        sol.LM_vec = np.zeros(par.wF_vec.size)
        sol.HM_vec = np.zeros(par.wF_vec.size)
        sol.LF_vec = np.zeros(par.wF_vec.size)
        sol.HF_vec = np.zeros(par.wF_vec.size)

        sol.beta0 = np.nan
        sol.beta1 = np.nan
    
    def calc_utility(self,LM,HM,LF,HF):
        """calculate utility"""

        par = self.par
        sol = self.sol

        # consumption of market goods
        C = par.wM*LM + par.wF*LF

        # home production with different sigmas
        if par.sigma == 1:
            H = HM**(1-par.alpha)*HF**par.alpha
        
        if par.sigma == 0:
            min(HM, HF)
        
        else:
            ((1-par.alpha)*HM**((par.sigma-1)/par.sigma)+par.alpha*HF**((par.sigma-1)/par.sigma))**(par.sigma/(par.sigma-1))
        
        # total consumption utility
        Q = C**par.omega*H**(1-par.omega)
        utility = np.fmax(Q,1e-8)**(1-par.rho)/(1-par.rho)

        #disutility of work
        epsilon_ = 1+1/par.epsilon
        TM = LM+HM
        TF = LF+HF
        disutility = par.nu*(TM**epsilon_/epsilon_+TF**epsilon_/epsilon_)

        return utility - disutility
    
    def solve_discrete(self,do_print=False):
        """solve model discretely"""

        par = self.par
        sol = self.sol
        opt = SimpleNamespace()

        # All possible choices
        x = np.linspace(0,24,49)
        LM,HM,LF,HF = np.meshgrid(x,x,x,x) # All combinations

        LM = LM.ravel()
        HM = HM.ravel()
        LF = LF.ravel()
        HF = HF.ravel()

        # Calculate utility
        u = self.calc_utility(LM,HM,LF,HF)
        
        # Utility minus infinity if constraint is broken
        I = (LM+HM > 24) | (LF+HF > 24)
        u[I] = -np.inf

        # find maximizing argument
        j = np.argmax(u)

        opt.LM = LM[j]
        opt.HM = HM[j]
        opt.LF = LF[j]
        opt.HF = HF[j]

        # Print
        if do_print:
            for k,v in opt.__dict__.items():
                print(f'{k} = {v:6.4f}')
        return opt



