from types import SimpleNamespace

from scipy import optimize
import numpy as np
import sympy as sm
from scipy.optimize import fsolve

import pandas as pd 
import matplotlib.pyplot as plt
import ipywidgets 


class SolowModelClass:
    def __init__(self):
        
        # We create the namespaces: 
        par = self.par = SimpleNamespace()
        sim = self.sim =SimpleNamespace()
        
        # Sim parameters
        sim.alpha = 0.333
        sim.phi = 0.333
        sim.delta = 0.02
        sim.n = 0.014
        sim.g = 0.016
        sim.s_K = 0.25
        sim.s_H = 0.129



        # We name our parameters:
        par.alpha = sm.symbols('alpha')
        par.phi  = sm.symbols('phi')
        par.delta = sm.symbols('delta')
        par.n = sm.symbols('n')
        par.g = sm.symbols('g')
        par.s_K = sm.symbols('s_{K}')
        par.s_H = sm.symbols('s_{H}')


        # We name our variables

        par.K_t = sm.symbols('K_{t}')
        par.K_t1 = sm.symbols('K_{t+1}')
        par.H_t = sm.symbols('H_{t}')
        par.H_t1 = sm.symbols('H_{t+1}')
        par.Y_t = sm.symbols('Y_{t}')
        par.A_t = sm.symbols('A_{t}')
        par.A_t1 = sm.symbols('A_{t+1}')
        par.L_t = sm.symbols('L_{t}')
        par.L_t1 = sm.symbols('L_{t+1}')
        par.r_t = sm.symbols('r_{t}')
        par.w_t = sm.symbols('w_{t}')
        
        # We name our per effective worker variables
        par.ktilde_t = sm.symbols('\tilde{k_{t}}')
        par.ktilde_t1 = sm.symbols('\tilde{k}_{t+1}')
        par.htilde_t = sm.symbols('\tilde{h_{t}}')
        par.htilde_t1 = sm.symbols('\tilde{h}_{t+1}')
        par.ytilde_t = sm.symbols('\tilde{y_{t}}')


    # we will define our functions as given by the book:
    def ProductionFunction(self):
        par = self.par
        Production = sm.Eq(par.Y_t,par.K_t**par.alpha * par.H_t**par.phi * (par.A_t * par.L_t)^(1-par.alpha - par.phi))

        return

    def LabourSupply(self):
        par = self.par
        Labor_supply = sm.Eq(par.L_t1,(1+par.n)*par.L_t )

        return

    def LabourProductivity(self):
        par = self.par
        Productivity = sm.Eq(par.A_t1 , (1+ par.g)*par.A_t)

        return

    def RentalrateOfCapital(self):
        par = self.par
        RentalRate = sm.Eq(par.r_t,par.alpha * (par.K_t/(par.A_t * par.L_t))^(par.alpha -1) * (par.H_t/(par.A_t * par.L_t))^par.phi)

        return

    def RealWageRate(self):
        par = self.par
        WageRate = sm.Eq(par.w_t, par.alpha * (par.K_t/(par.A_t * par.L_t))^par.alpha * (par.H_t/(par.A_t * par.L_t))^par.phi * par.A_t)

        return

    def HumanCapitalAccumulation(self):
        par = self.par
        HumanCapital = sm.Eq(par.H_t1 , par.s_H * par.Y_t + par.delta * par.H_t)

        return

    def PhysicalCapitalAccumulation(self):
        par = self.par
        PhysicalCapital = sm.Eq(par.K_t1 , par.s_K * par.Y_t + par.delta * par.K_t)

        return    

    def SteadyState(self , do_print=False):
        par = self.par
        alpha = par.alpha
        phi = par.phi
        s_K = par.s_K
        s_H = par.s_H
        delta = par.delta
        n = par.n
        g = par.g

        # In steady state we know that \tilde{k}_{t+1} = \tilde{k}_{t} = k
        # and that \tilde{h}_{t+1} = \tilde{h}_{t} = h
        # this we are defining now
        k = par.ktilde_t = par.ktilde_t1 
        h = par.htilde_t = par.htilde_t1 

        y =  k**alpha * h**phi

        # We define the function for which we are calculating the ss-value 
        ss_k = sm.Eq(k, 1/((1+n)*(1+g)) * (s_K)*y + (1-delta)*k ) 
        # We find the steady state for k, by putting the lef hand side equal to 0
        kss = sm.solve(ss_k,k)[0]
        
        # We will now do the same for h
        ss_h = sm.Eq(h, 1/((1+n)*(1+g)) * (s_H)*y + (1-delta)*h ) 
        hss = sm.solve(ss_h,h)[0]

        print('We now have these two values', 'k*=', kss , 'and h*=' , hss)

        print('We now need to substitute to find the real steady state values')

        # We will now do the substitution for h in kss and solve for k
        k_ss = sm.solve(kss.subs(h,hss),k)
        k_ss 

        # now we do the substitution for k i hss and solve for h
        h_ss = sm.solve(hss.subs(k,kss),h)
        h_ss

        print ('We now have the steady State values for h and k' , k_ss , h_ss)
        return
    from scipy.optimize import fsolve

    # ...

    def SteadyState2(self, do_print=False):
        par = self.sim

        # Define the steady-state equations for transformed variables
        ktilde_ss_eq = ((par.s_K**(1-par.phi) * par.s_H**par.phi) / (par.n + par.g + par.delta + par.n * par.g)) ** (1 / (1 - par.phi - par.alpha))
        htilde_ss_eq = ((par.s_K**par.alpha * par.s_H**(1-par.alpha)) / (par.n + par.g + par.delta + par.n * par.g)) ** (1 / (1 - par.phi - par.alpha))

        k_ss = ktilde_ss_eq
        h_ss = htilde_ss_eq

        if do_print:
            print('Steady state values:')
            print(f'k* = {k_ss:.4f}')
            print(f'h* = {h_ss:.4f}')

        return k_ss, h_ss



    def calculate_values(self, periods=10):
        sim = self.sim
        htilde_values = []
        ktilde_values = []

        for i in range(periods):
            htilde_val = i
            ktilde_val = i

            htilde_obj = ((sim.n + sim.g + sim.delta + sim.n * sim.g) / sim.s_K) ** (1 / sim.alpha) * htilde_val ** ((1 - sim.alpha) / sim.phi)
            ktilde_obj = ((sim.n + sim.g + sim.delta + sim.n * sim.g) / sim.s_H) ** (1 / sim.alpha) * ktilde_val ** ((1 - sim.phi) / sim.alpha)

            htilde_values.append(htilde_obj)
            ktilde_values.append(ktilde_obj)

        return htilde_values, ktilde_values

        
            
            

        




        






def solve_ss(alpha, c):
    """ Example function. Solve for steady state k. 

    Args:
        c (float): costs
        alpha (float): parameter

    Returns:
        result (RootResults): the solution represented as a RootResults object.

    """ 
    
    # a. Objective function, depends on k (endogenous) and c (exogenous).
    f = lambda k: k**alpha - c
    obj = lambda kss: kss - f(kss)

    #. b. call root finder to find kss.
    result = optimize.root_scalar(obj,bracket=[0.1,100],method='bisect')
    
    return result