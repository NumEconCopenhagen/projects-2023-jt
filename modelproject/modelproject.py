from types import SimpleNamespace
from scipy import optimize
import numpy as np
import sympy as sm


import pandas as pd 
import matplotlib.pyplot as plt
import ipywidgets 


class SolowModel:
    def solow(self):
        
        # We create the namespaces: 
        par = self.par = SimpleNamespace()
        sol = self.sol = SimpleNamespace()


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
        

        # We now define our equations as given by the book

        par.Y_t = par.K_t**par.alpha * par.H_t**par.phi * (par.A_t * parL_t)^(1-par.alpha - par.phi)
        par.r_t = par.alpha * (par.K_t/(par.A_t * par.L_t))^(par.alpha -1) * (par.H_t/(par.A_t * par.L_t))^par.phi
        par.w_t = par.alpha * (par.K_t/(par.A_t * par.L_t))^par.alpha * (par.H_t/(par.A_t * par.L_t))^par.phi * par.A_t
        par.K_t1 = par.s_K * par.Y_t - par.delta * par.K_t + par.K_t
        par.H_t1 = par.s_H * par.Y_t - par.delta * par.H_t + par.H_t
        par.L_t1 = (1+par.n)*par.L_t
        par.A_t1 = (1+ par.g)*par.A_t

    
        









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