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

        # We name our varibles:


        # We name our parameters:
        par.alpha = sm.symbols('alpha')
        par.phi  = sm.symbols('phi')
        par.delta = sm.symbols('delta')
        par.n = sm.symbols('n')
        par.g = sm.symbols('g')


        # We name our variables

        par.K_t = sm.symbols('K_{t}')
        par.K_t1 = sm.symbols('K_{t+1}')
        par.H_t = sm.symbols('H_{t}')
        par.H_t1 = sm.symbols('H_{t+1}')







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