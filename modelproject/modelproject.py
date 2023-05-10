from types import SimpleNamespace
from scipy import optimize
import numpy as np
import sympy as sm


import pandas as pd 
import matplotlib.pyplot as plt
import ipywidgets 


class SolowModel
    def solow(self)
        
        # We create the namespaces: 
        par = self.par = SimpleNamespace()

        # We name our varibles:


        # We name our parameters:
        par.alpha = sm.symbols('alpha')
        par.phi  = sm.symbols('phi')





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