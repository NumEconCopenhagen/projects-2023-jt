from types import SimpleNamespace

from scipy import optimize
import numpy as np
import sympy as sm

import pandas as pd 
import matplotlib.pyplot as plt
import ipywidgets 


class SolowModelClass:
    def parameters(self):
        
        # We create the namespaces: 
        par = self.par = SimpleNamespace()


        # We name our parameters:
        par.alpha = 0.33
        par.phi  = 0.5
        par.delta =  0.05
        par.n = 0.05
        par.g = 0.05
        par.s_K = 0.1
        par.s_H = 0.1


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

        P_C_Function = sm.lambdify((par.s_K,par.Y_t,par.delta,par.K_t),PhysicalCapital)



        return    

    def PhysicalCapitalFunction(self):
        par = self.par


    def SteadyStateValues(k,h,alpha,delta,s_K,s_H,g,n,phi, do_print=False):
        k = sm.symbols('k')
        h = sm.symbols('h')
        alpha = sm.symbols('alpha')
        delta = sm.symbols('delta')
        s_K = sm.symbols('s_K')
        s_H = sm.symbols('s_H')
        g = sm.symbols('g')
        n = sm.symbols('n')
        phi = sm.symbols('phi')

        y = k**alpha * h**phi

        # We define the function for which we are calculating the ss-value 
        ss_k = sm.Eq(k, 1/((1+n)*(1+g))*((s_K)*y+(1-delta)*k)) 
        # We find the steady state for k, by putting the lef hand side equal to 0
        kss = sm.solve(ss_k,k)[0]
                
        # We will now do the same for h
        ss_h = sm.Eq(h, 1/((1+n)*(1+g)) * ((s_H)*y+(1-delta)*h) ) 
        hss = sm.solve(ss_h,h)[0]

        ## print('We now have these two values', 'k*=', kss , 'and h*=' , hss)

        ## print('We now need to substitute to find the real steady state values')

        # We will now do the substitution for h in kss and solve for k
        k_ss = kss.subs(h,hss)

        # now we do the substitution for k i hss and solve for h
        h_ss = hss.subs(k,kss)

        print('k_ss = ' , sm.latex(k_ss) ,'h_ss = ' , sm.latex(h_ss))
        ##print ('We now have the steady State values for h and k' , k_ss , h_ss)
    
        return 

    def SteadyStateFunctions(k,h,alpha,delta,s_K,s_H,g,n,phi):

        kss_function = sm.lambdify((alpha,delta,s_K,s_H,g,n,phi),k_ss)
        hss_function = sm.lambdify((alpha,delta,s_K,s_H,g,n,phi),h_ss) 

        return kss_function and hss_function 
   

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