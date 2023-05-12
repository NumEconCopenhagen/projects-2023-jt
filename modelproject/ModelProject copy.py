from types import SimpleNamespace
from scipy import optimize
import numpy as np
import sympy as sm
import pandas as pd 
import matplotlib.pyplot as plt
import ipywidgets 


class SolowModelClass:
    def __init__(self):
   
        par = self.par = SimpleNamespace()

        #We start by setting our parameters (these are assumed constant)
        par.alpha = 1/3         # Share of capital
        par.phi = 1/3           # 
        par.delta = 0.02        # Depreciation rate
        par.n = 0.014           # Population growth rate
        par.g = 0.016           # Productivity
        par.s_K = 0.25          # Savings rate for physical kapital
        par.s_H = 0.129         # Savings rate for human capital
        par.A_t = 1             # The technological devolopment
        
        # We now define the variables
        par.L_t = 2             #current labour Force

        par.ktilde = 1
        par.htilde = 1

        par.alpha , par.phi , par.delta , par.n , par.g , par.s_K , par.s_H = sm.symbols('alpha') , sm.symbols('phi') , sm.symbols('delta') ,sm.symbols('n') , sm.symbols('g') ,sm.symbols('s_{K}') ,sm.symbols('s_{H}')   

        par.K_t = par.ktilde * (par.A_t * par.L_t)
        par.H_t = par.htilde * (par.A_t * par.L_t)

        # We define the production function
        par.Y = par.K_t**par.alpha * par.H_t**par.phi * (par.A_t * par.L_t)**(1-par.alpha - par.phi)
        # We define the per effective worker production function
        par.ytilde = par.Y /(par.A_t * par.L_t)

    "We now define all the function for the next period values, and then we update the current state"

    def next_period_physicalcapital(self):
        " We are calculating it as per effective worker physical capital."
        par = self.par
        # We get get rid of the par. notation so simplify the notation
        n , g, s_K,y_tilde,delta,k_tilde =  par.n,par.g,par.s_K,par.ytilde,par.delta,par.ktilde

        # We define the function
        k_next = (s_K * y_tilde + (1-delta)*k_tilde) / ((1+n)*(1+g))

        # We turn it in to af python function
        k_tildenext = sm.lambdify((n , g, s_K,y_tilde,delta,k_tilde),k_next)
        
        return k_tildenext(n , g, s_K,y_tilde,delta,k_tilde)

    def next_period_humancapital(self):
        " We are calculating it as per effective worker physical capital."
        par = self.par
        # We get get rid of the par. notation so simplify the notation
        n , g, s_H,y_tilde,delta,h_tilde = par.n,par.g,par.s_H,par.ytilde,par.delta,par.htilde

        # We define the function
        h_next = (s_H * y_tilde + (1-delta)*h_tilde) / ((1+n)*(1+g))

        # We turn it in to af python function
        h_tildenext = sm.lambdify((n , g, s_H,y_tilde,delta,h_tilde),h_next)

        return h_tildenext(n , g, s_H,y_tilde,delta,h_tilde)
    
    def next_period_labourforce(self):
        "We are calculating the next period labour force"
        par = self.par
        # We get get rid of the par. notation so simplify the notation
        L_t , n = par.L_t , par.n 

        # We define the function
        Lt_next = (1+n)*L_t 

        # We turn it in to a python function
        L_next = sm.lambdify((n,L_t),Lt_next)

        return L_next(n,L_t)

    def update(self):
        "We are updating the current states"
        par = self.par

        par.ktilde = self.next_period_physicalcapital()
        par.htilde = self.next_period_humancapital()
        par.L_t = self.next_period_labourforce()

        par.K_t = par.ktilde * (par.A_t * par.L_t)
        par.H_t = par.htilde * (par.A_t * par.L_t)
        par.Y = par.K_t**par.alpha * par.H_t**par.phi * (par.A_t * par.L_t)^(1-par.alpha - par.phi)
        # We compute the per effective worker production function
        par.ytilde = par.Y /(par.A_t * par.L_t)

    
    def SteadyStateValues(k,h,alpha,delta,s_K,s_H,g,n,phi, do_print=False):
        "This function is only used to compute the steady states expressions"
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

        # We will now do the substitution for h in kss and solve for k
        k_ss = kss.subs(h,hss)

        # now we do the substitution for k i hss and solve for h
        h_ss = hss.subs(k,kss)

        print('k_ss = ' , sm.latex(k_ss) ,'h_ss = ' , sm.latex(h_ss))
        ##print ('We now have the steady State values for h and k' , k_ss , h_ss)
    
        return 

    "We make functions for the two steady state functions, to make them easier to evaluate"
    def SteadyStatek_function(self):
        par = self.par

        # We get rid of the par. notation for easier readability
        alpha,delta,phi,n,g,s_K,s_H = par.alpha ,par.delta,par.phi,par.n,par.g,par.s_K,par.s_H

        # We define the seady state function
        k_tilde = ((s_K**(1-phi) * s_H**phi)/(n+g+delta +n*g))**(1/(1-phi-alpha))

        # We turn it in to a python function
        kss_function = sm.lambdify((alpha,phi,delta,n,g,s_K,s_H),k_tilde)

        return  kss_function(alpha,phi,delta,n,g,s_K,s_H)

    def SteadyStateh_funcion(self): 
        par = self.par

        # We get rid of the par. notation for easier readability
        alpha,delta,phi,n,g,s_K,s_H = par.alpha ,par.delta,par.phi,par.n,par.g,par.s_K,par.s_H

        # We define the steady state function
        h_tilde = ( (s_K**(alpha) * s_H**(1-alpha))/(n+g+delta +n*g))**(1/(1-phi-alpha))

        # We turn it in to a python function
        hss_function = sm.lambdify((alpha,phi,delta,n,g,s_K,s_H),h_tilde) 

        return hss_function(alpha,phi,delta,n,g,s_K,s_H)

   


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