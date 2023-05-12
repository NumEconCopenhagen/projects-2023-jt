from scipy import optimize
import numpy as np
import sympy as sm
from scipy.optimize import fsolve
from types import SimpleNamespace
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



    def Nullclines(self, do_sim=False, do_plot=True, periods=100, steady_state=True):
        par = self.par
        sim = self.sim
        periods = periods

        if do_sim:
            # Create the lambdified functions
            ncht_expr = ((par.n + par.g + par.delta + par.n * par.g) / par.s_K) ** (1 / par.phi) * par.ktilde_t ** ((1 - par.alpha) / par.phi)
            nckt_expr = (par.s_H/(par.n+par.g+par.delta+par.n*par.g))**(1/(1-par.phi))*par.ktilde_t**(par.alpha/(1-par.alpha))
            ncht_func = sm.lambdify(par.ktilde_t, ncht_expr.subs({par.alpha: sim.alpha, par.phi: sim.phi, par.delta: sim.delta, par.n: sim.n, par.g: sim.g, par.s_K: sim.s_K}))
            nckt_func = sm.lambdify(par.ktilde_t, nckt_expr.subs({par.alpha: sim.alpha, par.phi: sim.phi, par.delta: sim.delta, par.n: sim.n, par.g: sim.g, par.s_H: sim.s_H}))

            # Evaluate the functions for different t_values
            ktilde_vals = np.linspace(0, periods-1, periods)
            ncht_vals = ncht_func(ktilde_vals)
            htilde_vals = np.linspace(0, periods-1, periods)
            nckt_vals = nckt_func(htilde_vals)

            if do_plot:
                plt.plot(ncht_vals, label="ht=0")
                plt.plot(nckt_vals, label="kt=0")
                plt.xlim(0, periods-1)
                plt.ylim(0, periods-1)
                plt.xlabel('ktilde')
                plt.ylabel('htilde')
                plt.title('Nullclines')

            if steady_state:
                ktilde_expr = ((par.s_K**(1-par.phi) * par.s_H**par.phi)/(par.n + par.g + par.delta + par.n*par.g))**(1/(1-par.phi-par.alpha))
                htilde_expr = ((par.s_K**par.alpha * par.s_H**(1-par.alpha))/(par.n + par.g + par.delta + par.n*par.g))**(1/(1-par.phi-par.alpha))
                ktilde_func = sm.lambdify([], [ktilde_expr.subs({par.phi: sim.phi, par.alpha: sim.alpha, par.s_K: sim.s_K, par.s_H: sim.s_H, par.n: sim.n, par.g: sim.g, par.delta: sim.delta})])
                htilde_func = sm.lambdify([], [htilde_expr.subs({par.phi: sim.phi, par.alpha: sim.alpha, par.s_K: sim.s_K, par.s_H: sim.s_H, par.n: sim.n, par.g: sim.g, par.delta: sim.delta})])
                ktilde_steady_state = ktilde_func()[0]
                htilde_steady_state = htilde_func()[0]
                print(ktilde_steady_state,htilde_steady_state)

                if do_plot:
                    plt.plot(ktilde_steady_state, htilde_steady_state, 'ro', label='Steady State')
    

