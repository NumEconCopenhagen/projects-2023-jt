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
        Production = sm.Eq(par.Y_t, par.K_t**par.alpha * par.H_t**par.phi * (par.A_t * par.L_t)**(1 - par.alpha - par.phi))

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
        RentalRate = sm.Eq(par.r_t,par.alpha * (par.K_t/(par.A_t * par.L_t))**(par.alpha -1) * (par.H_t/(par.A_t * par.L_t))^par.phi)

        return

    def RealWageRate(self):
        par = self.par
        WageRate = sm.Eq(par.w_t, par.alpha * (par.K_t/(par.A_t * par.L_t))**par.alpha * (par.H_t/(par.A_t * par.L_t))**par.phi * par.A_t)

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
    

#class SimulationParameters:
    def __init__(self):
        # We create the namespaces:
        par = self.par = SimpleNamespace()
        sim = self.sim = SimpleNamespace()

        # Sim parameters
        sim.alpha = 0.333
        sim.phi = 0.333
        sim.delta = 0.02
        sim.n = 0.014
        sim.g = 0.016
        sim.s_K = 0.25
        sim.s_H = 0.129

        self.K0 = 1
        self.H0 = 1
        self.L0 = 1
        self.A0 = 1

    def simulate(self, periods):
        # Initialize the necessary variables
        K_vals = [self.K0]
        H_vals = [self.H0]
        L_vals = [self.L0]
        A_vals = [self.A0]
        Y_vals = []

        # Loop through the desired number of periods and calculate the variables
        for _ in range(periods):
            # Calculate output per capita (Y)
            Y = K_vals[-1]**self.sim.alpha * H_vals[-1]**self.sim.phi * (A_vals[-1] * L_vals[-1])**(1 - self.sim.alpha - self.sim.phi)
            Y_per_capita = Y / (L_vals[-1] * A_vals[-1])  # Output per capita

            Y_vals.append(Y_per_capita)

            # Calculate the updated values for capital per capita (K), human capital per capita (H), and technology (A)
            K_per_capita = K_vals[-1] / (L_vals[-1] * A_vals[-1])  # Capital per capita
            H_per_capita = H_vals[-1] / (L_vals[-1] * A_vals[-1])  # Human capital per capita
            A = (1 + self.sim.g) * A_vals[-1]  # Technology (A)

            # Calculate the updated values for labor (L)
            L = (1 + self.sim.n) * L_vals[-1]

            # Append the updated values to the respective lists
            K_vals.append(K_per_capita)
            H_vals.append(H_per_capita)
            L_vals.append(L)
            A_vals.append(A)

        # Plot the simulated values
        periods_range = range(periods + 1)
        plt.plot(periods_range, K_vals[:periods + 1], label="Physical Capital per Capita (~K)")  # Limit K_vals to match the specified number of periods
        plt.plot(periods_range, H_vals[:periods + 1], label="Human Capital per Capita (~H)")  # Limit H_vals to match the specified number of periods
        plt.plot(periods_range[:-1], Y_vals, label="Output per Capita (~Y)")  # Remove the last element from periods_range
        plt.xlabel("Periods")
        plt.ylabel("Value")
        plt.legend()
        plt.show()
        
class Simulation:
    def __init__(self):
        # Defining namespaces
        par = self.par = SimpleNamespace()
        sim = self.sim = SimpleNamespace()

        # Defining our sim parameters
        sim.alpha = 0.333
        sim.phi = 0.333
        sim.delta = 0.05
        sim.n = 0.014
        sim.g = 0.016
        sim.s_K = 0.25
        sim.s_H = 0.129

        # Defining our starting values
        sim.K0 = 1
        sim.H0 = 1
        sim.L0 = 1
        sim.A0 = 1
        
    # Defining our production function for our simulation
    
    def Productionfunction(self, K, H, A, L):
        sim = self.sim
        Yt = K**sim.alpha * H**sim.phi * (A * L)**(1 - sim.alpha - sim.phi)
        return Yt

    # Defining physical capital accumulation equation
    
    def Knextperiod(self, K, Yt):
        sim = self.sim
        Knext = sim.s_K * Yt - sim.delta * K + K
        return Knext

    # Defining human capital accumulation equation
    
    def Hnextperiod(self, H, Yt):
        sim = self.sim
        Hnext = sim.s_H * Yt - sim.delta * H + H
        return Hnext

        #Defining our population growth
        
    def Lnextperiod(self, L):
        sim = self.sim
        Lnext = (1 + sim.n) * L
        return Lnext
    
    # Defining our technology growth

    def Anextperiod(self, A):
        sim = self.sim
        Anext = (1 + sim.g) * A
        return Anext

    # Setting up our simulation
    def simulate(self, periods=100):
        
        # Creating empty arrays
        Yvalues = np.zeros(periods)
        Kvalues = np.zeros(periods)
        Hvalues = np.zeros(periods)
        Avalues = np.zeros(periods)
        Lvalues = np.zeros(periods)

        K = self.sim.K0
        H = self.sim.H0
        L = self.sim.L0
        A = self.sim.A0

        # Looping over our model equations
        for t in range(periods):
            Y = self.Productionfunction(K, H, A, L)
            K = self.Knextperiod(K, Y)
            H = self.Hnextperiod(H, Y)
            L = self.Lnextperiod(L)
            A = self.Anextperiod(A)

            # Updating our arrays with the calculated values
            Yvalues[t] = Y
            Kvalues[t] = K
            Hvalues[t] = H
            Avalues[t] = A
            Lvalues[t] = L

        # Calculating our tilde variables
        
        Y_per_capita = Yvalues / Lvalues
        K_per_capita = Kvalues / Lvalues
        H_per_capita = Hvalues / Lvalues

        Ytilde = Y_per_capita / Avalues
        Ktilde = K_per_capita / Avalues
        Htilde = H_per_capita / Avalues
        
        periods_range = range(periods)

        # Setting up the plot
        
        plt.plot(periods_range, Ytilde, label='Ytilde')
        plt.plot(periods_range, Ktilde, label='Ktilde')
        plt.plot(periods_range, Htilde, label='Htilde')
        plt.xlabel('Periods')
        plt.ylabel('Level')
        plt.title(f'Simulation of tilde variables for {periods} periods')
        plt.legend()
        plt.show()

        #return Kvalues, Hvalues, Yvalues, Avalues, Lvalues, Y_per_capita, K_per_capita, H_per_capita, Ytilde, Ktilde, Htilde

