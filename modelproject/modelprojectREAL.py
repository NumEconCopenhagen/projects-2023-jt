
# Importing the packages we need
from types import SimpleNamespace
from scipy import optimize
import numpy as np
import sympy as sm
from ipywidgets import interact
import pandas as pd 
import matplotlib.pyplot as plt
import ipywidgets as widgets
from IPython.display import display

# Creating a class for the general Solow model

class SolowGeneral:
    def __init__(self):

        # Setting up parameters
        # Defining namespaces
        par = self.par = SimpleNamespace()
        sim = self.sim = SimpleNamespace()

        # We start by setting up our algebraic parameters
        par.K = sm.symbols('k')
        par.A = sm.symbols('a')
        par.L = sm.symbols('L')
        par.Y = sm.symbols('Y')
        par.S = sm.symbols('S')
        
        par.alpha = sm.symbols('alpha')
        par.delta = sm.symbols('delta')
        par.s = sm.symbols('s')
        par.g = sm.symbols('g')
        par.n = sm.symbols('n')

        # We defien our capital per effective worker
        par.k = sm.symbols('k')

        # Defining our sim variables
        sim.alpha = 1/3
        sim.delta = 0.02
        sim.s = 0.2  
        sim.g = 0.02
        sim.n = 0.01
        self.num_periods = 100
        self.k_0 = 1e-7


    
    def steadystate_analytical(self):
        par = self.par
        # We define our transition equation
        f = par.k**par.alpha
        tranisition = sm.Eq(par.k, 1/((1+par.n)*(1+par.g))*((1-par.delta)*par.k + par.s*f))

        # We solve for the steady state
        steady_state = sm.solve(tranisition, par.k)[0]

        return steady_state
    
    def steadystate_value(self): 

        par = self.par
        sim = self.sim

        # We turn our symbolic steady state into a function
        steady_state_func = sm.lambdify((par.s,par.g,par.n,par.delta,par.alpha), self.steadystate_analytical())

        return steady_state_func(sim.s,sim.g,sim.n,sim.delta,sim.alpha)
    
    def steadystate_numerical(self) :
        par = self.par
        sim = self.sim

        # We define the equation we want to solve
        f = lambda k: par.k**par.alpha
        tranisition = lambda k: k - 1/((1+par.n)*(1+par.g))*((1-par.delta)*k + par.s*f(k))

        ss = optimize.root_scalar(tranisition, bracket=[0.1, 100], method='brentq')

        return ss.root
    
    def plot_transition_diagram(self):
        """
        Plot the transition diagram for the capital stock.
        """
        k_values = [self.k_0]  # List to store the capital stock values

        # Calculate the capital stock values for each period
        for _ in range(self.num_periods):
            k_t = k_values[-1]  # Current capital stock
            k_t1 = 1 / ((1 + self.sim.n) * (1 + self.sim.g)) * (
                    self.sim.s * k_t ** self.sim.alpha + (1 - self.sim.delta) * k_t)
            k_values.append(k_t1)

        # Transition diagram plotting
        plt.figure(figsize=(10, 6)) 
        plt.plot(k_values[:-1], k_values[1:], 'bo-', markersize=0, linewidth=1, label='k_t1')
        plt.plot(k_values[:-1], k_values[:-1], 'ro-', markersize=0, linewidth=1, label='k_t')
        plt.xlabel('Capital Stock at t')
        plt.ylabel('Capital Stock at t+1')
        plt.title('Transition Diagram: Capital Stock')
        plt.legend()
        plt.grid(True)
        plt.show()

    def display_interactive_plot(self):
        # FloatSliders for simulation variables
        num_periods_slider = widgets.IntSlider(min=1, max=500, step=1, value=300, description='Num Periods')
        alpha_slider = widgets.FloatSlider(min=0, max=0.5, step=0.01, value=1/3, description='Alpha')
        delta_slider = widgets.FloatSlider(min=0, max=0.5, step=0.01, value=0.02, description='Delta')
        n_slider = widgets.FloatSlider(min=0, max=0.5, step=0.01, value=0.01, description='n')
        g_slider = widgets.FloatSlider(min=0, max=0.5, step=0.01, value=0.02, description='g')
        s_slider = widgets.FloatSlider(min=0, max=0.5, step=0.01, value=0.20, description='s')
        num_periods_slider = widgets.IntSlider(min=0, max=500, step=1, value=200, description='Num Periods')

        # Interactive update of simulation variables
        def update_sliders(num_periods, alpha, delta, n, g, s):
            self.num_periods = num_periods
            self.sim.alpha = alpha
            self.sim.delta = delta
            self.sim.n = n
            self.sim.g = g
            self.sim.s = s
            self.plot_transition_diagram()

        interactive_plot = widgets.interactive(update_sliders,
                                               alpha=alpha_slider, delta=delta_slider, n=n_slider, g=g_slider, s=s_slider, num_periods=num_periods_slider)
        display(interactive_plot)

        
    def simulation(self, periods=1000):
        t_values = range(periods)

        sim = self.sim

        Y = np.zeros(periods)
        K = np.zeros(periods)
        L = np.zeros(periods)
        A = np.zeros(periods)

        alpha_slider = widgets.FloatSlider(min=0, max=0.5, step=0.01, value=1/3, description='Alpha:')
        delta_slider = widgets.FloatSlider(min=0, max=0.5, step=0.01, value=0.02, description='Delta:')
        s_slider = widgets.FloatSlider(min=0, max=1, step=0.01, value=0.2, description='S:')
        n_slider = widgets.FloatSlider(min=0, max=0.5, step=0.01, value=0.01, description='N:')
        g_slider = widgets.FloatSlider(min=0, max=0.5, step=0.01, value=0.02, description='G:')
        periods_slider = widgets.IntSlider(value=periods, min=100, max=1000, step=100, description='Periods:')

        @widgets.interact(alpha=alpha_slider, delta=delta_slider, s=s_slider, n=n_slider, g=g_slider, periods=periods_slider)
        def update_simulation(alpha, delta, s, n, g, periods):
            nonlocal Y, K, L, A

            sim.alpha = alpha
            sim.delta = delta
            sim.s = s
            sim.n = n
            sim.g = g

            K_t = 1
            L_t = 1
            A_t = 1
            S_t = 1

            for t in t_values:
                Y[t] = K_t ** sim.alpha * (A_t * L_t) ** (1 - sim.alpha)
                K_tnext = S_t + (1 - sim.delta) * K_t
                S_t = sim.s * Y[t]
                L_tnext = (1 + sim.n) * L_t
                A_tnext = (1 + sim.g) * A_t

                K[t] = K_t
                L[t] = L_t
                A[t] = A_t

                # Update the variables for the next time step
                K_t = K_tnext
                L_t = L_tnext
                A_t = A_tnext

            kpercap = K / L
            ktilde = kpercap / A

            plt.plot(ktilde, label="Technology adjusted capital per worker")
            plt.ylabel("Level of capital")
            plt.xlabel("Periods")
            plt.grid(True)
            plt.legend(loc=4)
            plt.title(f"Simulation of capital for {periods} periods")
            plt.xlim(0, periods)
        
        


# Creating a class for the extended Solowmodel with human capital

class SolowModelClass:
    def __init__(self):
                # We create the namespaces: 
        par = self.par = SimpleNamespace()
        sim = self.sim = SimpleNamespace()

        # We define our sim parameters, these will create the ground for our simulation.
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
        par.H_t = sm.symbols('H_{t}')
        par.Y_t = sm.symbols('Y_{t}')
        par.A_t = sm.symbols('A_{t}')
        par.L_t = sm.symbols('L_{t}')
        
        # # We name our per effective worker variables
        par.ktilde_t = sm.symbols('\tilde{k_{t}}')
        par.htilde_t = sm.symbols('\tilde{h_{t}}')
        par.ytilde_t = sm.symbols('\tilde{y_{t}}')


        
    def SteadyStateValues_k(k,h,alpha,delta,s_K,s_H,g,n,phi, do_print=False):
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
        return k_ss
    
    def SteadyStateValues_h(k,h,alpha,delta,s_K,s_H,g,n,phi, do_print=False):
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

        return h_ss
 

    def SteadyStateFunctions(self,alpha,phi,delta,n,g,s_K,s_H):

        par = self.sim
        alpha = par.alpha
        phi = par.phi
        delta = par.delta
        n = 0.014
        g= 0.016
        s_K = 0.25
        s_H = 0.129

        # We define the steady state functions
        k_tilde = ((s_K**(1-phi) * s_H**phi)/(n+g+delta +n*g))**(1/(1-phi-alpha))
        h_tilde = ( (s_K**(alpha) * s_H**(1-alpha))/(n+g+delta +n*g))**(1/(1-phi-alpha))
        
        # Now we turn them in to pyhton function, using sympy lambdify.
        kss_function = sm.lambdify((alpha,phi,delta,n,g,s_K,s_H),k_tilde)
        hss_function = sm.lambdify((alpha,phi,delta,n,g,s_K,s_H),h_tilde) 

        #Now we calculate the steady states
        kss_function(alpha,phi,delta,n,g,s_K,s_H)
        hss_function(alpha,phi,delta,n,g,s_K,s_H)

        return 'The steady state for k is ', kss_function(alpha,phi,delta,n,g,s_K,s_H) ,'and the steady state for h is',  hss_function(alpha,phi,delta,n,g,s_K,s_H)


    # We now define our method for simulating the Nullclines for our extended model
    def Nullclines(self, periods=500):
        par = self.par
        sim = self.sim
        periods = periods

        # Define the interactive function
        def plot_function(s_K, s_H, alpha_phi, delta, periods):
            alpha = alpha_phi
            phi = alpha_phi

            # Create the lambdified functions with updated s_K, s_H, alpha, phi, and delta values
            ncht_expr = (((par.n + par.g + delta + par.n * par.g) / s_K) ** (1 / phi)) * (par.ktilde_t ** ((1 - alpha) / phi))
            nckt_expr = (s_H/(par.n+par.g+delta+par.n*par.g))**(1/(1-phi))*par.ktilde_t**(alpha/(1-alpha))
            ncht_func = sm.lambdify(par.ktilde_t, ncht_expr.subs({par.alpha: alpha, par.phi: phi, par.delta: delta, par.n: sim.n, par.g: sim.g, par.s_K: s_K}))
            nckt_func = sm.lambdify(par.ktilde_t, nckt_expr.subs({par.alpha: alpha, par.phi: phi, par.delta: delta, par.n: sim.n, par.g: sim.g, par.s_H: s_H}))

            # Evaluate the functions for different t_values
            ktilde_vals = np.linspace(0, periods-1, periods)
            ncht_vals = ncht_func(ktilde_vals)
            htilde_vals = np.linspace(0, periods-1, periods)
            nckt_vals = nckt_func(htilde_vals)

            # Create the plot
                       
            plt.plot(ncht_vals, label="Δht=0")
            plt.plot(nckt_vals, label="Δkt=0")
            plt.xlim(0, periods-1)
            plt.ylim(0, periods-1)
            plt.xlabel('Level of physical capital')
            plt.ylabel('Level of human capital')
            plt.title('Phasediagram')

            # Calculate and display steady state
            
            try:
                ktilde_expr = ((s_K**(1-phi) * s_H**phi)/(par.n + par.g + delta + par.n*par.g))**(1/(1-phi-alpha))
                htilde_expr = ((s_K**alpha * s_H**(1-alpha))/(par.n + par.g + delta + par.n*par.g))**(1/(1-phi-alpha))
                ktilde_func = sm.lambdify([], [ktilde_expr.subs({par.phi: phi, par.alpha: alpha, par.delta: delta, par.s_K: s_K, par.s_H: s_H, par.n: sim.n, par.g: sim.g})])
                htilde_func = sm.lambdify([], [htilde_expr.subs({par.phi: phi, par.alpha: alpha, par.delta: delta, par.s_K: s_K, par.s_H: s_H, par.n: sim.n, par.g: sim.g})])
                ktilde_steady_state = ktilde_func()[0]
                htilde_steady_state = htilde_func()[0]
                plt.plot(ktilde_steady_state, htilde_steady_state, 'ro', label='Steady State')
            except Exception as e:
                print(f"Error calculating steady state: {e}")

            # Display the legend
            
            plt.legend()

            # Show the plot
            
            plt.show()

        # Create FloatSliders for s_K, s_H, alpha_phi, delta, and periods
        s_K_slider = widgets.FloatSlider(min=0.01, max=1.0, step=0.01, value=0.25, description='s_K')
        s_H_slider = widgets.FloatSlider(min=0.01, max=1.0, step=0.01, value=0.129, description='s_H')
        alpha_phi_slider = widgets.FloatSlider(min=0.001, max=0.5, step=0.001, value=1/3, description='alpha/phi')
        delta_slider = widgets.FloatSlider(min=0.001, max=0.1, step=0.001, value=0.02, description='delta')
        periods_dropdown = widgets.Dropdown(options=list(range(100, 1001, 100)), value=100, description='periods')

        # Call the interactive function with the sliders as arguments
        widgets.interact(plot_function, s_K=s_K_slider, s_H=s_H_slider, alpha_phi=alpha_phi_slider, delta=delta_slider, periods=periods_dropdown)




# Creating a class for the simulations that we need to do
class SimulationClass:
    def __init__(self):
        # Defining namespaces
        par = self.par = SimpleNamespace()
        sim = self.sim = SimpleNamespace()

        # Defining our sim parameters
        sim.alpha = 1/3
        sim.phi = 1/3
        sim.delta = 0.02
        sim.n = 0.014
        sim.g = 0.016
        sim.s_K = 0.25
        sim.s_H = 0.129
        sim.tau = 0  # Removed the enable_tax parameter and set tau to 0 by default

        # Defining our starting values
        sim.K0 = 1
        sim.H0 = 1
        sim.L0 = 1
        sim.A0 = 1

    def Productionfunction(self, K, H, A, L):
        sim = self.sim
        Yt = K**sim.alpha * H**sim.phi * (A * L)**(1 - sim.alpha - sim.phi)
        return Yt

    def Knextperiod(self, K, Yt):
        sim = self.sim
        Knext = sim.s_K * Yt - sim.delta * K + K
        return Knext

    def Hnextperiod(self, H, Yt, K):
        sim = self.sim
        Hnext = sim.s_H * Yt - sim.delta * H + H
        return Hnext

    def Lnextperiod(self, L):
        sim = self.sim
        Lnext = (1 + sim.n) * L
        return Lnext

    def Anextperiod(self, A):
        sim = self.sim
        Anext = (1 + sim.g) * A
        return Anext

    def simulate(self, periods=100, interactive=True):
        Yvalues = np.zeros(periods)
        Kvalues = np.zeros(periods)
        Hvalues = np.zeros(periods)
        Avalues = np.zeros(periods)
        Lvalues = np.zeros(periods)

        K = self.sim.K0
        H = self.sim.H0
        L = self.sim.L0
        A = self.sim.A0

        for t in range(periods):
            Y = self.Productionfunction(K, H, A, L)
            K = self.Knextperiod(K, Y)
            H = self.Hnextperiod(H, Y, K)
            L = self.Lnextperiod(L)
            A = self.Anextperiod(A)

            Yvalues[t] = Y
            Kvalues[t] = K
            Hvalues[t] = H
            Avalues[t] = A
            Lvalues[t] = L

        periods_range = range(periods)

        if interactive:
            self._create_s_H_plot(periods)
        else:
            Y_per_capita = Yvalues / Lvalues
            K_per_capita = Kvalues / Lvalues
            H_per_capita = Hvalues / Lvalues

            Ytilde = Y_per_capita / Avalues
            Ktilde = K_per_capita / Avalues
            Htilde = H_per_capita / Avalues

            fig, ax = plt.subplots()
            ax.plot(periods_range, Ytilde, label='Ytilde')
            ax.plot(periods_range, Ktilde, label='Ktilde')
            ax.plot(periods_range, Htilde, label='Htilde')
            ax.set_xlabel('Periods')
            ax.set_ylabel('Level')
            ax.set_title(f'Simulation of tilde variables for {periods} periods')
            ax.legend()
            plt.show()

    def _create_s_H_plot(self, periods):
        s_K_slider = widgets.FloatSlider(value=self.sim.s_K, min=0.01, max=0.5, step=0.01, description='s_K')
        s_H_slider = widgets.FloatSlider(value=self.sim.s_H, min=0.01, max=0.5, step=0.01, description='s_H')

        def update_simulation(s_K, s_H):
            self.sim.s_K = s_K
            self.sim.s_H = s_H
            self.simulate(periods, interactive=False)

        interact_plot = widgets.interact(update_simulation, s_K=s_K_slider, s_H=s_H_slider)
        display(interact_plot)