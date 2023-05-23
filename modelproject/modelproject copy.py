from types import SimpleNamespace

from scipy import optimize
import numpy as np
import sympy as sm
from ipywidgets import interact
import pandas as pd 
import matplotlib.pyplot as plt
import ipywidgets as widgets
#from IPython.display import display



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
        par = self.sim
        Production = sm.Eq(par.Y_t,par.K_t**par.alpha * par.H_t**par.phi * (par.A_t * par.L_t)^(1-par.alpha - par.phi))

        return

    def LabourSupply(self):
        par = self.sim
        Labor_supply = sm.Eq(par.L_t1,(1+par.n)*par.L_t )

        return

    def LabourProductivity(self):
        par = self.sim
        Productivity = sm.Eq(par.A_t1 , (1+ par.g)*par.A_t)

        return

    def RentalrateOfCapital(self):
        par = self.sim
        RentalRate = sm.Eq(par.r_t,par.alpha * (par.K_t/(par.A_t * par.L_t))^(par.alpha -1) * (par.H_t/(par.A_t * par.L_t))^par.phi)

        return

    def RealWageRate(self):
        par = self.sim
        WageRate = sm.Eq(par.w_t, par.alpha * (par.K_t/(par.A_t * par.L_t))^par.alpha * (par.H_t/(par.A_t * par.L_t))^par.phi * par.A_t)

        return

    def HumanCapitalAccumulation(self):
        par = self.sim
        HumanCapital = sm.Eq(par.H_t1 , par.s_H * par.Y_t + par.delta * par.H_t)

        return

    def PhysicalCapitalAccumulation(self):
        par = self.sim
        PhysicalCapital = sm.Eq(par.K_t1 , par.s_K * par.Y_t + par.delta * par.K_t)
 
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
        
        return 

    def SteadyStateFunctions(alpha,phi,delta,n,g,s_K,s_H ,do_print=True):

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
                print(f"ktilde in steady state: {ktilde_steady_state}, htilde in steady state: {htilde_steady_state}")

                if do_plot:
                    plt.plot(ktilde_steady_state, htilde_steady_state, 'ro', label='Steady State')
                    
    def Nullclines2(self, do_sim=False, do_plot=True, periods=500, steady_state=True):
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
            if do_plot:
                plt.plot(ncht_vals, label="Δht=0")
                plt.plot(nckt_vals, label="Δkt=0")
                plt.xlim(0, periods-1)
                plt.ylim(0, periods-1)
                plt.xlabel('Level of physical capital')
                plt.ylabel('Level of human capital')
                plt.title('Phasediagram')

            # Calculate and display steady state if enabled
            if steady_state:
                try:
                    ktilde_expr = ((s_K**(1-phi) * s_H**phi)/(par.n + par.g + delta + par.n*par.g))**(1/(1-phi-alpha))
                    htilde_expr = ((s_K**alpha * s_H**(1-alpha))/(par.n + par.g + delta + par.n*par.g))**(1/(1-phi-alpha))
                    ktilde_func = sm.lambdify([], [ktilde_expr.subs({par.phi: phi, par.alpha: alpha, par.delta: delta, par.s_K: s_K, par.s_H: s_H, par.n: sim.n, par.g: sim.g})])
                    htilde_func = sm.lambdify([], [htilde_expr.subs({par.phi: phi, par.alpha: alpha, par.delta: delta, par.s_K: s_K, par.s_H: s_H, par.n: sim.n, par.g: sim.g})])
                    ktilde_steady_state = ktilde_func()[0]
                    htilde_steady_state = htilde_func()[0]
                    #print(f"ktilde in steady state: {ktilde_steady_state}, htilde in steady state: {htilde_steady_state}")
                    if do_plot:
                        plt.plot(ktilde_steady_state, htilde_steady_state, 'ro', label='Steady State')
                except Exception as e:
                    print(f"Error calculating steady state: {e}")

            # Display the legend
            if do_plot:
                plt.legend()

            # Show the plot
            if do_plot:
                plt.show()

        # Create FloatSliders for s_K, s_H, alpha_phi, delta, and periods
        s_K_slider = widgets.FloatSlider(min=0.01, max=1.0, step=0.01, value=0.25, description='s_K')
        s_H_slider = widgets.FloatSlider(min=0.01, max=1.0, step=0.01, value=0.129, description='s_H')
        alpha_phi_slider = widgets.FloatSlider(min=0.001, max=0.5, step=0.001, value=1/3, description='alpha/phi')
        delta_slider = widgets.FloatSlider(min=0.001, max=0.1, step=0.001, value=0.02, description='delta')
        periods_dropdown = widgets.Dropdown(options=list(range(100, 1001, 100)), value=100, description='periods')

        # Call the interactive function with the sliders as arguments
        widgets.interact(plot_function, s_K=s_K_slider, s_H=s_H_slider, alpha_phi=alpha_phi_slider, delta=delta_slider, periods=periods_dropdown)

   
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

    # Defining our population growth
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
    def simulate(self, periods=100, interactive=False, destroy_period=None):
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
            if destroy_period is not None and t == destroy_period:
                K *= 0.5  # Simulate half the capital being destroyed

            # Updating our arrays with the calculated values
            Yvalues[t] = Y
            Kvalues[t] = K
            Hvalues[t] = H
            Avalues[t] = A
            Lvalues[t] = L

        periods_range = range(periods)

        if interactive:
            # Create interactive sliders for s_K and s_H
            s_K_slider = widgets.FloatSlider(value=self.sim.s_K, min=0.1, max=0.5, step=0.05, description='s_K')
            s_H_slider = widgets.FloatSlider(value=self.sim.s_H, min=0.1, max=0.5, step=0.05, description='s_H')

            # Create update function for the sliders
            def update_simulation(s_K, s_H):
                self.sim.s_K = s_K
                self.sim.s_H = s_H
                self.simulate(periods, interactive=False)

            # Create interactive plot
            interact_plot = widgets.interact(update_simulation, s_K=s_K_slider, s_H=s_H_slider)

            # Display the interactive plot
            display(interact_plot)
        else:
            # Calculate the tilde variables without interactive plot
            Y_per_capita = Yvalues / Lvalues
            K_per_capita = Kvalues / Lvalues
            H_per_capita = Hvalues / Lvalues

            Ytilde = Y_per_capita / Avalues
            Ktilde = K_per_capita / Avalues
            Htilde = H_per_capita / Avalues
            
            # Plotting
            fig, ax = plt.subplots()
            ax.plot(periods_range, Ytilde, label='Ytilde')
            ax.plot(periods_range, Ktilde, label='Ktilde')
            ax.plot(periods_range, Htilde, label='Htilde')
            ax.set_xlabel('Periods')
            ax.set_ylabel('Level')
            ax.set_title(f'Simulation of tilde variables for {periods} periods')
            ax.legend()
            plt.show()



class ExtensionClass:
    def __init__(self):
        pass

    def steady_state(k,alpha,phi,delta,n,g,tau,s_K):
        k = sm.symbols('k')
        h = sm.symbols('h')
        alpha = sm.symbols('alpha')
        delta = sm.symbols('delta')
        tau = sm.symbols('tau')
        s_K = sm.symbols('s_K')
        g = sm.symbols('g')
        n = sm.symbols('n')
        phi = sm.symbols('phi')
        y = k**alpha * h**phi

        # We define the two stedy state-functions:
        # during steady state k_t+1 = k_t = k and h_t+1 = h_t = h
        ss_k = sm.Eq(k, 1/((1+n)*(1+g))* (s_K* k**alpha * h**phi + (1-delta-tau)*k) )
        ss_h = sm.Eq(h, 1/((1+n)*(1+g))* (s_H* k + (1-delta)*h) )
        # We solve the steady state functions for k and h, using sympy
        ss_k = sm.solve(ss_k,k)[0]
        ss_h = sm.solve(ss_h,h)[0]

        # We will now do substitutions, to get the steady state values for k and h
        kss = ss_k.subs(h,ss_h)
        hss = ss_h.subs(k,ss_k)

        return kss,hss 
    
    def ss_functions(alpha,phi,delta,n,g,tau,s_K) : 

        alpha = 0.33
        phi = 0.33
        delta = 0.02
        n = 0.014
        g = 0.016
        s_K = 0.25
        tau = 0.2

        # We define the two stedy state-functions:
        k_tilde = (s_K / (delta + tau + n + g + n*g)**(1/(1-alpha-phi))) * (tau/(delta + n + g + n*g))**((phi)/(1-alpha-phi))
        h_tilde = (s_K / (delta + tau + n + g + n*g)**(1/(1-alpha-phi))) * (tau/(delta + n + g + n*g))**((1-alpha)/(1-alpha-phi))

        # now we turn them in to pyhton function, using sympy lambdify
        kss_function = sm.lambdify((alpha,phi,delta,n,g,tau,s_K),k_tilde)
        hss_function = sm.lambdify((alpha,phi,delta,n,g,tau,s_K),h_tilde)

        # Now we call on the functions, to get the steady state values for k and h
        kss_function(alpha,phi,delta,n,g,tau,s_K)
        hss_function(alpha,phi,delta,n,g,tau,s_K)

        return 'The Steady State value for k is',kss_function,'and the Steady State value for h is',hss_function







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