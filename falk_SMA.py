"""
--------------------------------------------------------------------------------------------------------------
This script implements the one-dimensional shape memory alloy (SMA) model proposed by Falk [1, 2]. 

The model is based on a Landau-type free energy function that describes the phase transformation behavior of
SMA materials.

The free energy is defined as a function of strain and temperature, such that:
    1. High temperatures (T > Ta): The energy has a single minimum at zero strain associated with stable 
                                   austenite (A) phase.
    2. Intermedia temperatures (Tm < T < Ta): The energy exhibits three distinct minima. One at zero strain 
                                              corresponding to the stable austenite (A), and two symmetric 
                                              minima at non-zero strains representing detwinned martensite 
                                              induced by compression (M-) and tension (M+), respectively.
    3. Low temperatures (T < Tm): The energy has two symmetric minima at non-zero strains associated with
                                  stable martensite phases (M- and M+), while the austenite phase becomes 
                                  unstable.

The parameter "alpha" from reference [1] is modified to "b^2/(4*a*(Ta - Tm))", as done by in Savi & Braga [2]. 
The modification is performed to define the temperature "Ta" above which the austenite phase is stable, and 
the free energy has only one minimum at zero strain.
--------------------------------------------------------------------------------------------------------------    
References:

[1]: Falk, F. (1980) - Model free energy, mechanics and thermodynamics of shape memory alloys. Acta 
     Metallurgica, 28, 1773-1780. (https://doi.org/10.1016/0001-6160(80)90030-9)
[2]: Falk, F. (1983) - One-dimensional model of shape memory alloys. Archives of Mechanics, 35(1), 63-84.
     (https://rcin.org.pl/ippt/publication/161795)
[3]: Savi & Braga (1993) - Chaotic vibrations of an oscillator with shape memory. Journal of the Brazilian
     Society of Mechanical Sciences, 15(1), 1-20. 
     (https://www.researchgate.net/publication/233783381_Chaotic_Vibrations_of_an_Oscillator_with_Shape_Memory)
--------------------------------------------------------------------------------------------------------------
Author: Luã G Costa [https://github.com/lguedesc]
Created: 02 Sep 2025
Last Update: 04 Sep 2025
--------------------------------------------------------------------------------------------------------------
"""

import matplotlib.pyplot as plt
import matplotlib.figure as mfig
import numpy as np
from typing import Callable

def newton_raphson(f: Callable[..., float], df: Callable[..., float], x0: float, 
                   tol: float = 1e-6, max_iter: int = 100,
                   *args, **kwargs) -> float:
    """
    ----------------------------------------------------------------------------
    This function implements the Newton-Raphson method for finding roots of a 
    function.
    ----------------------------------------------------------------------------
    Description of function arguments:
    - f: Callable[..., float]
        The function for which the roots are to be found.
    - df: Callable[..., float]
        The derivative of the function f.
    - x0: float
        Initial guess for the root.
    - tol: float
        Tolerance for the error. Default is 1e-6.
    - max_iter: int
        Maximum number of iterations. Default is 100.
    - *args, **kwargs
        Additional arguments to be passed to the function, f, and its 
        derivative, df.
    ----------------------------------------------------------------------------
    Returns: float
        The estimated root of the function.
    ----------------------------------------------------------------------------
    """    
    x = x0
    for i in range(max_iter):
        function = f(x, *args, **kwargs)
        derivative = df(x, *args, **kwargs) 

        if abs(derivative) < 1e-12:
            raise ValueError("Derivative is too small; no convergence.")

        x_new = x - (function/derivative)
        
        if abs(x_new - x) < tol:
            return x_new
        
        x = x_new
    
    raise ValueError("Maximum iterations reached; no convergence.")

class FalkSMAModel:
    """
    ----------------------------------------------------------------------------
    This class implements the Falk model for shape memory alloys.
    ----------------------------------------------------------------------------
    Parameters:
    - a: float
        First material parameter [MPa/ºC] or [MPa/K].
    - b: float
        Second material parameter [MPa/ºC] or [K].
    - Ta: float
        Austenite start temperature [ºC] or [K].
    - Tm: float
        Martensite start temperature [ºC] or [K].
    - T: float
        Current temperature [ºC] or [K].
    ----------------------------------------------------------------------------
    """ 
    def __init__(self, a: float, b: float, Ta: float, Tm: float, T: float):
        self.a = a      #[MPa/ºC] 
        self.b = b      #[MPa/ºC]
        self.Ta = Ta    #[ºC]
        self.Tm = Tm    #[ºC]
        self.T = T      #[ºC]
        
    def free_energy(self, epsilon: np.ndarray[np.float64]) -> np.ndarray[np.float64]:
        """
        ------------------------------------------------------------------------
        Obtains the free energy array for a given strain array.
        ------------------------------------------------------------------------
        Parameters:
        - epsilon: np.ndarray[np.float64]
            Strain array.  
        ------------------------------------------------------------------------
        Returns: np.ndarray[np.float64]
            The free energy array corresponding to the given strain array.
        ------------------------------------------------------------------------
        """
        psi = (self.a/2.0)*(self.T - self.Tm)*(epsilon**2) - (self.b/4.0)*(epsilon**4) + ((self.b**2)/(24*self.a*(self.Ta - self.Tm)))*(epsilon**6)
        return psi
    
    def stress(self, epsilon: np.ndarray[np.float64]) -> np.ndarray[np.float64]:
        """
        ------------------------------------------------------------------------
        Obtains the stress array for a given strain array.
        ------------------------------------------------------------------------
        Parameters:
        - epsilon: np.ndarray[np.float64]
            Strain array.  
        ------------------------------------------------------------------------
        Returns: np.ndarray
            The stress array corresponding to the given strain array.
        ------------------------------------------------------------------------
        """
        sigma = self.a*(self.T - self.Tm)*epsilon - self.b*(epsilon**3) + ((self.b**2)/(4*self.a*(self.Ta - self.Tm)))*(epsilon**5)
        return sigma
    
    def stress_for_root(self, epsilon: float) -> float:
        """
        ------------------------------------------------------------------------
        Auxiliary function to find the root of the stress equation for a given
        target stress using the Newton-Raphson method.
        ------------------------------------------------------------------------
        Parameters:
        - epsilon: np.ndarray[np.float64]
            Strain array.  
        ------------------------------------------------------------------------
        Returns: float
            The value of the stress equation minus the target stress.
        ------------------------------------------------------------------------
        """
        value = self.a*(self.T - self.Tm)*epsilon - self.b*(epsilon**3) + ((self.b**2)/(4*self.a*(self.Ta - self.Tm)))*(epsilon**5) - self.sigma_target
        return value
    
    def stress_derivative(self, epsilon: float) -> np.ndarray[np.float64]:
        """
        ------------------------------------------------------------------------
        Obtain the time derivative of the stress array for a given strain array.
        ------------------------------------------------------------------------
        Parameters:
        - epsilon: np.ndarray[np.float64]
            Strain array.  
        ------------------------------------------------------------------------
        Returns: np.ndarray[np.float64]
            The time derivative of the stress array corresponding to the given 
            strain array.
        ------------------------------------------------------------------------
        """
        dsigma = self.a*(self.T - self.Tm) - 3*self.b*(epsilon**2) + ((5*(self.b**2))/(4*self.a*(self.Ta - self.Tm)))*(epsilon**4)
        return dsigma
    
    def solution_stress_driven(self, sigma, epsilon_0 = 0.0115) -> np.ndarray[np.float64]:
        """
        ------------------------------------------------------------------------
        Obtain the numerical solution for the strain for a stress-driven 
        simulation using the Newton-Raphson method.
        ------------------------------------------------------------------------
        Parameters:
        - sigma: np.ndarray[np.float64]
            Stress array.  
        - epsilon_0: float
            Initial guess for the root.
        ------------------------------------------------------------------------
        Returns: np.ndarray[np.float64]
            The strain array corresponding to the given stress array.
        ------------------------------------------------------------------------
        """
        # Allocate array for strain results
        epsilon = np.zeros_like(sigma)
        # Loop thorugh all stress values and compute corresponding strain using Newton-Raphson method
        for j, s in enumerate(sigma):
            self.sigma_target = s
            try:
                epsilon[j] = newton_raphson(self.stress_for_root, 
                                            self.stress_derivative, 
                                            epsilon_0, 
                                            tol = 1e-6, 
                                            max_iter = 10000)
                epsilon_0 = epsilon[j]
            except ValueError:
                print(f"Did not converge for sigma={s}. Initial guess = {epsilon_0}")
                epsilon[j] = np.nan
        
        return epsilon    

def plot_single_behavior(axs: np.ndarray, i: int, T: float, 
                         epsilon: np.ndarray, sigma: np.ndarray, 
                         psi: np.ndarray, xlim: list = [-0.1, 0.1]) -> None:
    """
    ------------------------------------------------------------------------
    Auxiliary function to plot the free energy and stress-strain behavior
    for a single simulation type (strain-driven or stress-driven).
    ------------------------------------------------------------------------
    Parameters:
    - axs: np.ndarray
        Array of matplotlib axes objects for subplots.
    - i: int
        Index of the current subplot.
    - T: float
        Current temperature [ºC] or [K].
    - epsilon: np.ndarray
        Strain array.
    - sigma: np.ndarray
        Stress array.
    - psi: np.ndarray
        Free energy array.
    - xlim: list
        Limits for the x-axis in the plots. Default is [-0.1, 0.1].
    ------------------------------------------------------------------------
    Returns: None
    ------------------------------------------------------------------------
    """
    axs[0, i].plot(epsilon, psi, color = 'Green')
    axs[0, i].set_title(fr'$T={T}^{{\circ}}$C')
    axs[0, i].set_xlabel(r'$\epsilon$')
    axs[0, i].set_ylabel(r'$\psi \,[(J \times 10^6)/m^3]$')
    axs[0, i].set_xlim(xlim)
    
    axs[1, i].plot(epsilon, sigma, color = 'red')
    axs[1, i].set_title(fr'$T={T}^{{\circ}}$C')
    axs[1, i].set_xlabel(r'$\epsilon$')
    axs[1, i].set_ylabel(r'$\sigma \,[MPa]$')
    axs[1, i].set_xlim(xlim)
    
    for ax in axs.flatten():
        ax.axhline(0, color='gray', linestyle='-', linewidth=0.7)
        ax.axvline(0, color='gray', linestyle='-', linewidth=0.7)

def plot_both_behaviors(axs: np.ndarray, i: int, T: float, 
                        epsilon_1: np.ndarray, sigma_1: np.ndarray, 
                        epsilon_2: np.ndarray, sigma_2: np.ndarray, 
                        psi: np.ndarray, xlim: list = [-0.1, 0.1]) -> None:
    """
    ------------------------------------------------------------------------
    Auxiliary function to plot the free energy and stress-strain behavior
    for a single simulation type (strain-driven or stress-driven).
    ------------------------------------------------------------------------
    Parameters
    - axs: np.ndarray
        Array of matplotlib axes objects for subplots.
    - i: int
        Index of the current subplot.
    - T: float
        Current temperature [ºC] or [K].
    - epsilon_1: np.ndarray
        Strain array for case 1.
    - sigma_1: np.ndarray
        Stress array for case 1.
    - epsilon_2: np.ndarray
        Strain array for case 2.
    - sigma_2: np.ndarray
        Stress array for case 2.
    - psi: np.ndarray
        Free energy array.
    - xlim: list
        Limits for the x-axis in the plots. Default is [-0.1, 0.1].
    ------------------------------------------------------------------------
    Returns: None
    ------------------------------------------------------------------------
    """
    axs[0, i].plot(epsilon_1, psi, color = 'darkorange', lw = 1.7)
    axs[0, i].set_title(fr'$T={T}^{{\circ}}$C')
    axs[0, i].set_xlabel(r'$\epsilon$')
    axs[0, i].set_ylabel(r'$\psi \,[J] \times 10^6$')
    axs[0, i].set_xlim(xlim)
    
    axs[1, i].plot(epsilon_1, sigma_1, color = 'red', ls = '--', lw = 1.7, label='Strain-driven')
    axs[1, i].plot(epsilon_2, sigma_2, color = 'dimgrey', ls = '-', lw = 1.7, label='Stress-driven')
    axs[1, i].set_title(fr'$T={T}^{{\circ}}$C')
    axs[1, i].set_xlabel(r'$\epsilon$')
    axs[1, i].set_ylabel(r'$\sigma \,[MPa]$')
    axs[1, i].set_xlim(xlim)
    axs[1, i].legend()
    
    for ax in axs.flatten():
        ax.axhline(0, color='gray', linestyle='-', linewidth=0.7)
        ax.axvline(0, color='gray', linestyle='-', linewidth=0.7)

def save_example_fig(fig: mfig.Figure, axs: np.ndarray, dpi: int = 300) -> None: 
    """
    ----------------------------------------------------------------------------
    Adjusts subplots and saves a figure for example visualization.
    ----------------------------------------------------------------------------
    Parameters:
    - fig: matplotlib.figure.Figure (mfig.Figure)
        The figure object containing the subplots.
    - axs: np.ndarray
        Array of matplotlib axes objects corresponding to the subplots.
    - dpi: int, optional
        Resolution of the saved figure in dots per inch. Default is 300.
    ----------------------------------------------------------------------------
    Returns:
    - None
        The function modifies axes limits and saves the figure as 'example.png'.
    ----------------------------------------------------------------------------
    """
    axs[0, 0].set_ylim([-30, 40])
    axs[0, 1].set_ylim([-12, 40])
    axs[0, 2].set_ylim([-12, 40])
    
    xlim = 0.09
    for j in range(3):
        axs[0, j].set_xlim([-xlim, xlim])
        axs[1, j].set_xlim([-xlim, xlim])
            
    axs[1, 0].set_ylim([-2000, 2000])
    axs[1, 1].set_ylim([-1500, 1500])
    axs[1, 2].set_ylim([-750, 750])
    
    fig.savefig('example.png', dpi=dpi)

if __name__ == "__main__":
    
    # Program Input
    solution_type = "both" # Options: "strain-driven", "stress-driven", "both"
    
    # Define material parameters
    a = 1e3               #[MPa/ºC] 
    b = 40e6              #[MPa/ºC]
    Tm = 14.0             #[ºC]
    Ta = 40.0             #[ºC]
    T_list = [10, 25, 50] #[ºC]
            
    # Create 3 subplots with different temperatures
    fig, axs = plt.subplots(2, 3, figsize=(12, 7), layout='constrained')        
    
    for (i, T), ax in zip(enumerate(T_list), axs.flatten()):
        
        model = FalkSMAModel(a, b, Ta, Tm, T)
        
        if solution_type == "strain-driven":
            # Create loading profile
            epsilon_lim = 0.057
            epsilon = np.linspace(-epsilon_lim, epsilon_lim, 1000)   
            # Simulate and plot results
            sigma = model.stress(epsilon)
            psi = model.free_energy(epsilon)
            plot_single_behavior(axs, i, T, epsilon, sigma, psi)
            
        elif solution_type == "stress-driven":
            # Create loading profile
            sigma_lim = 1500
            sigma_up_1 = np.linspace(0, sigma_lim, 500)
            sigma_down = np.linspace(sigma_lim, -sigma_lim, 1000)
            sigma_up = np.linspace(-sigma_lim, sigma_lim, 1000)
            sigma = np.concatenate([sigma_up_1, sigma_down, sigma_up])   
            # Simulate and plot results
            epsilon = model.solution_stress_driven(sigma)            
            psi = model.free_energy(epsilon)
            plot_single_behavior(axs, i, T, epsilon, sigma, psi)
        else:
            # Strain-driven part
            epsilon_lim = 0.08
            epsilon_1 = np.linspace(-epsilon_lim, epsilon_lim, 1000)   
            sigma_1 = model.stress(epsilon_1)
            psi = model.free_energy(epsilon_1)
            
            # Stress-driven part
            sigma_lim = 1500
            sigma_up_1 = np.linspace(0, sigma_lim, 500)
            sigma_down = np.linspace(sigma_lim, -sigma_lim, 1000)
            sigma_up = np.linspace(-sigma_lim, sigma_lim, 1000)
            sigma_2 = np.concatenate([sigma_up_1, sigma_down, sigma_up])   
            
            epsilon_2 = model.solution_stress_driven(sigma_2, epsilon_0 = 0.2)                  
              
            plot_both_behaviors(axs, i, T, epsilon_1, sigma_1, epsilon_2, sigma_2, psi)
            
    save_example_fig(fig, axs, dpi = 600)    
    plt.show()