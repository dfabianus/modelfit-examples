import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
import random
import os

# 1. Model function - ODE for biomass growth
def biomass_growth_model(t, y, params):
    """
    ODE model for biomass growth: dX/dt = mu * X
    where mu = mu_max * S / (k + S)
    
    Parameters:
    -----------
    t : float
        Time point
    y : array-like
        State variables [X, S] where X is biomass and S is substrate
    params : tuple
        (mu_max, k) where mu_max is maximum growth rate and k is half-saturation constant
    
    Returns:
    --------
    array-like
        Derivatives [dX/dt, dS/dt]
    """
    X, S = y
    mu_max, k = params
    
    # Monod equation for growth rate
    mu = mu_max * S / (k + S)
    
    # Biomass growth rate
    dXdt = mu * X
    
    # Substrate consumption rate (assume yield coefficient = 1 for simplicity)
    dSdt = -dXdt
    
    return [dXdt, dSdt]

# 2. Simulation function
def simulate_growth(params, t_span, initial_conditions, t_eval=None):
    """
    Solve the ODE system for biomass growth over time
    
    Parameters:
    -----------
    params : tuple
        (mu_max, k) where mu_max is maximum growth rate and k is half-saturation constant
    t_span : tuple
        (t_start, t_end) time span for simulation
    initial_conditions : array-like
        [X0, S0] initial biomass and substrate concentrations
    t_eval : array-like, optional
        Specific time points at which to evaluate the solution
        
    Returns:
    --------
    dict
        Solution containing time points and state variables
    """
    solution = solve_ivp(
        lambda t, y: biomass_growth_model(t, y, params),
        t_span,
        initial_conditions,
        t_eval=t_eval,
        method='RK45'
    )
    
    return {
        'time': solution.t,
        'biomass': solution.y[0],
        'substrate': solution.y[1]
    }

# 3. Optimization function
def fit_model(data_time, data_biomass, initial_substrate, initial_guess=(0.5, 1.0)):
    """
    Fit the model parameters to experimental data
    
    Parameters:
    -----------
    data_time : array-like
        Time points of experimental data
    data_biomass : array-like
        Biomass measurements at each time point
    initial_substrate : float
        Initial substrate concentration
    initial_guess : tuple
        Initial guess for (mu_max, k)
        
    Returns:
    --------
    dict
        Optimized parameters and fit statistics
    """
    initial_biomass = data_biomass[0]
    initial_conditions = [initial_biomass, initial_substrate]
    
    # Define the cost function (sum of squared errors)
    def objective_function(params):
        # Simulate with current parameters
        simulation = simulate_growth(
            params, 
            (data_time[0], data_time[-1]), 
            initial_conditions,
            t_eval=data_time
        )
        
        # Calculate sum of squared errors
        simulated_biomass = simulation['biomass']
        squared_errors = (simulated_biomass - data_biomass) ** 2
        sse = np.sum(squared_errors)
        
        return sse
    
    # Optimize parameters
    bounds = [(0.01, 10.0), (0.01, 10.0)]  # Bounds for mu_max and k
    result = minimize(
        objective_function,
        initial_guess,
        bounds=bounds,
        method='L-BFGS-B'
    )
    
    # Calculate NRMSE (Normalized Root Mean Square Error)
    optimal_params = result.x
    simulation = simulate_growth(
        optimal_params, 
        (data_time[0], data_time[-1]), 
        initial_conditions,
        t_eval=data_time
    )
    simulated_biomass = simulation['biomass']
    mse = np.mean((simulated_biomass - data_biomass) ** 2)
    rmse = np.sqrt(mse)
    nrmse = rmse / (np.max(data_biomass) - np.min(data_biomass))
    
    return {
        'mu_max': optimal_params[0],
        'k': optimal_params[1],
        'sse': result.fun,
        'nrmse': nrmse,
        'success': result.success,
        'message': result.message
    }

# Generate synthetic data
def generate_example_data(true_params=(0.5, 2.0), noise_level=0.05):
    """
    Generate synthetic biomass growth data with noise
    
    Parameters:
    -----------
    true_params : tuple
        True values of (mu_max, k) to generate data
    noise_level : float
        Relative noise level to add to the data
        
    Returns:
    --------
    tuple
        (time_points, biomass_data, initial_substrate)
    """
    # Set initial conditions and time span
    initial_biomass = 0.1  # g/L
    initial_substrate = 10.0  # g/L
    initial_conditions = [initial_biomass, initial_substrate]
    t_span = (0, 24)  # hours
    
    # Generate clean data
    time_points = np.linspace(t_span[0], t_span[1], 20)
    clean_simulation = simulate_growth(
        true_params,
        t_span,
        initial_conditions,
        t_eval=time_points
    )
    
    # Add noise to biomass data
    clean_biomass = clean_simulation['biomass']
    np.random.seed(42)  # For reproducibility
    noise = np.random.normal(0, noise_level * clean_biomass.mean(), len(clean_biomass))
    noisy_biomass = clean_biomass + noise
    
    # Ensure no negative values
    noisy_biomass = np.maximum(noisy_biomass, 0.001)
    
    return time_points, noisy_biomass, initial_substrate

# Main function
def main():
    # Generate example data
    true_params = (0.3, 1.5)  # mu_max = 0.3 h^-1, k = 1.5 g/L
    print(f"True parameters: mu_max = {true_params[0]}, k = {true_params[1]}")
    
    time_points, biomass_data, initial_substrate = generate_example_data(true_params, noise_level=0.05)
    
    # Fit the model
    initial_guess = (0.1, 1.0)  # Initial guess for optimization
    fit_results = fit_model(time_points, biomass_data, initial_substrate, initial_guess)
    
    print("\nFit results:")
    print(f"mu_max = {fit_results['mu_max']:.4f} h^-1")
    print(f"k = {fit_results['k']:.4f} g/L")
    print(f"NRMSE = {fit_results['nrmse']:.4f}")
    
    # Simulate with the fitted parameters
    fitted_params = (fit_results['mu_max'], fit_results['k'])
    initial_conditions = [biomass_data[0], initial_substrate]
    fitted_simulation = simulate_growth(
        fitted_params,
        (time_points[0], time_points[-1]),
        initial_conditions,
        t_eval=np.linspace(time_points[0], time_points[-1], 100)
    )
    
    # Plot the results
    plt.figure(figsize=(12, 8))
    
    # Plot biomass data and fit
    plt.subplot(2, 1, 1)
    plt.scatter(time_points, biomass_data, color='blue', label='Experimental data')
    plt.plot(fitted_simulation['time'], fitted_simulation['biomass'], 'r-', label='Model fit')
    plt.xlabel('Time (h)')
    plt.ylabel('Biomass (g/L)')
    plt.legend()
    plt.title('Biomass Growth Model Fitting')
    plt.grid(True)
    
    # Plot substrate consumption
    plt.subplot(2, 1, 2)
    plt.plot(fitted_simulation['time'], fitted_simulation['substrate'], 'g-', label='Substrate')
    plt.xlabel('Time (h)')
    plt.ylabel('Substrate (g/L)')
    plt.legend()
    plt.title('Substrate Consumption')
    plt.grid(True)
    
    plt.tight_layout()
    
    # Get the project root directory (2 levels up from this script)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    figures_dir = os.path.join(project_root, 'figures')
    
    # Create figures directory if it doesn't exist
    os.makedirs(figures_dir, exist_ok=True)
    
    # Save the plot to the figures directory using an absolute path
    plt.savefig(os.path.join(figures_dir, 'monod_growth_fit.png'), dpi=300)
    plt.show()

if __name__ == "__main__":
    main() 