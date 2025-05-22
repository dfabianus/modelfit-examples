import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
import os
import warnings

# Filter out specific RuntimeWarnings related to numerical issues
warnings.filterwarnings("ignore", category=RuntimeWarning, message="overflow encountered")
warnings.filterwarnings("ignore", category=RuntimeWarning, message="invalid value encountered")

# 1. Model function - ODE for biomass growth with death phase
def biomass_growth_death_model(t, y, params):
    """
    ODE model for biomass growth with death phase:
    dX/dt = mu * X - kd * X
    where mu = mu_max * S / (k + S)
    
    Parameters:
    -----------
    t : float
        Time point
    y : array-like
        State variables [X, S] where X is biomass and S is substrate
    params : tuple
        (mu_max, k, kd) where:
        - mu_max is maximum growth rate
        - k is half-saturation constant
        - kd is death rate constant
    
    Returns:
    --------
    array-like
        Derivatives [dX/dt, dS/dt]
    """
    X, S = y
    mu_max, k, kd = params
    
    # Protect against division by zero or negative values
    S = max(S, 1e-10)
    
    # Monod equation for growth rate
    mu = mu_max * S / (k + S)
    
    # Biomass growth rate with death term
    dXdt = mu * X - kd * X
    
    # Substrate consumption rate (assume yield coefficient = 1 for simplicity)
    # Only consumed for growth, not for maintenance during death phase
    dSdt = -mu * X
    
    return [dXdt, dSdt]

# 2. Simulation function
def simulate_growth_death(params, t_span, initial_conditions, t_eval=None):
    """
    Solve the ODE system for biomass growth with death phase over time
    
    Parameters:
    -----------
    params : tuple
        (mu_max, k, kd) parameters for the model
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
    try:
        solution = solve_ivp(
            lambda t, y: biomass_growth_death_model(t, y, params),
            t_span,
            initial_conditions,
            t_eval=t_eval,
            method='RK45',
            rtol=1e-6,
            atol=1e-8
        )
        
        return {
            'time': solution.t,
            'biomass': solution.y[0],
            'substrate': solution.y[1],
            'success': solution.success
        }
    except Exception as e:
        # Return a failed simulation with high error
        if t_eval is None:
            t_eval = np.linspace(t_span[0], t_span[1], 100)
        
        return {
            'time': t_eval,
            'biomass': np.ones_like(t_eval) * 1e6,  # Large error
            'substrate': np.ones_like(t_eval) * 1e6,
            'success': False
        }

# 3. Optimization function
def fit_model(data_time, data_biomass, initial_substrate, initial_guess=(0.5, 1.0, 0.05)):
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
        Initial guess for (mu_max, k, kd)
        
    Returns:
    --------
    dict
        Optimized parameters and fit statistics
    """
    initial_biomass = data_biomass[0]
    initial_conditions = [initial_biomass, initial_substrate]
    
    # Define the cost function (sum of squared errors)
    def objective_function(params):
        # Enforce parameter constraints
        if any(p <= 0 for p in params):
            return 1e10  # Large penalty for negative parameters
        
        # Simulate with current parameters
        simulation = simulate_growth_death(
            params, 
            (data_time[0], data_time[-1]), 
            initial_conditions,
            t_eval=data_time
        )
        
        # Check if simulation was successful
        if not simulation.get('success', False):
            return 1e10  # Large penalty for failed simulations
        
        # Calculate sum of squared errors with numerical stability
        simulated_biomass = simulation['biomass']
        
        # Handle potential NaN or Inf values
        if np.any(np.isnan(simulated_biomass)) or np.any(np.isinf(simulated_biomass)):
            return 1e10
        
        # Calculate errors safely
        errors = simulated_biomass - data_biomass
        sse = np.sum(errors**2)
        
        # Add regularization to avoid extreme parameter values
        regularization = 0.01 * (params[0]**2 + params[1]**2 + 10*params[2]**2)
        
        return sse + regularization
    
    # Try multiple optimization attempts with different initial guesses
    best_result = None
    best_score = float('inf')
    
    # Define bounds for parameters - more constrained to avoid numerical issues
    bounds = [(0.1, 2.0), (0.1, 5.0), (0.001, 0.2)]  # Bounds for mu_max, k, and kd
    
    # Try a few different initial guesses
    initial_guesses = [
        initial_guess,
        (0.4, 0.8, 0.03),
        (0.6, 1.2, 0.07)
    ]
    
    for guess in initial_guesses:
        try:
            result = minimize(
                objective_function,
                guess,
                bounds=bounds,
                method='L-BFGS-B',
                options={'ftol': 1e-8, 'gtol': 1e-6, 'maxiter': 500}
            )
            
            if result.fun < best_score:
                best_score = result.fun
                best_result = result
        except Exception as e:
            print(f"Optimization attempt failed: {e}")
            continue
    
    if best_result is None:
        # If all optimizations failed, return the initial guess with high error
        return {
            'mu_max': initial_guess[0],
            'k': initial_guess[1],
            'kd': initial_guess[2],
            'sse': 1e10,
            'nrmse': 1.0,
            'success': False,
            'message': "All optimization attempts failed"
        }
    
    # Calculate NRMSE (Normalized Root Mean Square Error) with the best parameters
    optimal_params = best_result.x
    simulation = simulate_growth_death(
        optimal_params, 
        (data_time[0], data_time[-1]), 
        initial_conditions,
        t_eval=data_time
    )
    simulated_biomass = simulation['biomass']
    
    # Safe calculation of errors
    squared_errors = np.square(simulated_biomass - data_biomass)
    mse = np.mean(squared_errors)
    rmse = np.sqrt(mse)
    nrmse = rmse / (np.max(data_biomass) - np.min(data_biomass))
    
    return {
        'mu_max': optimal_params[0],
        'k': optimal_params[1],
        'kd': optimal_params[2],
        'sse': best_result.fun,
        'nrmse': nrmse,
        'success': best_result.success,
        'message': best_result.message
    }

# Generate synthetic data
def generate_example_data(true_params=(0.5, 1.0, 0.05), noise_level=0.05):
    """
    Generate synthetic biomass growth data with death phase and noise
    
    Parameters:
    -----------
    true_params : tuple
        True values of (mu_max, k, kd) to generate data
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
    
    # Use longer time span to observe death phase
    t_span = (0, 48)  # hours
    
    # Generate clean data
    time_points = np.linspace(t_span[0], t_span[1], 30)
    clean_simulation = simulate_growth_death(
        true_params,
        t_span,
        initial_conditions,
        t_eval=time_points
    )
    
    # Add noise to biomass data
    clean_biomass = clean_simulation['biomass']
    np.random.seed(42)  # For reproducibility
    
    # Scale noise by local biomass value to avoid negative values
    relative_noise = np.random.normal(0, noise_level, len(clean_biomass))
    noisy_biomass = clean_biomass * (1 + relative_noise)
    
    # Ensure no negative values
    noisy_biomass = np.maximum(noisy_biomass, 0.001)
    
    return time_points, noisy_biomass, initial_substrate

# Main function
def main():
    print("Biomass Growth Model with Death Phase")
    print("=====================================")
    
    # Generate example data with death phase
    # Adjusted parameters for better numerical stability
    true_params = (0.5, 0.8, 0.05)  # mu_max = 0.5 h^-1, k = 0.8 g/L, kd = 0.05 h^-1
    print(f"True parameters: mu_max = {true_params[0]}, k = {true_params[1]}, kd = {true_params[2]}")
    
    time_points, biomass_data, initial_substrate = generate_example_data(true_params, noise_level=0.05)
    
    # Fit the model
    initial_guess = (0.3, 0.5, 0.02)  # Initial guess for optimization
    print("\nOptimizing parameters, this may take a moment...")
    fit_results = fit_model(time_points, biomass_data, initial_substrate, initial_guess)
    
    print("\nFit results:")
    print(f"mu_max = {fit_results['mu_max']:.4f} h^-1")
    print(f"k = {fit_results['k']:.4f} g/L")
    print(f"kd = {fit_results['kd']:.4f} h^-1")
    print(f"NRMSE = {fit_results['nrmse']:.4f}")
    print(f"Success: {fit_results['success']}")
    
    # Simulate with the fitted parameters
    fitted_params = (fit_results['mu_max'], fit_results['k'], fit_results['kd'])
    initial_conditions = [biomass_data[0], initial_substrate]
    
    # Use more points for smoother curves
    fitted_simulation = simulate_growth_death(
        fitted_params,
        (time_points[0], time_points[-1]),
        initial_conditions,
        t_eval=np.linspace(time_points[0], time_points[-1], 200)
    )
    
    # Also simulate with the true parameters for comparison
    true_simulation = simulate_growth_death(
        true_params,
        (time_points[0], time_points[-1]),
        initial_conditions,
        t_eval=np.linspace(time_points[0], time_points[-1], 200)
    )
    
    # Plot the results
    plt.figure(figsize=(12, 10))
    
    # Plot biomass data and fit
    plt.subplot(3, 1, 1)
    plt.scatter(time_points, biomass_data, color='blue', label='Experimental data')
    plt.plot(fitted_simulation['time'], fitted_simulation['biomass'], 'r-', label='Model fit')
    plt.plot(true_simulation['time'], true_simulation['biomass'], 'g--', label='True model')
    plt.xlabel('Time (h)')
    plt.ylabel('Biomass (g/L)')
    plt.legend()
    plt.title('Biomass Growth with Death Phase - Model Fitting')
    plt.grid(True)
    
    # Plot substrate consumption
    plt.subplot(3, 1, 2)
    plt.plot(fitted_simulation['time'], fitted_simulation['substrate'], 'r-', label='Fitted model')
    plt.plot(true_simulation['time'], true_simulation['substrate'], 'g--', label='True model')
    plt.xlabel('Time (h)')
    plt.ylabel('Substrate (g/L)')
    plt.legend()
    plt.title('Substrate Consumption')
    plt.grid(True)
    
    # Plot specific growth rate
    plt.subplot(3, 1, 3)
    
    # Calculate net specific growth rate (growth - death)
    fitted_mu_max, fitted_k, fitted_kd = fitted_params
    true_mu_max, true_k, true_kd = true_params
    
    # Calculate effective growth rates
    fitted_substrate = fitted_simulation['substrate']
    fitted_growth_rate = fitted_mu_max * fitted_substrate / (fitted_k + fitted_substrate) - fitted_kd
    
    true_substrate = true_simulation['substrate']
    true_growth_rate = true_mu_max * true_substrate / (true_k + true_substrate) - true_kd
    
    plt.plot(fitted_simulation['time'], fitted_growth_rate, 'r-', label='Fitted net growth rate')
    plt.plot(true_simulation['time'], true_growth_rate, 'g--', label='True net growth rate')
    plt.axhline(y=0, color='k', linestyle='-', alpha=0.3)
    plt.xlabel('Time (h)')
    plt.ylabel('Net Specific Growth Rate (h^-1)')
    plt.legend()
    plt.title('Net Specific Growth Rate (μ - kd)')
    plt.grid(True)
    
    plt.tight_layout()
    
    # Get the project root directory (2 levels up from this script)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    figures_dir = os.path.join(project_root, 'figures')
    
    # Create figures directory if it doesn't exist
    os.makedirs(figures_dir, exist_ok=True)
    
    # Save the plot to the figures directory using an absolute path
    plt.savefig(os.path.join(figures_dir, 'growth_death_model_fit.png'), dpi=300)
    plt.show()

if __name__ == "__main__":
    main() 