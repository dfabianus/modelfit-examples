import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
import os
import warnings

# Filter out specific RuntimeWarnings related to numerical issues
warnings.filterwarnings("ignore", category=RuntimeWarning, message="overflow encountered")
warnings.filterwarnings("ignore", category=RuntimeWarning, message="invalid value encountered")

# 1. Model function - ODE for biomass growth with dual substrates
def dual_substrate_growth_model(t, y, params):
    """
    ODE model for biomass growth with two competing substrates:
    dX/dt = mu1 * X + mu2 * X
    where:
    - mu1 = mu_max1 * S1 / (k1 + S1)  # High affinity substrate (e.g., glucose)
    - mu2 = mu_max2 * S2 / (k2 + S2)  # Low affinity substrate
    
    Parameters:
    -----------
    t : float
        Time point
    y : array-like
        State variables [X, S1, S2] where:
        - X is biomass
        - S1 is high-affinity substrate (e.g., glucose)
        - S2 is low-affinity substrate
    params : tuple
        (mu_max1, k1, mu_max2, k2) where:
        - mu_max1 is maximum growth rate on substrate 1
        - k1 is half-saturation constant for substrate 1
        - mu_max2 is maximum growth rate on substrate 2
        - k2 is half-saturation constant for substrate 2
    
    Returns:
    --------
    array-like
        Derivatives [dX/dt, dS1/dt, dS2/dt]
    """
    X, S1, S2 = y
    mu_max1, k1, mu_max2, k2 = params
    
    # Protect against division by zero or negative values
    S1 = max(S1, 1e-10)
    S2 = max(S2, 1e-10)
    
    # Monod equations for growth rates on each substrate
    mu1 = mu_max1 * S1 / (k1 + S1)
    mu2 = mu_max2 * S2 / (k2 + S2)
    
    # Total biomass growth rate (sum of growth on both substrates)
    dXdt = (mu1 + mu2) * X
    
    # Substrate consumption rates (assume yield coefficient = 1 for simplicity)
    dS1dt = -mu1 * X
    dS2dt = -mu2 * X
    
    return [dXdt, dS1dt, dS2dt]

# 2. Simulation function
def simulate_dual_substrate_growth(params, t_span, initial_conditions, t_eval=None):
    """
    Solve the ODE system for biomass growth with dual substrates over time
    
    Parameters:
    -----------
    params : tuple
        (mu_max1, k1, mu_max2, k2) parameters for the model
    t_span : tuple
        (t_start, t_end) time span for simulation
    initial_conditions : array-like
        [X0, S1_0, S2_0] initial biomass and substrate concentrations
    t_eval : array-like, optional
        Specific time points at which to evaluate the solution
        
    Returns:
    --------
    dict
        Solution containing time points and state variables
    """
    try:
        solution = solve_ivp(
            lambda t, y: dual_substrate_growth_model(t, y, params),
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
            'substrate1': solution.y[1],
            'substrate2': solution.y[2],
            'success': solution.success
        }
    except Exception as e:
        # Return a failed simulation with high error
        if t_eval is None:
            t_eval = np.linspace(t_span[0], t_span[1], 100)
        
        return {
            'time': t_eval,
            'biomass': np.ones_like(t_eval) * 1e6,  # Large error
            'substrate1': np.ones_like(t_eval) * 1e6,
            'substrate2': np.ones_like(t_eval) * 1e6,
            'success': False
        }

# 3. Optimization function
def fit_model(data_time, data_biomass, initial_substrate1, initial_substrate2, initial_guess=(0.5, 0.2, 0.3, 2.0)):
    """
    Fit the model parameters to experimental data
    
    Parameters:
    -----------
    data_time : array-like
        Time points of experimental data
    data_biomass : array-like
        Biomass measurements at each time point
    initial_substrate1 : float
        Initial concentration of high-affinity substrate
    initial_substrate2 : float
        Initial concentration of low-affinity substrate
    initial_guess : tuple
        Initial guess for (mu_max1, k1, mu_max2, k2)
        
    Returns:
    --------
    dict
        Optimized parameters and fit statistics
    """
    initial_biomass = data_biomass[0]
    initial_conditions = [initial_biomass, initial_substrate1, initial_substrate2]
    
    # Define the cost function (sum of squared errors)
    def objective_function(params):
        # Enforce parameter constraints
        if any(p <= 0 for p in params):
            return 1e10  # Large penalty for negative parameters
        
        # Simulate with current parameters
        simulation = simulate_dual_substrate_growth(
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
        regularization = 0.01 * (params[0]**2 + params[1]**2 + params[2]**2 + params[3]**2)
        
        return sse + regularization
    
    # Try multiple optimization attempts with different initial guesses
    best_result = None
    best_score = float('inf')
    
    # Define bounds for parameters
    # mu_max1, k1, mu_max2, k2
    bounds = [(0.1, 1.0), (0.01, 1.0), (0.05, 1.0), (0.5, 10.0)]
    
    # Try a few different initial guesses
    initial_guesses = [
        initial_guess,
        (0.4, 0.1, 0.25, 3.0),
        (0.6, 0.3, 0.35, 1.5)
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
            'mu_max1': initial_guess[0],
            'k1': initial_guess[1],
            'mu_max2': initial_guess[2],
            'k2': initial_guess[3],
            'sse': 1e10,
            'nrmse': 1.0,
            'success': False,
            'message': "All optimization attempts failed"
        }
    
    # Calculate NRMSE (Normalized Root Mean Square Error) with the best parameters
    optimal_params = best_result.x
    simulation = simulate_dual_substrate_growth(
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
        'mu_max1': optimal_params[0],
        'k1': optimal_params[1],
        'mu_max2': optimal_params[2],
        'k2': optimal_params[3],
        'sse': best_result.fun,
        'nrmse': nrmse,
        'success': best_result.success,
        'message': best_result.message
    }

# Generate synthetic data
def generate_example_data(true_params=(0.5, 0.2, 0.3, 2.0), noise_level=0.05):
    """
    Generate synthetic biomass growth data with dual substrates and noise
    
    Parameters:
    -----------
    true_params : tuple
        True values of (mu_max1, k1, mu_max2, k2) to generate data
    noise_level : float
        Relative noise level to add to the data
        
    Returns:
    --------
    tuple
        (time_points, biomass_data, initial_substrate1, initial_substrate2)
    """
    # Set initial conditions and time span
    initial_biomass = 0.1  # g/L
    initial_substrate1 = 10.0  # g/L (high affinity substrate)
    initial_substrate2 = 5.0  # g/L (low affinity substrate)
    initial_conditions = [initial_biomass, initial_substrate1, initial_substrate2]
    
    # Use appropriate time span to observe both substrate consumption phases
    t_span = (0, 30)  # hours
    
    # Generate clean data
    time_points = np.linspace(t_span[0], t_span[1], 40)
    clean_simulation = simulate_dual_substrate_growth(
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
    
    return time_points, noisy_biomass, initial_substrate1, initial_substrate2

# Main function
def main():
    print("Dual Substrate Growth Model")
    print("==========================")
    
    # Generate example data with dual substrates
    # Parameters: mu_max1, k1, mu_max2, k2
    # Set k1 << k2 to demonstrate substrate preference (diauxic growth)
    true_params = (0.5, 0.001, 0.1, 2.0)  # Lower k1 means higher affinity for substrate 1
    print(f"True parameters:")
    print(f"  Substrate 1 (high affinity): mu_max1 = {true_params[0]}, k1 = {true_params[1]}")
    print(f"  Substrate 2 (low affinity):  mu_max2 = {true_params[2]}, k2 = {true_params[3]}")
    
    time_points, biomass_data, initial_substrate1, initial_substrate2 = generate_example_data(true_params, noise_level=0.05)
    
    # Fit the model
    initial_guess = (0.4, 0.15, 0.25, 10.5)  # Initial guess for optimization
    print("\nOptimizing parameters, this may take a moment...")
    fit_results = fit_model(time_points, biomass_data, initial_substrate1, initial_substrate2, initial_guess)
    
    print("\nFit results:")
    print(f"  Substrate 1: mu_max1 = {fit_results['mu_max1']:.4f} h^-1, k1 = {fit_results['k1']:.4f} g/L")
    print(f"  Substrate 2: mu_max2 = {fit_results['mu_max2']:.4f} h^-1, k2 = {fit_results['k2']:.4f} g/L")
    print(f"  NRMSE = {fit_results['nrmse']:.4f}")
    print(f"  Success: {fit_results['success']}")
    
    # Simulate with the fitted parameters
    fitted_params = (fit_results['mu_max1'], fit_results['k1'], fit_results['mu_max2'], fit_results['k2'])
    initial_conditions = [biomass_data[0], initial_substrate1, initial_substrate2]
    
    # Use more points for smoother curves
    fitted_simulation = simulate_dual_substrate_growth(
        fitted_params,
        (time_points[0], time_points[-1]),
        initial_conditions,
        t_eval=np.linspace(time_points[0], time_points[-1], 200)
    )
    
    # Also simulate with the true parameters for comparison
    true_simulation = simulate_dual_substrate_growth(
        true_params,
        (time_points[0], time_points[-1]),
        initial_conditions,
        t_eval=np.linspace(time_points[0], time_points[-1], 200)
    )
    
    # Plot the results
    plt.figure(figsize=(12, 12))
    
    # Plot biomass data and fit
    plt.subplot(3, 1, 1)
    plt.scatter(time_points, biomass_data, color='blue', label='Experimental data')
    plt.plot(fitted_simulation['time'], fitted_simulation['biomass'], 'r-', label='Model fit')
    plt.plot(true_simulation['time'], true_simulation['biomass'], 'g--', label='True model')
    plt.xlabel('Time (h)')
    plt.ylabel('Biomass (g/L)')
    plt.legend()
    plt.title('Biomass Growth with Dual Substrates - Model Fitting')
    plt.grid(True)
    
    # Plot substrate consumption
    plt.subplot(3, 1, 2)
    plt.plot(fitted_simulation['time'], fitted_simulation['substrate1'], 'r-', label='Substrate 1 (high affinity) - Fitted')
    plt.plot(fitted_simulation['time'], fitted_simulation['substrate2'], 'b-', label='Substrate 2 (low affinity) - Fitted')
    plt.plot(true_simulation['time'], true_simulation['substrate1'], 'g--', label='Substrate 1 - True')
    plt.plot(true_simulation['time'], true_simulation['substrate2'], 'c--', label='Substrate 2 - True')
    plt.xlabel('Time (h)')
    plt.ylabel('Substrate (g/L)')
    plt.legend()
    plt.title('Substrate Consumption')
    plt.grid(True)
    
    # Plot specific growth rates on each substrate
    plt.subplot(3, 1, 3)
    
    # Calculate growth rates for fitted model
    fitted_mu_max1, fitted_k1, fitted_mu_max2, fitted_k2 = fitted_params
    fitted_s1 = fitted_simulation['substrate1']
    fitted_s2 = fitted_simulation['substrate2']
    
    fitted_mu1 = fitted_mu_max1 * fitted_s1 / (fitted_k1 + fitted_s1)
    fitted_mu2 = fitted_mu_max2 * fitted_s2 / (fitted_k2 + fitted_s2)
    fitted_mu_total = fitted_mu1 + fitted_mu2
    
    # Calculate growth rates for true model
    true_mu_max1, true_k1, true_mu_max2, true_k2 = true_params
    true_s1 = true_simulation['substrate1']
    true_s2 = true_simulation['substrate2']
    
    true_mu1 = true_mu_max1 * true_s1 / (true_k1 + true_s1)
    true_mu2 = true_mu_max2 * true_s2 / (true_k2 + true_s2)
    true_mu_total = true_mu1 + true_mu2
    
    plt.plot(fitted_simulation['time'], fitted_mu1, 'r-', label='Growth rate on S1 - Fitted')
    plt.plot(fitted_simulation['time'], fitted_mu2, 'b-', label='Growth rate on S2 - Fitted')
    plt.plot(fitted_simulation['time'], fitted_mu_total, 'k-', label='Total growth rate - Fitted')
    plt.plot(true_simulation['time'], true_mu_total, 'g--', label='Total growth rate - True')
    
    plt.xlabel('Time (h)')
    plt.ylabel('Specific Growth Rate (h^-1)')
    plt.legend()
    plt.title('Specific Growth Rates')
    plt.grid(True)
    
    plt.tight_layout()
    
    # Get the project root directory (2 levels up from this script)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    figures_dir = os.path.join(project_root, 'figures')
    
    # Create figures directory if it doesn't exist
    os.makedirs(figures_dir, exist_ok=True)
    
    # Save the plot to the figures directory using an absolute path
    plt.savefig(os.path.join(figures_dir, 'dual_substrate_growth_fit.png'), dpi=300)
    plt.show()

if __name__ == "__main__":
    main() 