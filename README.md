# Biomass Growth Model Fitting

This project demonstrates how to fit bioprocess models to time-series data of biomass growth. The implementation includes examples of simple Monod growth and a more sophisticated model that incorporates cell death.

## Overview

The implementation includes two examples, each with three core components:

### 1. Monod Growth Model

The basic model describes biomass growth following Monod kinetics:

- **Model Function**: Defines the ordinary differential equation (ODE) for biomass growth:
  - $\frac{dX}{dt} = \mu \cdot X$
  - Where $\mu = \mu_{max} \cdot \frac{S}{K_s + S}$ (Monod equation)
  - $X$ = biomass concentration
  - $S$ = substrate concentration
  - $\mu_{max}$ = maximum specific growth rate
  - $K_s$ = half-saturation constant

- **Simulation Function**: Uses `scipy.integrate.solve_ivp` to solve the ODE system over a specified time span.

- **Optimization Function**: Minimizes the cost function (sum of squared errors) by fitting the model parameters ($\mu_{max}$ and $K_s$) to experimental data.

### 2. Growth-Death Model

A more sophisticated model that accounts for cell death, showing biomass decline after peak growth:

- **Model Function**: Extends the basic model with a death term:
  - $\frac{dX}{dt} = \mu \cdot X - k_d \cdot X$
  - Where $\mu = \mu_{max} \cdot \frac{S}{K_s + S}$
  - $k_d$ = death rate constant

- **Simulation Function**: Solves the extended ODE system including the death kinetics.

- **Optimization Function**: Fits three parameters ($\mu_{max}$, $K_s$, and $k_d$) to experimental data showing growth and decline phases.

## Model Fitting Concept

The model fitting process follows an iterative optimization approach to find the best parameters that describe the experimental data:

```mermaid
flowchart TD
    A[Experimental Data] --> B[Initial Parameter Guess]
    B --> C[Simulate Model]
    C --> D[Calculate Error]
    D --> E{Error Acceptable?}
    E -->|No| F[Update Parameters]
    F --> C
    E -->|Yes| G[Final Parameters]
    
    H[Cost Function] --> D
    I[ODE Model] --> C
    J[Optimization Algorithm] --> F
    
    classDef data fill:#e1f5fe,stroke:#0288d1,stroke-width:2px
    classDef process fill:#fff9c4,stroke:#fbc02d,stroke-width:2px
    classDef decision fill:#f8bbd0,stroke:#c2185b,stroke-width:2px
    classDef model fill:#d1c4e9,stroke:#512da8,stroke-width:2px
    
    class A,G data
    class B,C,D,F process
    class E decision
    class H,I,J model
```

### Key Components

1. **Experimental Data**: Time-series measurements of biomass concentration and substrate levels.

2. **ODE Model**: Mathematical representation of the biological system:
   - Simple Monod model: $\frac{dX}{dt} = \mu_{max} \cdot \frac{S}{K_s + S} \cdot X$
   - Growth-death model: $\frac{dX}{dt} = \mu_{max} \cdot \frac{S}{K_s + S} \cdot X - k_d \cdot X$

3. **Simulation**: Numerical solution of the ODE system using initial conditions and current parameter estimates.

4. **Cost Function**: Quantifies the difference between model predictions and experimental data:
   - Sum of Squared Errors (SSE): $\sum_{i=1}^{n} (y_{pred,i} - y_{exp,i})^2$
   - Normalized Root Mean Square Error (NRMSE): $\frac{\sqrt{\frac{1}{n}\sum_{i=1}^{n} (y_{pred,i} - y_{exp,i})^2}}{y_{max} - y_{min}}$

5. **Optimization Algorithm**: Updates parameter values to minimize the cost function (e.g., L-BFGS-B algorithm).

6. **Convergence**: The process continues until the error reaches an acceptable level or the algorithm converges.

## Results

### Monod Growth Model

![Monod Growth Model Results](figures/monod_growth_fit.png)

The simple Monod growth model shows biomass accumulation as substrate is consumed, following classic Monod kinetics.

### Growth-Death Model

![Growth-Death Model Results](figures/growth_death_model_fit.png)

The growth-death model demonstrates a more complex behavior where biomass initially grows but then declines after substrate depletion, due to the death rate constant ($k_d$). The net specific growth rate plot shows when growth transitions to decline.

## Directory Structure

```
modelfit-examples/
├── examples/
│   ├── biomass_growth_model.py         # Simple Monod growth model
│   └── biomass_growth_death_model.py   # Growth model with death phase
├── figures/
│   ├── monod_growth_fit.png            # Output from simple model
│   └── growth_death_model_fit.png      # Output from death phase model
├── main.py                             # Runner for all examples
├── README.md
└── pyproject.toml
```

## Requirements

- Python >= 3.13
- NumPy
- SciPy
- Matplotlib

## Usage

This project uses `uv` for package and environment management. To get started:

### Setup

First, sync the dependencies:

```bash
# Install dependencies from pyproject.toml
uv sync
```

### Running Examples

#### Run All Examples at Once

The simplest way to run all examples is to use the main.py file:

```bash
# Run all examples sequentially
uv run main.py
```

This will execute all examples one after another and save their output plots to the `figures/` directory.

#### Run Individual Examples

Alternatively, you can run each example individually:

```bash
# Simple Monod growth model
uv run examples/biomass_growth_model.py

# Growth model with death phase
uv run examples/biomass_growth_death_model.py
```

If you prefer using Python directly:

```bash
# Simple Monod growth model
python examples/biomass_growth_model.py

# Growth model with death phase
python examples/biomass_growth_death_model.py
```

Each script will:
1. Generate synthetic biomass growth data with some noise
2. Fit the model parameters to this data
3. Plot the results showing both the original data points and the fitted model
4. Save the plot to the `figures/` directory

## Output

### Monod Growth Model
- The true parameters used to generate the synthetic data ($\mu_{max}$ and $K_s$)
- The fitted parameters ($\mu_{max}$ and $K_s$)
- The normalized root mean square error (NRMSE) of the fit
- A plot showing biomass growth and substrate consumption

### Growth-Death Model
- The true parameters used to generate the synthetic data ($\mu_{max}$, $K_s$, and $k_d$)
- The fitted parameters ($\mu_{max}$, $K_s$, and $k_d$)
- The normalized root mean square error (NRMSE) of the fit
- A plot showing biomass growth, decline, substrate consumption, and net growth rate

## Model Equations

### Monod Growth Model

$$\frac{dX}{dt} = \mu \cdot X$$
$$\frac{dS}{dt} = -\frac{dX}{dt}$$  (assuming yield coefficient = 1 for simplicity)

Where:
$$\mu = \mu_{max} \cdot \frac{S}{K_s + S}$$

### Growth-Death Model

$$\frac{dX}{dt} = \mu \cdot X - k_d \cdot X$$
$$\frac{dS}{dt} = -\mu \cdot X$$  (substrate consumed only for growth, not during death)

Where:
$$\mu = \mu_{max} \cdot \frac{S}{K_s + S}$$
