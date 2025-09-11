# Lagrange-Multiplier.py
A Python utility for solving **constrained optimization problems** using the **Lagrange Multiplier method**.  
This tool leverages [SymPy](https://www.sympy.org/) to symbolically compute solutions to optimization problems of the form:  

min f(x1, x2, ..., xn)   subject to   g(x1, x2, ..., xn) = 0

---

## ✨ Features
- Symbolic optimization with **Lagrange multipliers**
- Supports multiple variables
- Returns solutions including the **Lagrange multiplier**
- Easy to extend for more constraints

---

## ⚙️ Function Signature

```python
lagrange_multiplier(objective, constraint, variables)
```

# Least_Action_SHO
Principle of Least Action — Harmonic Oscillator:

This project applies the **principle of Least Action** to obtain the path that minimizes the action relative to all other possible paths between the same points, for a **harmonic oscillator system** with potential:

\[
V(x) = 1/2 kx^2
\]

The simulation uses a **Feynman Path Integral**–style approach:  
- Multiple possible paths between two fixed points are generated.  
- Each path contributes a **complex amplitude** proportional to \(\exp(iS/\hbar)\), where \(S\) is the classical action.  
- The **classical path** emerges as the path of least action, dominating the sum.  

---

## ✨ Features
- Simulates **100 perturbed paths** between fixed endpoints.  
- Calculates **kinetic, potential, and Lagrangian** for each path.  
- Evaluates the **action integral** using Simpson’s rule.  
- Computes **quantum probability amplitudes**.  
- Plots all sampled paths along with the **classical least-action path**.  

---

# 

## ⚙️ Requirements
Install the dependencies with:
```python
pip install numpy matplotlib scipy symmpy 
```

# Advanced Lagrange Multiplier Solver

A comprehensive Python implementation for solving constrained optimization problems using the method of Lagrange multipliers with advanced numerical techniques, visualization, and analysis capabilities.

## Features

- **Multiple Constraint Handling**: Solve problems with multiple equality constraints
- **Dual Solution Methods**: Both symbolic (exact) and numerical (approximate) solving approaches
- **Interactive Visualization**: 2D contour plots and 3D surface visualizations
- **Solution Verification**: Automatic verification against KKT conditions
- **Comprehensive Analysis**: Detailed solution analysis with constraint violation checking
- **Export Capabilities**: Export results to JSON and HTML report formats
- **Error Handling**: Robust error handling and logging
- **Batch Processing**: Solve multiple optimization problems in batch mode
- **Sensitivity Analysis**: Analyze how solutions change with constraint parameters
- **Performance Benchmarking**: Compare different solution methods

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/advanced-lagrange-solver.git
cd advanced-lagrange-solver

# Install required dependencies
pip install numpy sympy matplotlib scipy pandas seaborn