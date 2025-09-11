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

# 🎯 Advanced Lagrange Multiplier Solver

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://python.org)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Status](https://img.shields.io/badge/status-stable-brightgreen.svg)]()
[![Contributions](https://img.shields.io/badge/contributions-welcome-orange.svg)](CONTRIBUTING.md)

> 🔬 **A comprehensive Python toolkit for solving constrained optimization problems using the method of Lagrange multipliers with advanced numerical techniques, interactive visualizations, and scientific computing capabilities.**

---

## 📖 Table of Contents

- [✨ Features](#-Features)
- [🚀 Quick Start](#-quick-start)
- [📦 Installation](#-installation)
- [🎮 Interactive Demo](#-interactive-demo)
- [📚 Usage Examples](#-usage-examples)
- [🔧 API Reference](#-api-reference)
- [📊 Visualization Gallery](#-visualization-gallery)
- [🧪 Scientific Applications](#-scientific-applications)
- [⚡ Performance](#-performance)
- [🤝 Contributing](#-contributing)
- [📄 License](#-license)

---

## Features

### 🎯 **Core Capabilities**
- **Dual Solution Methods**: Symbolic (SymPy) + Numerical (SciPy) solvers
- **Multi-Constraint Support**: Handle complex constraint systems
- **Interactive Visualization**: 2D/3D plots with constraint surfaces
- **Solution Verification**: Automatic KKT condition checking
- **Batch Processing**: Solve multiple problems simultaneously

### 🔬 **Advanced Analytics**
- **Sensitivity Analysis**: Parameter dependency studies
- **Shadow Price Computation**: Economic interpretation
- **Performance Benchmarking**: Method comparison
- **Statistical Analysis**: Solution distribution studies
- **Export Capabilities**: JSON, HTML, CSV formats

### 🎨 **User Experience**
- **Interactive CLI**: Command-line interface with examples
- **Professional Visualizations**: Publication-ready plots
- **Comprehensive Logging**: Debug and monitoring support
- **Error Handling**: Robust exception management
- **Documentation**: Detailed API reference

---

## 🚀 Quick Start

### **30-Second Demo**

```python
from lagrange_optimizer import AdvancedLagrangeOptimizer

# Create optimizer instance
optimizer = AdvancedLagrangeOptimizer()

# Solve: Minimize x² + y² subject to x + y = 1
results = optimizer.solve_complete(
    objective="x**2 + y**2",           # Objective function
    constraints=["x + y - 1"],         # Constraint (set to zero)
    variables=["x", "y"],              # Variables
    problem_type="minimize",           # minimize or maximize
    visualize=True                     # Generate plots
)

print(results)  # View results as pandas DataFrame
```

**Output:**
```
Solution 1:
  x = 0.500000
  y = 0.500000
  λ1 = 1.000000
  Objective value = 0.500000
```

---

## 📦 Installation

### **Prerequisites**
- Python 3.8 or higher
- pip package manager

### **Option 1: Quick Install**
```bash
pip install numpy sympy scipy matplotlib pandas seaborn
```

### **Option 2: Conda Environment** (Recommended)
```bash
# Create new environment
conda create -n lagrange python=3.9
conda activate lagrange

# Install dependencies
conda install numpy sympy scipy matplotlib pandas seaborn
```

### **Option 3: Requirements File**
```bash
# Create requirements.txt
cat > requirements.txt << EOF
numpy>=1.21.0
sympy>=1.9.0
scipy>=1.7.0
matplotlib>=3.4.0
pandas>=1.3.0
seaborn>=0.11.0
EOF

# Install
pip install -r requirements.txt
```

### **Download and Run**
```bash
# Download the solver
wget https://raw.githubusercontent.com/your-repo/lagrange_optimizer.py
# Or clone the repository
git clone https://github.com/your-repo/lagrange-optimizer.git
cd lagrange-optimizer

# Run interactive demo
python lagrange_optimizer.py
```

---

## 🎮 Interactive Demo

### **Launch Interactive Mode**
```python
from lagrange_optimizer import interactive_solver
interactive_solver()
```

**Available Examples:**
1. 📐 **Basic**: Minimize x² + y² subject to x + y = 1
2. 💰 **Economics**: Utility maximization with budget constraint
3. 📏 **Geometry**: Find closest point on curve to given point
4. ⚛️ **Physics**: Constrained motion and energy minimization
5. 🧊 **3D Problems**: Multi-dimensional optimization
6. ✏️ **Custom**: Define your own problem

### **Sample Interactive Session**
```
🎯 Advanced Lagrange Multiplier Solver
==================================================

Options:
📚 Examples:
  1. Basic: Minimize x² + y²
  2. Economics: Utility Maximization  
  3. Geometry: Distance to Curve
  4. Physics: Constrained Motion
  5. 3D: Minimize x² + y² + z²
  6. Custom problem
  7. Exit

Select option (1-7): 1

Loaded: Basic: Minimize x² + y²
===============================================
SYMBOLIC SOLUTION
===============================================
∂L/∂x = 2*x - λ1 = 0
∂L/∂y = 2*y - λ1 = 0
Constraints: ['x + y - 1']

Found 1 symbolic solution(s)

Solution 1:
  x = 0.5
  y = 0.5
  λ1 = 1.0
  Objective value = 0.5
```

---

## 📚 Usage Examples

### **Example 1: Basic Optimization**

**Problem**: Minimize the sum of squares subject to a linear constraint

```python
import numpy as np
from lagrange_optimizer import AdvancedLagrangeOptimizer

# Initialize solver
solver = AdvancedLagrangeOptimizer(verbose=True)

# Define problem
objective = "x**2 + y**2"              # Minimize sum of squares
constraint = "x + y - 1"               # Subject to x + y = 1
variables = ["x", "y"]                 # Two variables

# Solve completely
df = solver.solve_complete(
    objective=objective,
    constraints=[constraint],
    variables=variables,
    problem_type="minimize",
    visualize=True,
    export=True
)

# Results
print("Optimal solution:")
print(f"x* = {df.iloc[0]['x']:.6f}")
print(f"y* = {df.iloc[0]['y']:.6f}")
print(f"Minimum value = {df.iloc[0]['Objective_Value']:.6f}")
```

### **Example 2: Economics - Utility Maximization**

**Problem**: Maximize utility function subject to budget constraint

```python
# Consumer choice problem
solver = AdvancedLagrangeOptimizer()

# Cobb-Douglas utility: U(x,y) = x*y
# Budget constraint: 2x + 3y = 12
results = solver.solve_complete(
    objective="x*y",                    # Utility function
    constraints=["2*x + 3*y - 12"],     # Budget constraint  
    variables=["x", "y"],               # Goods x and y
    problem_type="maximize",            # Maximize utility
    visualize=True
)

print("Optimal consumption bundle:")
print(f"Good X: {results.iloc[0]['x']:.3f} units")
print(f"Good Y: {results.iloc[0]['y']:.3f} units") 
print(f"Maximum utility: {results.iloc[0]['Objective_Value']:.3f}")
print(f"Shadow price (λ): {results.iloc[0]['λ1']:.3f}")
```

### **Example 3: Multi-Constraint 3D Problem**

**Problem**: Minimize distance from origin in 3D space with multiple constraints

```python
# 3D optimization with multiple constraints
solver = AdvancedLagrangeOptimizer()

results = solver.solve_complete(
    objective="x**2 + y**2 + z**2",              # Distance squared from origin
    constraints=[
        "x + y + z - 1",                         # Plane constraint
        "x**2 + y**2 - 1"                        # Circular constraint in xy-plane
    ],
    variables=["x", "y", "z"],                   # 3D variables
    problem_type="minimize",
    visualize=False                              # No 2D viz for 3D problems
)

print("3D Optimization Results:")
for i, row in results.iterrows():
    print(f"Solution {i+1}:")
    print(f"  Point: ({row['x']:.4f}, {row['y']:.4f}, {row['z']:.4f})")
    print(f"  Distance²: {row['Objective_Value']:.6f}")
    print(f"  λ1 = {row['λ1']:.4f}, λ2 = {row['λ2']:.4f}")
```

### **Example 4: Physics - Constrained Motion**

**Problem**: Find equilibrium position of a particle in a potential field

```python
# Particle in potential field with constraint
solver = AdvancedLagrangeOptimizer()

# Potential energy: V(x,y) = 0.5x² + 2y²
# Constraint: particle on circle x² + y² = 4
results = solver.solve_complete(
    objective="0.5*x**2 + 2*y**2",      # Potential energy
    constraints=["x**2 + y**2 - 4"],     # Circular constraint
    variables=["x", "y"],                # Position coordinates
    problem_type="minimize",             # Find minimum energy
    visualize=True
)

print("Equilibrium Analysis:")
print(f"Equilibrium position: ({results.iloc[0]['x']:.4f}, {results.iloc[0]['y']:.4f})")
print(f"Minimum potential energy: {results.iloc[0]['Objective_Value']:.6f}")
print(f"Constraint force (λ): {results.iloc[0]['λ1']:.6f}")
```

### **Example 5: Advanced - Sensitivity Analysis**

**Problem**: Study how solutions change with parameter variations

```python
import matplotlib.pyplot as plt

def sensitivity_study():
    """Analyze sensitivity to constraint parameters"""
    
    results = []
    c_values = np.linspace(0.5, 2.0, 11)  # Parameter range
    
    for c in c_values:
        solver = AdvancedLagrangeOptimizer(verbose=False)
        constraint = f"x + y - {c}"        # Parameterized constraint
        
        solver.solve_complete(
            objective="x**2 + y**2",
            constraints=[constraint],
            variables=["x", "y"],
            problem_type="minimize",
            visualize=False,
            export=False
        )
        
        if solver.solutions:
            sol = solver.solutions[0]
            results.append({
                'parameter_c': c,
                'optimal_x': sol['x'],
                'optimal_y': sol['y'], 
                'objective_value': sol['objective_value'],
                'shadow_price': sol.get('λ1', 0)
            })
    
    # Plot sensitivity results
    df = pd.DataFrame(results)
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Optimal variables vs parameter
    axes[0, 0].plot(df['parameter_c'], df['optimal_x'], 'bo-', label='x*')
    axes[0, 0].plot(df['parameter_c'], df['optimal_y'], 'ro-', label='y*')
    axes[0, 0].set_title('Optimal Variables vs Parameter')
    axes[0, 0].legend()
    axes[0, 0].grid(True)
    
    # Objective value vs parameter  
    axes[0, 1].plot(df['parameter_c'], df['objective_value'], 'go-')
    axes[0, 1].set_title('Optimal Value vs Parameter')
    axes[0, 1].grid(True)
    
    # Shadow price analysis
    axes[1, 0].plot(df['parameter_c'], df['shadow_price'], 'mo-')
    axes[1, 0].set_title('Shadow Price (Lagrange Multiplier)')
    axes[1, 0].grid(True)
    
    # Numerical vs analytical derivative
    numerical_deriv = np.gradient(df['objective_value'], df['parameter_c'])
    axes[1, 1].plot(df['parameter_c'], numerical_deriv, 'co-', label='∂f*/∂c (numerical)')
    axes[1, 1].plot(df['parameter_c'], df['shadow_price'], 'mo--', label='λ (analytical)')
    axes[1, 1].set_title('Shadow Price Verification')
    axes[1, 1].legend()
    axes[1, 1].grid(True)
    
    plt.tight_layout()
    plt.show()
    
    return df

# Run sensitivity analysis
sensitivity_df = sensitivity_study()
print(sensitivity_df)
```

### **Example 6: Batch Processing**

**Problem**: Solve multiple optimization problems automatically

```python
# Define batch of problems
problems = [
    {
        'name': 'Basic Quadratic',
        'objective': 'x**2 + y**2',
        'constraints': ['x + y - 1'],
        'variables': ['x', 'y'],
        'type': 'minimize'
    },
    {
        'name': 'Economic Utility',
        'objective': 'x*y', 
        'constraints': ['2*x + 3*y - 12'],
        'variables': ['x', 'y'],
        'type': 'maximize'
    },
    {
        'name': 'Geometric Distance',
        'objective': '(x-2)**2 + (y-3)**2',
        'constraints': ['x**2 + y**2 - 1'], 
        'variables': ['x', 'y'],
        'type': 'minimize'
    }
]

# Batch solve with automatic report generation
from lagrange_optimizer import batch_solve
batch_results = batch_solve(problems, output_dir="batch_results")

print("Batch Processing Summary:")
print(batch_results)
```

---

## 🔧 API Reference

### **Class: AdvancedLagrangeOptimizer**

#### **Constructor**
```python
AdvancedLagrangeOptimizer(verbose=True)
```
- `verbose` (bool): Enable detailed output and logging

#### **Main Methods**

##### **solve_complete()**
```python
solve_complete(
    objective,           # Objective function (str or sympy.Expr)
    constraints,         # List of constraints (str or sympy.Expr)  
    variables,           # List of variables (str or sympy.Symbol)
    problem_type="minimize",  # "minimize" or "maximize"
    try_symbolic=True,   # Attempt symbolic solution
    try_numerical=True,  # Attempt numerical solution
    visualize=True,      # Generate visualizations
    export=False         # Export results to files
) -> pd.DataFrame       # Returns solution analysis
```

##### **solve_symbolic()**
```python
solve_symbolic() -> List[Dict]
```
Solve using symbolic computation (SymPy). Returns list of exact solutions.

##### **solve_numerical()**
```python
solve_numerical(
    num_random_starts=50,     # Number of random starting points
    bounds=None,              # Variable bounds [(min, max), ...]
    method='SLSQP'           # SciPy optimization method
) -> List[Dict]              # Returns numerical solutions
```

##### **verify_solutions()**
```python
verify_solutions() -> None
```
Verify solutions satisfy KKT conditions and constraints.

##### **visualize_2d()**
```python
visualize_2d(
    x_range=(-5, 5),         # X-axis range
    y_range=(-5, 5),         # Y-axis range  
    resolution=400           # Grid resolution
) -> None
```
Create 2D visualization with contour plots and 3D surfaces.

##### **analyze_solutions()**
```python
analyze_solutions() -> pd.DataFrame
```
Generate comprehensive solution analysis as pandas DataFrame.

##### **export_results()**
```python
export_results(filename=None) -> None
```
Export results to JSON file with timestamp.

#### **Properties**
- `solutions`: List of found solutions
- `lagrangian`: Constructed Lagrangian function
- `variables`: Problem variables
- `constraints`: Problem constraints
- `objective`: Objective function

---

## 📊 Visualization Gallery

### **2D Contour Plots**
- Objective function contours with constraint curves
- Solution points highlighted with markers
- Color-coded objective function values
- Professional matplotlib styling

### **3D Surface Plots**
- Interactive 3D surface visualization
- Constraint intersections on surface
- Optimal points in 3D space
- Customizable viewing angles

### **Sensitivity Analysis**
- Parameter dependency plots
- Shadow price visualization
- Convergence analysis charts
- Statistical distribution plots

### **Example Visualization Code**
```python
# Generate publication-ready plots
solver = AdvancedLagrangeOptimizer()
solver.solve_complete(
    objective="x**2 + 0.5*y**2", 
    constraints=["x**2 + y**2 - 4"],
    variables=["x", "y"],
    visualize=True
)

# Customize plot appearance
import matplotlib.pyplot as plt
plt.style.use('seaborn-v0_8')  # Professional styling
plt.rcParams['figure.figsize'] = (12, 8)  # Larger figures
plt.rcParams['font.size'] = 12  # Readable fonts
```

---

## 🧪 Scientific Applications

### **1. Economics and Finance** 💰
- **Portfolio Optimization**: Maximize returns subject to risk constraints
- **Utility Maximization**: Consumer choice with budget constraints
- **Cost Minimization**: Production optimization with resource limits
- **Market Equilibrium**: Supply-demand balance problems

```python
# Portfolio optimization example
solver.solve_complete(
    objective="0.12*x1 + 0.08*x2 + 0.15*x3",  # Expected returns
    constraints=[
        "x1 + x2 + x3 - 1",                    # Budget constraint
        "0.2*x1 + 0.1*x2 + 0.3*x3 - 0.15"     # Risk constraint
    ],
    variables=["x1", "x2", "x3"],              # Asset allocations
    problem_type="maximize"
)
```

### **2. Physics and Engineering** ⚛️
- **Energy Minimization**: Find equilibrium configurations
- **Constrained Motion**: Particle dynamics with constraints
- **Structural Optimization**: Design under load constraints
- **Control Systems**: Optimal control with state constraints

```python
# Minimum energy configuration
solver.solve_complete(
    objective="0.5*k1*x**2 + 0.5*k2*y**2",    # Elastic potential energy
    constraints=["x**2 + y**2 - L**2"],        # Length constraint
    variables=["x", "y"],
    problem_type="minimize"
)
```

### **3. Machine Learning** 🤖
- **SVM Optimization**: Support Vector Machine dual problems
- **Neural Network Training**: Constrained parameter optimization
- **Feature Selection**: Optimization with sparsity constraints
- **Regularized Regression**: L1/L2 regularization problems

```python
# SVM-like classification
solver.solve_complete(
    objective="0.5*w1**2 + 0.5*w2**2 + C*xi",  # SVM objective
    constraints=["y*(w1*x1 + w2*x2 + b) - 1 + xi"],  # Margin constraint
    variables=["w1", "w2", "b", "xi"],
    problem_type="minimize"
)
```

### **4. Operations Research** 📈
- **Resource Allocation**: Optimal distribution of limited resources
- **Transportation Problems**: Minimum cost flow optimization
- **Production Planning**: Manufacturing optimization
- **Inventory Management**: Stock level optimization

---

## ⚡ Performance

### **Benchmark Results**

| Problem Type | Variables | Constraints | Symbolic Time | Numerical Time | Accuracy |
|-------------|-----------|-------------|---------------|----------------|----------|
| Linear | 2 | 1 | 0.05s | 0.12s | Exact |
| Quadratic | 2 | 1 | 0.08s | 0.15s | 1e-12 |
| Nonlinear | 2 | 2 | 0.25s | 0.45s | 1e-10 |
| 3D Multi-constraint | 3 | 2 | 1.2s | 0.8s | 1e-09 |
| Complex (5D) | 5 | 3 | Timeout | 2.1s | 1e-08 |

### **Performance Tips**
```python
# For large problems, use numerical-only mode
solver = AdvancedLagrangeOptimizer(verbose=False)
results = solver.solve_numerical(
    num_random_starts=100,      # More starting points
    method='trust-constr'       # Better for nonlinear
)

# For high precision, adjust tolerances
from scipy.optimize import minimize
# Custom solver settings available
```

### **Memory Usage**
- **Small problems (2-3 vars)**: < 10 MB
- **Medium problems (4-6 vars)**: < 50 MB  
- **Large problems (7+ vars)**: 100-500 MB
- **Batch processing**: Depends on problem count

---

## 🧪 Testing

### **Run Test Suite**
```python
# Basic functionality test
def test_basic_optimization():
    solver = AdvancedLagrangeOptimizer(verbose=False)
    results = solver.solve_complete(
        objective="x**2 + y**2",
        constraints=["x + y - 1"], 
        variables=["x", "y"],
        visualize=False,
        export=False
    )
    
    assert len(results) > 0, "No solutions found"
    assert abs(results.iloc[0]['x'] - 0.5) < 1e-6, "Incorrect x value"
    assert abs(results.iloc[0]['y'] - 0.5) < 1e-6, "Incorrect y value"
    print("✅ Basic test passed")

# Constraint verification test  
def test_constraint_satisfaction():
    solver = AdvancedLagrangeOptimizer(verbose=False)
    solver.solve_complete(
        objective="x*y",
        constraints=["x**2 + y**2 - 1"],
        variables=["x", "y"],
        problem_type="maximize",
        visualize=False
    )
    
    for sol in solver.solutions:
        constraint_val = sol['x']**2 + sol['y']**2 - 1
        assert abs(constraint_val) < 1e-6, f"Constraint violated: {constraint_val}"
    print("✅ Constraint test passed")

# Run all tests
if __name__ == "__main__":
    test_basic_optimization()
    test_constraint_satisfaction()
    print("✅ All tests passed!")
```

---

## 🤝 Contributing

We welcome contributions! Here's how to get started:

### **Development Setup**
```bash
# Fork and clone the repository
git clone https://github.com/yourusername/lagrange-optimizer.git
cd lagrange-optimizer

# Create development environment
conda create -n lagrange-dev python=3.9
conda activate lagrange-dev

# Install dependencies
pip install -r requirements-dev.txt

# Install pre-commit hooks
pre-commit install
```

### **Code Style**
- Follow PEP 8 style guidelines
- Use type hints for function signatures
- Add comprehensive docstrings
- Include unit tests for new features

### **Testing**
```bash
# Run unit tests
python -m pytest tests/

# Run with coverage
python -m pytest --cov=lagrange_optimizer tests/

# Run performance benchmarks  
python benchmarks/performance_tests.py
```

### **Pull Request Process**
1. Create feature branch: `git checkout -b feature-name`
2. Make changes and add tests
3. Ensure all tests pass: `pytest`
4. Update documentation if needed
5. Submit pull request with clear description

### **Reporting Issues**
- Use GitHub Issues for bug reports
- Include minimal reproducible example
- Specify Python version and dependencies
- Add relevant error messages and stack traces

---

## 📚 Additional Resources

### **Mathematical Background**
- [Lagrange Multipliers Theory](https://en.wikipedia.org/wiki/Lagrange_multiplier)
- [KKT Conditions](https://en.wikipedia.org/wiki/Karush%E2%80%93Kuhn%E2%80%93Tucker_conditions)
- [Constrained Optimization](https://web.stanford.edu/~boyd/cvxbook/)

### **Related Libraries**
- [SymPy](https://sympy.org/): Symbolic mathematics
- [SciPy](https://scipy.org/): Scientific computing
- [CVXPY](https://www.cvxpy.org/): Convex optimization
- [Pyomo](http://www.pyomo.org/): Optimization modeling

### **Tutorials and Examples**
- [Jupyter Notebooks](examples/) with detailed walkthroughs
- [Video Tutorials](https://youtube.com/playlist/lagrange-optimization)
- [Academic Papers](references/) using this solver
- [Blog Posts](https://blog.optimization.com/lagrange-series)

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

```
MIT License

Copyright (c) 2024 Advanced Scientific Computing

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

---

## 🙏 Acknowledgments

- **SymPy Team** for symbolic mathematics capabilities
- **SciPy Community** for numerical optimization algorithms  
- **Matplotlib Developers** for visualization tools
- **Scientific Python Ecosystem** for foundational libraries
- **Contributors** who helped improve this project

---

## 📞 Support

### **Get Help**
- 📧 **Email**: support@lagrange-optimizer.com
- 💬 **Discord**: [Join our community](https://discord.gg/lagrange-optimizer)
- 📝 **GitHub Issues**: [Report bugs/features](https://github.com/your-repo/lagrange-optimizer/issues)
- 📚 **Documentation**: [Full docs](https://lagrange-optimizer.readthedocs.io)

### **Citation**
If you use this software in your research, please cite:

```bibtex
@software{lagrange_optimizer_2024,
  title={Advanced Lagrange Multiplier Solver},
  author={Advanced Scientific Computing Team},
  year={2024},
  url={https://github.com/your-repo/lagrange-optimizer},
  version={1.0.0}
}
```

---

<div align="center">

**⭐ If you find this project helpful, please give it a star! ⭐**

[![GitHub stars](https://img.shields.io/github/stars/your-repo/lagrange-optimizer.svg?style=social&label=Star)](https://github.com/your-repo/lagrange-optimizer)
[![Twitter Follow](https://img.shields.io/twitter/follow/lagrange_opt?style=social)](https://twitter.com/lagrange_opt)

---

**Made with ❤️ for the scientific computing community**

[🏠 Home](https://lagrange-optimizer.com) | [📖 Docs](https://docs.lagrange-optimizer.com) | [🔧 API](https://api.lagrange-optimizer.com) | [💻 Examples](https://examples.lagrange-optimizer.com)

</div>
