"""
Advanced Lagrange Multiplier Solver for Scientific Computing
===========================================================

A comprehensive Python implementation for solving constrained optimization problems
using the method of Lagrange multipliers with advanced numerical techniques,
visualization, and analysis capabilities.

Features:
- Multiple constraint handling
- Numerical and symbolic solutions
- Interactive visualization
- Solution verification and analysis
- Export capabilities
- Comprehensive error handling

"""

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from scipy.optimize import minimize, fsolve
import warnings
import json
from datetime import datetime
import seaborn as sns
from typing import List, Dict, Union, Tuple, Optional
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# Set plotting style
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")
warnings.filterwarnings('ignore', category=RuntimeWarning)

class AdvancedLagrangeOptimizer:
    """
    Advanced Lagrange multiplier solver with comprehensive features for scientific computing.
    """
    
    def __init__(self, verbose: bool = True):
        """
        Initialize the optimizer.
        
        Args:
            verbose (bool): Enable detailed output
        """
        self.verbose = verbose
        self.solutions = []
        self.lagrangian = None
        self.variables = None
        self.multipliers = None
        self.objective = None
        self.constraints = []
        self.problem_type = 'minimize'
        
    def setup_problem(self, 
                     objective: Union[str, sp.Expr],
                     constraints: Union[List[str], List[sp.Expr]],
                     variables: Union[List[str], List[sp.Symbol]],
                     problem_type: str = 'minimize') -> None:
        """
        Set up the optimization problem.
        
        Args:
            objective: Objective function as string or SymPy expression
            constraints: List of constraint equations (set to zero)
            variables: List of variables as strings or SymPy symbols
            problem_type: 'minimize' or 'maximize'
        """
        self.problem_type = problem_type.lower()
        
        # Convert to SymPy expressions if needed
        if isinstance(variables[0], str):
            self.variables = [sp.Symbol(var) for var in variables]
        else:
            self.variables = variables
            
        if isinstance(objective, str):
            self.objective = sp.sympify(objective)
        else:
            self.objective = objective
            
        self.constraints = []
        for constraint in constraints:
            if isinstance(constraint, str):
                self.constraints.append(sp.sympify(constraint))
            else:
                self.constraints.append(constraint)
        
        # Create Lagrange multipliers
        self.multipliers = [sp.Symbol(f'λ{i+1}') for i in range(len(self.constraints))]
        
        # Construct Lagrangian
        self._construct_lagrangian()
        
        if self.verbose:
            print(f"Problem setup complete:")
            print(f"Type: {self.problem_type.capitalize()}")
            print(f"Variables: {[str(v) for v in self.variables]}")
            print(f"Objective: {self.objective}")
            print(f"Constraints: {[str(c) for c in self.constraints]}")
    
    def _construct_lagrangian(self) -> None:
        """Construct the Lagrangian function."""
        lagrangian = self.objective.copy()
        
        for i, constraint in enumerate(self.constraints):
            lagrangian -= self.multipliers[i] * constraint
            
        self.lagrangian = lagrangian
        
        if self.verbose:
            print(f"\nLagrangian: {self.lagrangian}")
    
    def solve_symbolic(self) -> List[Dict]:
        """
        Solve the optimization problem symbolically using SymPy.
        
        Returns:
            List of solution dictionaries
        """
        if self.verbose:
            print("\n" + "="*60)
            print("SYMBOLIC SOLUTION")
            print("="*60)
        
        try:
            # Compute gradients
            gradient_eqs = []
            
            # ∇L = 0 with respect to all variables
            for var in self.variables:
                grad_eq = sp.diff(self.lagrangian, var)
                gradient_eqs.append(grad_eq)
                if self.verbose:
                    print(f"∂L/∂{var} = {grad_eq} = 0")
            
            # Add constraint equations
            gradient_eqs.extend(self.constraints)
            if self.verbose:
                print(f"Constraints: {[str(c) for c in self.constraints]}")
            
            # Solve the system
            all_vars = self.variables + self.multipliers
            solutions = sp.solve(gradient_eqs, all_vars, dict=True)
            
            if not solutions:
                if self.verbose:
                    print("No symbolic solutions found. Trying numerical methods...")
                return []
            
            # Process solutions
            processed_solutions = []
            for sol in solutions:
                processed_sol = self._process_solution(sol)
                if processed_sol:
                    processed_solutions.append(processed_sol)
            
            self.solutions = processed_solutions
            
            if self.verbose:
                print(f"\nFound {len(processed_solutions)} symbolic solution(s)")
                for i, sol in enumerate(processed_solutions):
                    print(f"\nSolution {i+1}:")
                    for key, value in sol.items():
                        if key != 'objective_value' and key != 'constraint_values':
                            print(f"  {key} = {value}")
                    print(f"  Objective value = {sol['objective_value']}")
            
            return processed_solutions
            
        except Exception as e:
            logger.error(f"Symbolic solution failed: {e}")
            if self.verbose:
                print(f"Symbolic solution failed: {e}")
            return []
    
    def solve_numerical(self, 
                       num_random_starts: int = 50,
                       bounds: Optional[List[Tuple]] = None,
                       method: str = 'SLSQP') -> List[Dict]:
        """
        Solve using numerical optimization methods.
        
        Args:
            num_random_starts: Number of random starting points
            bounds: Variable bounds as list of tuples
            method: Optimization method ('SLSQP', 'trust-constr', etc.)
        
        Returns:
            List of solution dictionaries
        """
        if self.verbose:
            print("\n" + "="*60)
            print("NUMERICAL SOLUTION")
            print("="*60)
        
        # Convert to numerical functions
        objective_func = sp.lambdify(self.variables, self.objective, 'numpy')
        constraint_funcs = [sp.lambdify(self.variables, c, 'numpy') for c in self.constraints]
        
        # Set default bounds if not provided
        if bounds is None:
            bounds = [(-10, 10)] * len(self.variables)
        
        solutions = []
        
        # Multiple random starting points
        np.random.seed(42)  # For reproducibility
        
        for i in range(num_random_starts):
            # Generate random starting point
            x0 = np.array([np.random.uniform(b[0], b[1]) for b in bounds])
            
            try:
                # Define constraints for scipy
                constraints = []
                for j, cf in enumerate(constraint_funcs):
                    constraints.append({
                        'type': 'eq',
                        'fun': lambda x, func=cf: func(*x)
                    })
                
                # Objective function (negate for maximization)
                if self.problem_type == 'maximize':
                    obj_func = lambda x: -objective_func(*x)
                else:
                    obj_func = lambda x: objective_func(*x)
                
                # Solve
                result = minimize(obj_func, x0, method=method, 
                                constraints=constraints, bounds=bounds,
                                options={'ftol': 1e-9, 'disp': False})
                
                if result.success:
                    # Check if solution satisfies constraints
                    constraint_violations = [abs(cf(*result.x)) for cf in constraint_funcs]
                    if all(cv < 1e-6 for cv in constraint_violations):
                        
                        # Create solution dictionary
                        sol_dict = {}
                        for j, var in enumerate(self.variables):
                            sol_dict[str(var)] = float(result.x[j])
                        
                        # Calculate objective value (correct sign)
                        obj_val = float(objective_func(*result.x))
                        sol_dict['objective_value'] = obj_val
                        
                        # Estimate Lagrange multipliers
                        if hasattr(result, 'multipliers') and result.multipliers is not None:
                            for j, mult in enumerate(self.multipliers):
                                if j < len(result.multipliers):
                                    sol_dict[str(mult)] = float(result.multipliers[j])
                        else:
                            # Estimate multipliers using gradients
                            multipliers = self._estimate_multipliers(result.x)
                            for j, mult in enumerate(self.multipliers):
                                if j < len(multipliers):
                                    sol_dict[str(mult)] = float(multipliers[j])
                        
                        sol_dict['constraint_values'] = constraint_violations
                        solutions.append(sol_dict)
                
            except Exception as e:
                logger.debug(f"Numerical optimization attempt {i+1} failed: {e}")
                continue
        
        # Remove duplicates
        unique_solutions = self._remove_duplicate_solutions(solutions)
        
        # Sort by objective value
        unique_solutions.sort(key=lambda x: x['objective_value'])
        
        self.solutions.extend(unique_solutions)
        
        if self.verbose:
            print(f"\nFound {len(unique_solutions)} numerical solution(s)")
            for i, sol in enumerate(unique_solutions[:5]):  # Show first 5
                print(f"\nSolution {i+1}:")
                for key, value in sol.items():
                    if key not in ['constraint_values']:
                        if isinstance(value, float):
                            print(f"  {key} = {value:.6f}")
                        else:
                            print(f"  {key} = {value}")
        
        return unique_solutions
    
    def _estimate_multipliers(self, x: np.ndarray) -> List[float]:
        """Estimate Lagrange multipliers using gradient information."""
        try:
            # Compute gradients numerically
            obj_grad = self._numerical_gradient(self.objective, x)
            constraint_grads = [self._numerical_gradient(c, x) for c in self.constraints]
            
            # Solve: ∇f + Σ(λᵢ∇gᵢ) = 0 for λᵢ
            # This is a least squares problem: A*λ = -∇f where A = [∇g₁, ∇g₂, ...]
            A = np.column_stack(constraint_grads)
            b = -obj_grad
            
            # Solve using least squares
            multipliers, residuals, rank, s = np.linalg.lstsq(A.T, b, rcond=None)
            
            return multipliers.tolist()
            
        except Exception as e:
            logger.debug(f"Multiplier estimation failed: {e}")
            return [0.0] * len(self.constraints)
    
    def _numerical_gradient(self, expr: sp.Expr, x: np.ndarray, h: float = 1e-8) -> np.ndarray:
        """Compute numerical gradient of expression at point x."""
        func = sp.lambdify(self.variables, expr, 'numpy')
        grad = np.zeros_like(x)
        
        for i in range(len(x)):
            x_plus = x.copy()
            x_minus = x.copy()
            x_plus[i] += h
            x_minus[i] -= h
            
            grad[i] = (func(*x_plus) - func(*x_minus)) / (2 * h)
        
        return grad
    
    def _process_solution(self, sol: Dict) -> Optional[Dict]:
        """Process and validate a solution."""
        try:
            processed = {}
            
            # Extract variable values
            for var in self.variables:
                if var in sol:
                    value = complex(sol[var])
                    if abs(value.imag) < 1e-10:  # Essentially real
                        processed[str(var)] = float(value.real)
                    else:
                        return None  # Skip complex solutions
                else:
                    return None
            
            # Extract multiplier values
            for mult in self.multipliers:
                if mult in sol:
                    value = complex(sol[mult])
                    if abs(value.imag) < 1e-10:
                        processed[str(mult)] = float(value.real)
                    else:
                        processed[str(mult)] = 0.0
                else:
                    processed[str(mult)] = 0.0
            
            # Calculate objective value
            var_values = [processed[str(var)] for var in self.variables]
            obj_func = sp.lambdify(self.variables, self.objective, 'numpy')
            processed['objective_value'] = float(obj_func(*var_values))
            
            # Verify constraints
            constraint_values = []
            for constraint in self.constraints:
                const_func = sp.lambdify(self.variables, constraint, 'numpy')
                constraint_values.append(float(const_func(*var_values)))
            
            processed['constraint_values'] = constraint_values
            
            return processed
            
        except Exception as e:
            logger.debug(f"Solution processing failed: {e}")
            return None
    
    def _remove_duplicate_solutions(self, solutions: List[Dict], tolerance: float = 1e-6) -> List[Dict]:
        """Remove duplicate solutions based on variable values."""
        unique_solutions = []
        
        for sol in solutions:
            is_duplicate = False
            for unique_sol in unique_solutions:
                # Compare variable values
                differences = []
                for var in self.variables:
                    var_str = str(var)
                    if var_str in sol and var_str in unique_sol:
                        differences.append(abs(sol[var_str] - unique_sol[var_str]))
                
                if differences and max(differences) < tolerance:
                    is_duplicate = True
                    break
            
            if not is_duplicate:
                unique_solutions.append(sol)
        
        return unique_solutions
    
    def verify_solutions(self) -> None:
        """Verify that solutions satisfy the KKT conditions."""
        if not self.solutions:
            print("No solutions to verify.")
            return
        
        if self.verbose:
            print("\n" + "="*60)
            print("SOLUTION VERIFICATION")
            print("="*60)
        
        for i, sol in enumerate(self.solutions):
            if self.verbose:
                print(f"\nVerifying Solution {i+1}:")
            
            # Check constraint satisfaction
            constraint_satisfied = True
            for j, constraint_val in enumerate(sol.get('constraint_values', [])):
                satisfied = abs(constraint_val) < 1e-6
                constraint_satisfied = constraint_satisfied and satisfied
                if self.verbose:
                    status = "✓" if satisfied else "✗"
                    print(f"  Constraint {j+1}: {constraint_val:.2e} {status}")
            
            # Check gradient condition (∇L = 0)
            var_values = [sol[str(var)] for var in self.variables]
            mult_values = [sol.get(str(mult), 0) for mult in self.multipliers]
            
            gradient_satisfied = True
            for var in self.variables:
                grad_val = sp.diff(self.lagrangian, var)
                grad_func = sp.lambdify(self.variables + self.multipliers, grad_val, 'numpy')
                grad_at_sol = grad_func(*(var_values + mult_values))
                
                satisfied = abs(grad_at_sol) < 1e-6
                gradient_satisfied = gradient_satisfied and satisfied
                if self.verbose:
                    status = "✓" if satisfied else "✗"
                    print(f"  ∂L/∂{var} = {grad_at_sol:.2e} {status}")
            
            overall_status = "VALID" if (constraint_satisfied and gradient_satisfied) else "INVALID"
            if self.verbose:
                print(f"  Overall: {overall_status}")
    
    def visualize_2d(self, 
                    x_range: Tuple[float, float] = (-5, 5),
                    y_range: Tuple[float, float] = (-5, 5),
                    resolution: int = 400) -> None:
        """
        Create 2D visualization for 2-variable problems.
        
        Args:
            x_range: Range for x-axis
            y_range: Range for y-axis  
            resolution: Grid resolution
        """
        if len(self.variables) != 2:
            print("2D visualization only available for 2-variable problems.")
            return
        
        x_var, y_var = self.variables
        
        # Create meshgrid
        x = np.linspace(x_range[0], x_range[1], resolution)
        y = np.linspace(y_range[0], y_range[1], resolution)
        X, Y = np.meshgrid(x, y)
        
        # Evaluate objective function
        obj_func = sp.lambdify([x_var, y_var], self.objective, 'numpy')
        Z = obj_func(X, Y)
        
        # Create subplots
        fig, axes = plt.subplots(1, 2, figsize=(15, 6))
        
        # Contour plot
        ax1 = axes[0]
        contour = ax1.contour(X, Y, Z, levels=20, alpha=0.6)
        ax1.clabel(contour, inline=True, fontsize=8)
        contour_fill = ax1.contourf(X, Y, Z, levels=50, alpha=0.3, cmap='viridis')
        plt.colorbar(contour_fill, ax=ax1, label='Objective Function')
        
        # Plot constraints
        for i, constraint in enumerate(self.constraints):
            const_func = sp.lambdify([x_var, y_var], constraint, 'numpy')
            C = const_func(X, Y)
            ax1.contour(X, Y, C, levels=[0], colors=f'C{i+1}', linewidths=3, 
                       linestyles='--', alpha=0.8)
        
        # Plot solutions
        if self.solutions:
            sol_x = [sol[str(x_var)] for sol in self.solutions if str(x_var) in sol]
            sol_y = [sol[str(y_var)] for sol in self.solutions if str(y_var) in sol]
            ax1.scatter(sol_x, sol_y, color='red', s=100, marker='*', 
                       edgecolors='black', linewidth=2, label='Solutions', zorder=5)
        
        ax1.set_xlabel(f'{x_var}')
        ax1.set_ylabel(f'{y_var}')
        ax1.set_title('Contour Plot with Constraints')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        
        # 3D surface plot
        ax2 = fig.add_subplot(122, projection='3d')
        surface = ax2.plot_surface(X, Y, Z, cmap='viridis', alpha=0.7)
        
        # Plot constraint curves on surface
        for i, constraint in enumerate(self.constraints):
            const_func = sp.lambdify([x_var, y_var], constraint, 'numpy')
            C = const_func(X, Y)
            ax2.contour(X, Y, C, levels=[0], colors=f'C{i+1}', linewidths=3,
                       offset=np.min(Z), zdir='z', alpha=0.8)
        
        # Plot solutions on surface
        if self.solutions:
            sol_x = [sol[str(x_var)] for sol in self.solutions if str(x_var) in sol]
            sol_y = [sol[str(y_var)] for sol in self.solutions if str(y_var) in sol]
            sol_z = [sol['objective_value'] for sol in self.solutions]
            ax2.scatter(sol_x, sol_y, sol_z, color='red', s=100, marker='*',
                       edgecolors='black', linewidth=2, label='Solutions')
        
        ax2.set_xlabel(f'{x_var}')
        ax2.set_ylabel(f'{y_var}')
        ax2.set_zlabel('Objective Function')
        ax2.set_title('3D Surface Plot')
        
        plt.colorbar(surface, ax=ax2, shrink=0.5, aspect=5, label='Objective Function')
        plt.tight_layout()
        plt.show()
    
    def analyze_solutions(self) -> pd.DataFrame:
        """Create detailed analysis of solutions."""
        if not self.solutions:
            print("No solutions to analyze.")
            return pd.DataFrame()
        
        # Create DataFrame
        df_data = []
        for i, sol in enumerate(self.solutions):
            row = {'Solution': i + 1}
            
            # Variable values
            for var in self.variables:
                row[str(var)] = sol.get(str(var), np.nan)
            
            # Multiplier values
            for mult in self.multipliers:
                row[str(mult)] = sol.get(str(mult), np.nan)
            
            # Objective value
            row['Objective_Value'] = sol.get('objective_value', np.nan)
            
            # Constraint violations
            for j, const_val in enumerate(sol.get('constraint_values', [])):
                row[f'Constraint_{j+1}_Violation'] = abs(const_val)
            
            df_data.append(row)
        
        df = pd.DataFrame(df_data)
        
        if self.verbose:
            print("\n" + "="*60)
            print("SOLUTION ANALYSIS")
            print("="*60)
            print(df.round(6))
            
            # Summary statistics
            print("\nSummary Statistics:")
            print(f"Number of solutions: {len(self.solutions)}")
            if len(self.solutions) > 0:
                obj_values = [sol['objective_value'] for sol in self.solutions]
                print(f"Best objective value: {min(obj_values) if self.problem_type == 'minimize' else max(obj_values):.6f}")
                print(f"Worst objective value: {max(obj_values) if self.problem_type == 'minimize' else min(obj_values):.6f}")
                print(f"Mean objective value: {np.mean(obj_values):.6f}")
                print(f"Std objective value: {np.std(obj_values):.6f}")
        
        return df
    
    def export_results(self, filename: Optional[str] = None) -> None:
        """Export results to JSON file."""
        if not filename:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"lagrange_results_{timestamp}.json"
        
        export_data = {
            'problem': {
                'type': self.problem_type,
                'variables': [str(var) for var in self.variables],
                'objective': str(self.objective),
                'constraints': [str(const) for const in self.constraints],
                'lagrangian': str(self.lagrangian)
            },
            'solutions': self.solutions,
            'timestamp': datetime.now().isoformat(),
            'num_solutions': len(self.solutions)
        }
        
        with open(filename, 'w') as f:
            json.dump(export_data, f, indent=2, default=str)
        
        if self.verbose:
            print(f"\nResults exported to {filename}")
    
    def solve_complete(self, 
                      objective: Union[str, sp.Expr],
                      constraints: Union[List[str], List[sp.Expr]], 
                      variables: Union[List[str], List[sp.Symbol]],
                      problem_type: str = 'minimize',
                      try_symbolic: bool = True,
                      try_numerical: bool = True,
                      visualize: bool = True,
                      export: bool = False) -> pd.DataFrame:
        """
        Complete solution pipeline.
        
        Args:
            objective: Objective function
            constraints: List of constraints
            variables: List of variables
            problem_type: 'minimize' or 'maximize'
            try_symbolic: Try symbolic solution first
            try_numerical: Try numerical solution
            visualize: Create visualization (for 2D problems)
            export: Export results to file
        
        Returns:
            DataFrame with solution analysis
        """
        # Setup problem
        self.setup_problem(objective, constraints, variables, problem_type)
        
        # Clear previous solutions
        self.solutions = []
        
        # Try symbolic solution first
        if try_symbolic:
            symbolic_solutions = self.solve_symbolic()
            if not symbolic_solutions and self.verbose:
                print("Symbolic solution unsuccessful, proceeding with numerical methods...")
        
        # Try numerical solution
        if try_numerical:
            self.solve_numerical()
        
        # Verify solutions
        self.verify_solutions()
        
        # Visualize if 2D
        if visualize and len(self.variables) == 2 and self.solutions:
            self.visualize_2d()
        
        # Analyze solutions
        df = self.analyze_solutions()
        
        # Export if requested
        if export:
            self.export_results()
        
        return df


def interactive_solver():
    """Interactive command-line interface for the optimizer."""
    print("🎯 Advanced Lagrange Multiplier Solver")
    print("=" * 50)
    
    optimizer = AdvancedLagrangeOptimizer(verbose=True)
    
    # Load predefined examples
    examples = {
        '1': {
            'name': 'Basic: Minimize x² + y²',
            'objective': 'x**2 + y**2',
            'constraints': ['x + y - 1'],
            'variables': ['x', 'y'],
            'type': 'minimize'
        },
        '2': {
            'name': 'Economics: Utility Maximization',
            'objective': 'x*y',
            'constraints': ['2*x + 3*y - 12'],
            'variables': ['x', 'y'],
            'type': 'maximize'
        },
        '3': {
            'name': 'Geometry: Distance to Curve',
            'objective': '(x-1)**2 + (y-2)**2',
            'constraints': ['x**2 + y**2 - 1'],
            'variables': ['x', 'y'],
            'type': 'minimize'
        },
        '4': {
            'name': 'Physics: Constrained Motion',
            'objective': '0.5*x**2 + 2*y**2',
            'constraints': ['x**2 + y**2 - 4'],
            'variables': ['x', 'y'],
            'type': 'minimize'
        },
        '5': {
            'name': '3D: Minimize x² + y² + z²',
            'objective': 'x**2 + y**2 + z**2',
            'constraints': ['x + y + z - 1', 'x**2 + y**2 - 1'],
            'variables': ['x', 'y', 'z'],
            'type': 'minimize'
        }
    }
    
    while True:
        print("\nOptions:")
        print("📚 Examples:")
        for key, example in examples.items():
            print(f"  {key}. {example['name']}")
        print("  6. Custom problem")
        print("  7. Exit")
        
        choice = input("\nSelect option (1-7): ").strip()
        
        if choice == '7':
            break
        elif choice in examples:
            example = examples[choice]
            print(f"\nLoaded: {example['name']}")
            
            df = optimizer.solve_complete(
                objective=example['objective'],
                constraints=example['constraints'],
                variables=example['variables'],
                problem_type=example['type'],
                visualize=True,
                export=False
            )
            
        elif choice == '6':
            print("\n🔧 Custom Problem Setup:")
            
            # Get problem type
            prob_type = input("Problem type (minimize/maximize) [minimize]: ").strip().lower()
            if prob_type not in ['minimize', 'maximize']:
                prob_type = 'minimize'
            
            # Get variables
            vars_input = input("Variables (space-separated, e.g., 'x y z'): ").strip()
            variables = vars_input.split()
            
            if not variables:
                print("❌ No variables provided!")
                continue
            
            # Get objective function
            objective = input(f"Objective function f({','.join(variables)}): ").strip()
            if not objective:
                print("❌ No objective function provided!")
                continue
            
            # Get constraints
            print("Enter constraints (one per line, empty line to finish):")
            constraints = []
            while True:
                constraint = input("  Constraint (=0): ").strip()
                if not constraint:
                    break
                constraints.append(constraint)
            
            if not constraints:
                print("❌ No constraints provided!")
                continue
            
            try:
                df = optimizer.solve_complete(
                    objective=objective,
                    constraints=constraints,
                    variables=variables,
                    problem_type=prob_type,
                    visualize=True,
                    export=True
                )
                
            except Exception as e:
                print(f"❌ Error solving problem: {e}")
        
        else:
            print("❌ Invalid choice!")
        
        input("\nPress Enter to continue...")


# Example usage and testing
if __name__ == "__main__":
    # Example 1: Basic optimization
    print("Example 1: Basic Constrained Optimization")
    print("-" * 50)
    
    optimizer = AdvancedLagrangeOptimizer()
    
    # Minimize x² + y² subject to x + y = 1
    df = optimizer.solve_complete(
        objective="x**2 + y**2",
        constraints=["x + y - 1"],
        variables=["x", "y"],
        problem_type="minimize",
        visualize=True,
        export=False
    )
    
    print("\n" + "="*70)
    print("Example 2: Economics - Utility Maximization")
    print("-" * 70)
    
    # Maximize utility x*y subject to budget constraint 2x + 3y = 12
    optimizer2 = AdvancedLagrangeOptimizer()
    df2 = optimizer2.solve_complete(
        objective="x*y",
        constraints=["2*x + 3*y - 12"],
        variables=["x", "y"],
        problem_type="maximize",
        visualize=True,
        export=False
    )
    
    print("\n" + "="*70)
    print("Example 3: Multi-constraint 3D Problem")
    print("-" * 70)
    
    # Minimize x² + y² + z² subject to x + y + z = 1 and x² + y² = 1
    optimizer3 = AdvancedLagrangeOptimizer()
    df3 = optimizer3.solve_complete(
        objective="x**2 + y**2 + z**2",
        constraints=["x + y + z - 1", "x**2 + y**2 - 1"],
        variables=["x", "y", "z"],
        problem_type="minimize",
        visualize=False,  # No 2D visualization for 3D problems
        export=False
    )
    
    # Advanced analysis and comparison
    print("\n" + "="*70)
    print("ADVANCED ANALYSIS AND COMPARISON")
    print("="*70)
    
    # Compare different optimization methods
    def compare_methods(objective, constraints, variables, problem_type="minimize"):
        """Compare symbolic vs numerical methods."""
        print(f"\nComparing methods for: {problem_type} {objective}")
        print(f"Subject to: {constraints}")
        
        optimizer_sym = AdvancedLagrangeOptimizer(verbose=False)
        optimizer_num = AdvancedLagrangeOptimizer(verbose=False)
        
        # Setup problems
        optimizer_sym.setup_problem(objective, constraints, variables, problem_type)
        optimizer_num.setup_problem(objective, constraints, variables, problem_type)
        
        # Solve symbolically
        sym_solutions = optimizer_sym.solve_symbolic()
        
        # Solve numerically
        num_solutions = optimizer_num.solve_numerical(num_random_starts=100)
        
        print(f"Symbolic solutions: {len(sym_solutions)}")
        print(f"Numerical solutions: {len(num_solutions)}")
        
        if sym_solutions and num_solutions:
            sym_obj = sym_solutions[0]['objective_value']
            num_obj = num_solutions[0]['objective_value']
            print(f"Best symbolic objective: {sym_obj:.8f}")
            print(f"Best numerical objective: {num_obj:.8f}")
            print(f"Difference: {abs(sym_obj - num_obj):.2e}")
        
        return sym_solutions, num_solutions
    
    # Test various problems
    test_problems = [
        ("x**2 + y**2", ["x + y - 1"], ["x", "y"], "minimize"),
        ("x*y", ["2*x + 3*y - 12"], ["x", "y"], "maximize"),
        ("(x-1)**2 + (y-2)**2", ["x**2 + y**2 - 1"], ["x", "y"], "minimize"),
        ("x**2 + 2*y**2", ["x**2 + y**2 - 4"], ["x", "y"], "minimize")
    ]
    
    for obj, const, vars, prob_type in test_problems:
        compare_methods(obj, const, vars, prob_type)
    
    # Sensitivity analysis
    print("\n" + "="*70)
    print("SENSITIVITY ANALYSIS")
    print("="*70)
    
    def sensitivity_analysis():
        """Analyze how solutions change with constraint parameters."""
        print("Analyzing sensitivity to constraint parameter changes...")
        
        # Base problem: minimize x² + y² subject to x + y = c
        results = []
        c_values = np.linspace(0.5, 2.0, 11)
        
        for c in c_values:
            optimizer = AdvancedLagrangeOptimizer(verbose=False)
            constraint = f"x + y - {c}"
            
            try:
                optimizer.solve_complete(
                    objective="x**2 + y**2",
                    constraints=[constraint],
                    variables=["x", "y"],
                    problem_type="minimize",
                    visualize=False,
                    export=False
                )
                
                if optimizer.solutions:
                    sol = optimizer.solutions[0]
                    results.append({
                        'c': c,
                        'x': sol['x'],
                        'y': sol['y'],
                        'objective': sol['objective_value'],
                        'lambda': sol.get('λ1', 0)
                    })
            except:
                continue
        
        if results:
            df_sensitivity = pd.DataFrame(results)
            print("\nSensitivity Analysis Results:")
            print(df_sensitivity.round(6))
            
            # Plot sensitivity
            fig, axes = plt.subplots(2, 2, figsize=(12, 10))
            
            axes[0, 0].plot(df_sensitivity['c'], df_sensitivity['x'], 'bo-', label='x*')
            axes[0, 0].plot(df_sensitivity['c'], df_sensitivity['y'], 'ro-', label='y*')
            axes[0, 0].set_xlabel('Constraint parameter c')
            axes[0, 0].set_ylabel('Optimal values')
            axes[0, 0].set_title('Optimal Variables vs Constraint Parameter')
            axes[0, 0].legend()
            axes[0, 0].grid(True, alpha=0.3)
            
            axes[0, 1].plot(df_sensitivity['c'], df_sensitivity['objective'], 'go-')
            axes[0, 1].set_xlabel('Constraint parameter c')
            axes[0, 1].set_ylabel('Optimal objective value')
            axes[0, 1].set_title('Optimal Objective vs Constraint Parameter')
            axes[0, 1].grid(True, alpha=0.3)
            
            axes[1, 0].plot(df_sensitivity['c'], df_sensitivity['lambda'], 'mo-')
            axes[1, 0].set_xlabel('Constraint parameter c')
            axes[1, 0].set_ylabel('Lagrange multiplier λ')
            axes[1, 0].set_title('Lagrange Multiplier vs Constraint Parameter')
            axes[1, 0].grid(True, alpha=0.3)
            
            # Shadow price interpretation
            axes[1, 1].plot(df_sensitivity['c'][:-1], 
                           np.diff(df_sensitivity['objective']) / np.diff(df_sensitivity['c']), 
                           'co-', label='Numerical derivative')
            axes[1, 1].plot(df_sensitivity['c'], df_sensitivity['lambda'], 'mo--', 
                           label='Lagrange multiplier', alpha=0.7)
            axes[1, 1].set_xlabel('Constraint parameter c')
            axes[1, 1].set_ylabel('Shadow price (∂f*/∂c)')
            axes[1, 1].set_title('Shadow Price Analysis')
            axes[1, 1].legend()
            axes[1, 1].grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.show()
    
    sensitivity_analysis()
    
    # Performance benchmark
    print("\n" + "="*70)
    print("PERFORMANCE BENCHMARK")
    print("="*70)
    
    def benchmark_methods():
        """Benchmark different solution methods."""
        import time
        
        problems = [
            ("Simple", "x**2 + y**2", ["x + y - 1"], ["x", "y"]),
            ("Nonlinear", "x**2 + y**4", ["x**2 + y**2 - 1"], ["x", "y"]),
            ("Multi-constraint", "x**2 + y**2 + z**2", 
             ["x + y + z - 1", "x**2 + y**2 - 1"], ["x", "y", "z"]),
        ]
        
        benchmark_results = []
        
        for name, obj, const, vars in problems:
            print(f"\nBenchmarking: {name}")
            
            # Symbolic method
            optimizer_sym = AdvancedLagrangeOptimizer(verbose=False)
            optimizer_sym.setup_problem(obj, const, vars)
            
            start_time = time.time()
            sym_solutions = optimizer_sym.solve_symbolic()
            sym_time = time.time() - start_time
            
            # Numerical method
            optimizer_num = AdvancedLagrangeOptimizer(verbose=False)
            optimizer_num.setup_problem(obj, const, vars)
            
            start_time = time.time()
            num_solutions = optimizer_num.solve_numerical(num_random_starts=50)
            num_time = time.time() - start_time
            
            benchmark_results.append({
                'Problem': name,
                'Symbolic_Time': sym_time,
                'Symbolic_Solutions': len(sym_solutions),
                'Numerical_Time': num_time,
                'Numerical_Solutions': len(num_solutions),
                'Speedup': num_time / sym_time if sym_time > 0 else np.inf
            })
            
            print(f"  Symbolic: {len(sym_solutions)} solutions in {sym_time:.4f}s")
            print(f"  Numerical: {len(num_solutions)} solutions in {num_time:.4f}s")
        
        df_benchmark = pd.DataFrame(benchmark_results)
        print("\nBenchmark Summary:")
        print(df_benchmark)
    
    benchmark_methods()
    
    # Start interactive mode
    print("\n" + "="*70)
    print("INTERACTIVE MODE")
    print("="*70)
    print("Starting interactive solver...")
    
    try:
        interactive_solver()
    except KeyboardInterrupt:
        print("\n\nExiting interactive solver. Thank you!")
    
    print("\n🎯 Advanced Lagrange Multiplier Solver - Session Complete!")
    print("="*70)
    
    # Additional utility functions
    def create_optimization_report(optimizer, filename=None):
        """Create a comprehensive optimization report."""
        if not filename:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"optimization_report_{timestamp}.html"
        
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>Optimization Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 40px; }}
                .header {{ background-color: #f0f0f0; padding: 20px; border-radius: 5px; }}
                .solution {{ background-color: #e8f4fd; padding: 15px; margin: 10px 0; border-radius: 5px; }}
                table {{ border-collapse: collapse; width: 100%; }}
                th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
                th {{ background-color: #f2f2f2; }}
                .success {{ color: green; }}
                .warning {{ color: orange; }}
                .error {{ color: red; }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>🎯 Lagrange Multiplier Optimization Report</h1>
                <p><strong>Generated:</strong> {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
                <p><strong>Problem Type:</strong> {optimizer.problem_type.capitalize()}</p>
                <p><strong>Objective:</strong> {optimizer.objective}</p>
                <p><strong>Constraints:</strong> {[str(c) for c in optimizer.constraints]}</p>
                <p><strong>Variables:</strong> {[str(v) for v in optimizer.variables]}</p>
            </div>
            
            <h2>📊 Solutions Found: {len(optimizer.solutions)}</h2>
        """
        
        for i, sol in enumerate(optimizer.solutions):
            html_content += f"""
            <div class="solution">
                <h3>Solution {i+1}</h3>
                <table>
                    <tr><th>Variable/Parameter</th><th>Value</th></tr>
            """
            
            for var in optimizer.variables:
                html_content += f"<tr><td>{var}</td><td>{sol.get(str(var), 'N/A'):.6f}</td></tr>"
            
            for mult in optimizer.multipliers:
                html_content += f"<tr><td>{mult}</td><td>{sol.get(str(mult), 'N/A'):.6f}</td></tr>"
            
            html_content += f"<tr><td><strong>Objective Value</strong></td><td><strong>{sol['objective_value']:.6f}</strong></td></tr>"
            html_content += "</table>"
            
            # Constraint verification
            html_content += "<h4>Constraint Verification:</h4><ul>"
            for j, const_val in enumerate(sol.get('constraint_values', [])):
                status = "✓ Satisfied" if abs(const_val) < 1e-6 else "✗ Violated"
                color = "success" if abs(const_val) < 1e-6 else "error"
                html_content += f'<li class="{color}">Constraint {j+1}: {const_val:.2e} {status}</li>'
            html_content += "</ul></div>"
        
        html_content += """
            </body>
            </html>
        """
        
        with open(filename, 'w') as f:
            f.write(html_content)
        
        print(f"📄 Comprehensive report saved to {filename}")
    
    def batch_solve(problem_list, output_dir="lagrange_results"):
        """Solve multiple optimization problems in batch."""
        import os
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        results = []
        
        for i, problem in enumerate(problem_list):
            print(f"\nSolving problem {i+1}/{len(problem_list)}: {problem.get('name', 'Unnamed')}")
            
            try:
                optimizer = AdvancedLagrangeOptimizer(verbose=False)
                df = optimizer.solve_complete(
                    objective=problem['objective'],
                    constraints=problem['constraints'],
                    variables=problem['variables'],
                    problem_type=problem.get('type', 'minimize'),
                    visualize=False,
                    export=False
                )
                
                # Save individual results
                problem_name = problem.get('name', f'problem_{i+1}').replace(' ', '_')
                optimizer.export_results(f"{output_dir}/{problem_name}_results.json")
                create_optimization_report(optimizer, f"{output_dir}/{problem_name}_report.html")
                
                results.append({
                    'name': problem.get('name', f'Problem {i+1}'),
                    'status': 'Success',
                    'num_solutions': len(optimizer.solutions),
                    'best_objective': min(s['objective_value'] for s in optimizer.solutions) if optimizer.solutions else None
                })
                
            except Exception as e:
                results.append({
                    'name': problem.get('name', f'Problem {i+1}'),
                    'status': f'Failed: {str(e)}',
                    'num_solutions': 0,
                    'best_objective': None
                })
        
        # Save batch summary
        batch_df = pd.DataFrame(results)
        batch_df.to_csv(f"{output_dir}/batch_summary.csv", index=False)
        print(f"\n📊 Batch processing complete! Results saved to {output_dir}/")
        print(batch_df)
        
        return batch_df

# Example batch processing
if __name__ == "__main__":
    # Define a batch of problems to solve
    batch_problems = [
        {
            'name': 'Basic Quadratic',
            'objective': 'x**2 + y**2',
            'constraints': ['x + y - 1'],
            'variables': ['x', 'y'],
            'type': 'minimize'
        },
        {
            'name': 'Utility Maximization',
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
    
    # Uncomment to run batch processing
    # batch_results = batch_solve(batch_problems)