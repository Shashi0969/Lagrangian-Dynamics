import sympy as sp

def lagrange_multiplier_method(objective, constraint, variables):
    """
    Solves a constrained optimization problem using the Lagrange multiplier method.
    
    Parameters:
    objective (sympy expression): The objective function to optimize.
    constraint (sympy expression): The constraint equation (set to zero).
    variables (list of sympy symbols): List of variables in the objective function.
    
    Returns:
    list of dicts: Solutions (including Lagrange multipliers) satisfying the conditions.
    """
    # Create Lagrange multiplier symbol
    lmbda = sp.symbols('λ')
    # Construct the Lagrangian: L = objective - λ * constraint
    lagrangian = objective - lmbda * constraint
    
    # Compute the gradient of the Lagrangian
    gradient = [sp.diff(lagrangian, var) for var in variables]
    # Add the constraint equation to the system
    equations = gradient + [constraint]
    
    # Solve the system of equations
    solutions = sp.solve(equations, variables + [lmbda], dict=True)
    return solutions

# Example usage
if __name__ == "__main__":
    """ Call the function with an example """
    # Define variables
    x, y = sp.symbols('x y ')
    
    # Objective function: f(x, y) = x^2 + y^2
    objective_function = input("Enter the objective function to be optiimized: ") # x**2 + y**2
    
    # Constraint: g(x, y) = x + y - 1 = 0
    constraint_equation = x + y - 1
    
    # Apply the Lagrange multiplier method
    solutions = lagrange_multiplier_method(objective_function, constraint_equation, [x, y])
    
    # Print the solutions
    print("Solutions:")
    for sol in solutions:
        print(sol)
