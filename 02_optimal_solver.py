import numpy as np
from scipy.optimize import minimize
import time, warnings

# Start the timer
start_time = time.time()

def cos(x):
    """Convert angle from degrees to radians and return its cosine."""
    return np.cos(np.deg2rad(x, dtype = 'float64'))

def sin(x):
    """Convert angle from degrees to radians and return its sine."""
    return np.sin(np.deg2rad(x, dtype = 'float64'))

def tan(x):
    """Convert angle from degrees to radians and return its tangent."""
    return np.tan(np.deg2rad(x, dtype = 'float64'))

def coefficients(delta, x, y):
    """
    Calculate coefficients A and B based on delta, x, and y.

    Parameters:
    delta (float): Angle in degrees.
    x, y (float): Parameters.

    Returns:
    tuple: A1, A2, A3, A4, B1, B2, B3, B4
    """
    # Calculate intermediate values
    B2 = np.power(x, 2)
    B3 = np.power(y, 2)
    B1 = B2 * B3
    B4h = x * y
    
    # Set coefficients based on delta value
    if delta == 90:
        A1, A2, A3, A4 = 1, -1, -1, 1
    else:
        L = cos(delta)
        k = B4h / L 
        A1 = B1 - k
        A2 = B2 + k
        A3 = B3 + k
        A4 = 1 - k
    
    B4 = 2 * B4h
    
    return A1, A2, A3, A4, B1, B2, B3, B4

def generate_inequality_constraints(vars):
    """
    Generate the inequality constraints based on the input variables.

    Parameters:
    vars (list): A list of 9 variables (x1, x2, ..., x9).

    Returns:
    np.array: An array of 4 inequality values.
    """
    # Unpack the input variables
    x1, x2, x3, x4, x5, x6, x7, x8, x9 = vars

    # Define the inequality expressions
    ineq1 = x1 * x3 * x5 * (1 + x1) * (1 + x3) * (1 + x5)
    ineq2 = x1 * x4 * x6 * (1 + x1) * (1 + x4) * (1 + x6)
    ineq3 = x2 * x4 * x7 * (1 + x2) * (1 + x4) * (1 + x7)
    ineq4 = x2 * x3 * x8 * (1 + x2) * (1 + x3) * (1 + x8)

    # Return the array of constraints
    return np.array([ineq1, ineq2, ineq3, ineq4])

def inequality_constraints_wrapper():
    """
    Creates a list of inequality constraint dictionaries for the optimization.

    Returns:
    list: A list of 8 dictionaries, each representing an inequality constraint.
    """
    # Create a lambda function for each inequality constraint
    return [{'type': 'ineq', 'fun': lambda vars, idx=i: generate_inequality_constraints(vars)[idx]} for i in range(4)]

def random_with_inequality_constraints():
    """
    Generate random variables that satisfy the inequality constraints.

    Uses scipy's minimize function to enforce the constraints and random initialization.

    Returns:
    np.array: The optimized variables satisfying the constraints.
    """
    # Define a constant objective function (as we're only interested in satisfying the constraints)
    def objective(_):
        return 0  # Objective function always returns 0

    # Helper function to generate an initial guess for x1 to x8 from a normal distribution
    # and x9 from a uniform distribution based on the provided bounds
    def generate_initial_guess(bounds):
        """Generate initial guesses for x1 to x9 based on provided bounds."""
        # x1 to x8 are normally distributed around 0 with a standard deviation of 100
        # x9 is generated from a uniform distribution in the provided bounds
        return [np.random.normal(0, 100) for _ in range(8)] + [float(np.random.uniform(*bounds))]

    # Get the list of inequality constraints
    constraints = inequality_constraints_wrapper()

    # Mapping `switch` values to appropriate bounds for x9
    bounds_map = {
        11: (0, 100),   # For switch values 11 and 12, x9 is in the range [0, 100)
        12: (0, 100),
        21: (-100, 0),  # For switch values 21 and 22, x9 is in the range [-100, 0)
        22: (-100, 0)
    }

    # Check if the current switch value is one of the predefined values
    if switch in bounds_map:
        bounds = bounds_map[switch]  # Get the bounds for the current `switch` value

        # Continuously try generating variables until the constraints are satisfied
        while True:
            # Generate a random initial guess for x1 to x9 based on the bounds
            initial_guess = generate_initial_guess(bounds)
            
            # Run the optimization with the random initial guess, trying to satisfy the constraints
            result = minimize(objective, initial_guess, constraints=constraints)
            
            # Check if the generated solution satisfies the constraints
            # The function `generate_inequality_constraints` returns values for all constraints
            # We check if all the values are positive (which means the constraints are met)
            if np.all(generate_inequality_constraints(result.x) > 0):
                return result.x  # Return the result if constraints are satisfied

def generate_constraints(vars):
    """
    Generate the system of equations based on input variables.

    Parameters:
    vars (list): List of 9 variables.

    Returns:
    numpy array: Array of 9 constraint equations.
    """
    # Unpack variables
    x1, x2, x3, x4, x5, x6, x7, x8, x9 = vars
    
    # Define the equations using the variables and coefficients
    eq1 = A11 + A12 * x3 * x5 * x9 + A13 * x1 * x5 * x9 + A14 * x1 * x3 * x9
    eq2 = A21 + A22 * x4 * x6 * x9 + A23 * x1 * x6 * x9 + A24 * x1 * x4 * x9
    eq3 = A31 + A32 * x4 * x7 * x9 + A33 * x2 * x7 * x9 + A34 * x2 * x4 * x9
    eq4 = A41 + A42 * x3 * x8 * x9 + A43 * x2 * x8 * x9 + A44 * x2 * x3 * x9

    # More complex equations involving x1, x3, x4, x5, etc.
    eq5 = (B11 + B12 * x3 * x5 * x9 + B13 * x1 * x5 * x9 + x1 * x3 * x9)**2 - B14**2 * (1 + x5) * (1 + x5 * x9) * x1 * x3 * x9
    eq6 = (B21 + B22 * x4 * x6 * x9 + B23 * x1 * x6 * x9 + x1 * x4 * x9)**2 - B24**2 * (1 + x6) * (1 + x6 * x9) * x1 * x4 * x9
    eq7 = (B31 + B32 * x4 * x7 * x9 + B33 * x2 * x7 * x9 + x2 * x4 * x9)**2 - B34**2 * (1 + x7) * (1 + x7 * x9) * x2 * x4 * x9
    eq8 = (B41 + B42 * x3 * x8 * x9 + B43 * x2 * x8 * x9 + x2 * x3 * x9)**2 - B44**2 * (1 + x8) * (1 + x8 * x9) * x2 * x3 * x9

    # Constraints for M < 1 (x9 > 0)
    eq91 = x5 * x6 * x9 - 1
    eq92 = x7 * x8 * x9 - 1
    eq93 = (((x5 - x6) * (x7 * x8 * x9 - 1) + (x7 - x8) * (x5 * x6 * x9 - 1))**2 - 4 * (x5 * x7 * (1 + x6) * (1 + x8) * (1 + x6 * x9) * (1 + x8 * x9) + x6 * x8 * (1 + x5) * (1 + x7) * (1 + x5 * x9) * (1 + x7 * x9)))**2 - 64 * x5 * x6 * x7 * x8 * (1 + x5) * (1 + x6) * (1 + x7) * (1 + x8) * (1 + x5 * x9) * (1 + x6 * x9) * (1 + x7 * x9) * (1 + x8 * x9)

    # Constraints for M > 1 (x9 < 0)
    eq94 = x5 * x9 + x6 * x9 + x5 * x6 * x9 + 1
    eq95 = x7 * x9 + x8 * x9 + x7 * x8 * x9 + 1
    eq96 = (((x5 - x6) * (x7 * x9 + x8 * x9 + x7 * x8 * x9 + 1) + (x7 - x8) * (x5 * x9 + x6 * x9 + x5 * x6 * x9 + 1))**2 - 4 * (x5 * x7 * (1 + x6) * (1 + x8) * (1 + x5 * x9) * (1 + x7 * x9) + x6 * x8 * (1 + x5) * (1 + x7) * (1 + x6 * x9) * (1 + x8 * x9)))**2 - 64 * x5 * x6 * x7 * x8 * (1 + x5) * (1 + x6) * (1 + x7) * (1 + x8) * (1 + x5 * x9) * (1 + x6 * x9) * (1 + x7 * x9) * (1 + x8 * x9)

    # Mapping switch to the correct equation
    switch_map = {
        11: eq91**2 + eq92**2,
        12: eq93,
        21: eq94**2 + eq95**2,
        22: eq96
    }

    eq9 = switch_map.get(switch)  # Get the appropriate equation based on switch
    
    # Return all equations as an array
    return np.array([eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9])

def constraints_wrapper():
    """
    Create and return a list of constraint dictionaries for optimization.

    Returns:
    list: List of constraint dictionaries.
    """
    # Helper function to generate constraints
    def constraint_factory(index):
        return lambda vars: generate_constraints(vars)[index]
    
    # Create a list of constraints
    constraints_list = [{'type': 'eq', 'fun': constraint_factory(i)} for i in range(9)]
    return constraints_list

def objective(vars):
    """
    Define the objective function to be minimized.

    Parameters:
    vars (list): List of 9 variables.

    Returns:
    float: The value of the objective function.
    """
    # Unpack variables
    x1, x2, x3, x4, x5, x6, x7, x8, x9 = vars
    
    # Define the parts of the objective function
    p1 = cos(alpha1)**2 - (1 - x3 * x5 * x9 + x1 * x5 * x9 - x1 * x3 * x9)**2 / (4 * x1 * x5 * x9 * (1 + x3) * (1 + x3 * x9))
    p2 = cos(alpha2)**2 - (1 - x4 * x6 * x9 + x1 * x6 * x9 - x1 * x4 * x9)**2 / (4 * x1 * x6 * x9 * (1 + x4) * (1 + x4 * x9))
    p3 = cos(alpha3)**2 - (1 - x4 * x7 * x9 + x2 * x7 * x9 - x2 * x4 * x9)**2 / (4 * x2 * x7 * x9 * (1 + x4) * (1 + x4 * x9))
    p4 = cos(alpha4)**2 - (1 - x3 * x8 * x9 + x2 * x8 * x9 - x2 * x3 * x9)**2 / (4 * x2 * x8 * x9 * (1 + x3) * (1 + x3 * x9))

    # Define the parts involving gamma angles
    p5 = cos(gamma1)**2 - (1 + x3 * x5 * x9 - x1 * x5 * x9 - x1 * x3 * x9)**2 / (4 * x3 * x5 * x9 * (1 + x1) * (1 + x1 * x9))
    p6 = cos(gamma2)**2 - (1 + x4 * x6 * x9 - x1 * x6 * x9 - x1 * x4 * x9)**2 / (4 * x4 * x6 * x9 * (1 + x1) * (1 + x1 * x9))
    p7 = cos(gamma3)**2 - (1 + x4 * x7 * x9 - x2 * x7 * x9 - x2 * x4 * x9)**2 / (4 * x4 * x7 * x9 * (1 + x2) * (1 + x2 * x9))
    p8 = cos(gamma4)**2 - (1 + x3 * x8 * x9 - x2 * x8 * x9 - x2 * x3 * x9)**2 / (4 * x3 * x8 * x9 * (1 + x2) * (1 + x2 * x9))
    
    # Return the sum of squares of all parts
    return p1**2 + p2**2 + p3**2 + p4**2 + p5**2 + p6**2 + p7**2 + p8**2

def optimal_solution():
    """
    Find the optimal solution using a constraint-based minimization.

    Returns:
    dict: Dictionary containing optimal solutions.
    """
    # Get the list of constraints
    cons = constraints_wrapper()

    # Set k_cons based on the switch value
    k_cons = 1 if switch in {11, 12} else -1 if switch in {21, 22} else None

    dictionary = {}  # Store solutions
    eps = 1  # Error threshold
    i = 0
    
    while i < 2000 and eps > 0.01:  # Run loop until we hit error threshold or 2000 iterations
        np.random.seed(i)  # Seed randomness for reproducibility
        
        # Generate the initial guess for the variables
        initial_guess = random_with_inequality_constraints()
        
        # Perform the minimization using the 'trust-constr' method
        # solution = minimize(objective, initial_guess, method='trust-constr', constraints=cons)
        
        # Perform the minimization using the 'SLSQP' method
        solution = minimize(objective, initial_guess, method='SLSQP', constraints=cons)

        # Unpack solution variables
        x1, x2, x3, x4, x5, x6, x7, x8, x9 = solution.x

        # Create the list to compare and store
        listX = [objective(solution.x), x1, x2, x3, x4, x5, x6, x7, x8, x9]

        # Check if the solution is already in the dictionary
        exists = any(np.allclose(values, listX, atol=1e-5) for values in dictionary.values())  
              
        # Check if the constraints are satisfied and solution is new
        if np.all(np.isclose(generate_constraints([x1, x2, x3, x4, x5, x6, x7, x8, x9]), 0, atol=1e-5)) and x9 * k_cons > 0 and not exists:
            with warnings.catch_warnings():
                warnings.filterwarnings('error')
                try:
                    # Calculate cosSigma values for the solution
                    cosSigma1 = (1 - x9 * (x1 * x3 + x1 * x5 + x3 * x5 + 2 * x1 * x3 * x5)) / (x9 * np.sqrt(4 * x1 * x3 * x5 * (1 + x1) * (1 + x3) * (1 + x5)))
                    cosSigma2 = (1 - x9 * (x1 * x4 + x1 * x6 + x4 * x6 + 2 * x1 * x4 * x6)) / (x9 * np.sqrt(4 * x1 * x4 * x6 * (1 + x1) * (1 + x4) * (1 + x6)))
                    cosSigma3 = (1 - x9 * (x2 * x4 + x2 * x7 + x4 * x7 + 2 * x2 * x4 * x7)) / (x9 * np.sqrt(4 * x2 * x4 * x7 * (1 + x2) * (1 + x4) * (1 + x7)))
                    cosSigma4 = (1 - x9 * (x2 * x3 + x2 * x8 + x3 * x8 + 2 * x2 * x3 * x8)) / (x9 * np.sqrt(4 * x2 * x3 * x8 * (1 + x2) * (1 + x3) * (1 + x8)))
                    
                    # Check if cosSigma values are valid
                    if np.all(abs(cosSigma) < 1 for cosSigma in [cosSigma1, cosSigma2, cosSigma3, cosSigma4]):
                        dictionary[i] = listX  # Store solution
                except:
                    pass
            
        i += 1  # Increment iteration
    
    # Sort the dictionary by the objective function values
    sorted_dictionary = dict(sorted(dictionary.items(), key=lambda item: item[1][0]))

    # Get the current time for the file name
    now = time.localtime()
    current_time = time.strftime("%Y_%m_%d_%H_%M_%S", now)
    
    # Write results to a file
    with open(f"roots_{current_time}.txt", "a") as rootsFile:
        rootsFile.write(f"Switch = {switch}\n\n")
        for values in sorted_dictionary.values():
            eps, *x_vars = values
            text = (f"objective = {eps}\n"+
                    "\n".join([f"x{i+1} = {x_vars[i]}" for i in range(9)]) + "\n\n\n")
            rootsFile.write(text)

    return sorted_dictionary

def read_sample_angles_from_file(file_path, option_number):
    """
    Reads the file and extracts the lists Alphas, Betas, Gammas, Deltas, and DihedralAngles for the specified option.

    Parameters:
    - file_path (str): The path to the file containing the angles.
    - option_number (int): The option number to extract (1, 2, 3, etc.).

    Returns:
    - tuple: A tuple containing lists (Alphas, Betas, Gammas, Deltas, DihedralAngles).
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Initialize variables to hold the parsed data
    Alphas, Betas, Gammas, Deltas, DihedralAngles = [], [], [], [], []
    option_found = False  # Flag to track if we've found the right option

    # Look for the selected option and extract the values
    current_option = 0
    for i, line in enumerate(lines):
        line = line.strip()
        
        # Look for "Option X" to identify the correct section
        if line.startswith("Option"):
            current_option = int(line.split()[1])
            if current_option == option_number:
                option_found = True
                continue  # Move to the next lines

        if option_found:
            if line.startswith("Alphas"):
                Alphas = eval(line.split('=')[1].strip())
            elif line.startswith("Betas"):
                Betas = eval(line.split('=')[1].strip())
            elif line.startswith("Gammas"):
                Gammas = eval(line.split('=')[1].strip())
            elif line.startswith("Deltas"):
                Deltas = eval(line.split('=')[1].strip())
            elif line.startswith("DihedralAngles"):
                DihedralAngles = eval(line.split('=')[1].strip())
                break  # We found everything we need, no need to process further lines

    if not option_found:
        raise ValueError(f"Option {option_number} not found in the file.")

    # Ensure that all angles are found
    if not (Alphas and Betas and Gammas and Deltas and DihedralAngles):
        raise ValueError(f"Incomplete data for Option {option_number}: Alphas, Betas, Gammas, Deltas, or DihedralAngles are missing.")

    return Alphas, Betas, Gammas, Deltas, DihedralAngles

def main():
    """
    Main function to execute the optimization process.

    This function reads sample angles from a file, assigns them to global variables, 
    computes various parameters and coefficients, and finally finds the optimal solution 
    based on the defined equations and constraints.
    """
    # Set the variant of eq9 to be used in the optimization
    global switch
    switch = 22

    # Load sample values for Alphas, Betas, Gammas, Deltas, and DihedralAngles from the file
    # You can select which option to load by passing the option number to 'read_angles_from_file'
    AlphasS, BetasS, GammasS, DeltasS, DihedralAnglesS = read_sample_angles_from_file("sample_angles.txt", 1)
    
    # Print the loaded values for verification (you can remove this in production)
    print(AlphasS, BetasS, GammasS, DeltasS, DihedralAnglesS)
    
    # Unpack and assign the Alphas values to global variables
    global alpha1, alpha2, alpha3, alpha4
    alpha1, alpha2, alpha3, alpha4 = AlphasS
    
    # Unpack and assign the Betas values to global variables
    global beta1, beta2, beta3, beta4
    beta1, beta2, beta3, beta4 = BetasS
    
    # Unpack and assign the Gammas values to global variables
    global gamma1, gamma2, gamma3, gamma4
    gamma1, gamma2, gamma3, gamma4 = GammasS

    # Unpack and assign the Deltas values to global variables
    global delta1, delta2, delta3, delta4
    delta1, delta2, delta3, delta4 = DeltasS

    # Unpack and assign the DihedralAngles values to local variables
    phi, psi2, theta, psi1 = DihedralAnglesS
    
    # Compute parameters based on the dihedral angles
    # These parameters are the tangents of half-angles: tan(angle / 2)
    z, w1, w2, u = map(lambda angle: tan(angle / 2), [phi, psi1, psi2, theta])
    
    # Calculate the coefficients based on delta values and parameters (z, w1, w2, u)
    # The coefficients are used later for optimization
    global A11, A12, A13, A14, B11, B12, B13, B14
    global A21, A22, A23, A24, B21, B22, B23, B24
    global A31, A32, A33, A34, B31, B32, B33, B34
    global A41, A42, A43, A44, B41, B42, B43, B44
    
    # Compute the coefficients for each of the four delta values
    A11, A12, A13, A14, B11, B12, B13, B14 = coefficients(delta1, z, w1)
    A21, A22, A23, A24, B21, B22, B23, B24 = coefficients(delta2, z, w2)
    A31, A32, A33, A34, B31, B32, B33, B34 = coefficients(delta3, u, w2)
    A41, A42, A43, A44, B41, B42, B43, B44 = coefficients(delta4, u, w1)

    # Find the optimal solution by solving the system based on the computed coefficients
    dX = optimal_solution()

    # Get the key with the minimum objective value from the solution dictionary
    min_key = min(dX, key=dX.get)
    
    # Print the optimal solution and the number of solutions evaluated
    print(min_key, dX[min_key], len(dX))
    
    # End the timer and calculate the elapsed time in minutes
    end_time = time.time()
    elapsed_time = (end_time - start_time) / 60
    print(f"Elapsed time: {elapsed_time:.2f} minutes")

if __name__ == "__main__":
    main()