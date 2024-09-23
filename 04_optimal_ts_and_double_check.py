import mpmath as mp
import numpy as np
import ast

def sin(x):
    """Convert angle from degrees to radians and return its sine."""
    return np.sin(np.deg2rad(x, dtype = 'float64'))

def inverse_jacobi_dn_complex(x, m, tol=1e-6, max_attempts=20):
    """
    Numerically find u (real or complex) such that dn(u, m) = x, using the correct periodicity
    and random initial guesses if needed. This is the inverse Jacobi Dn function.

    Parameters:
        x (float or complex): The value of dn(u, m) where we are solving for u.
        m (float or complex): The parameter of the Jacobi elliptic function.
        tol (float): The tolerance for the root-finding process.
        max_attempts (int): Maximum number of random guesses.

    Returns:
        complex: The value of u such that dn(u, m) = x, potentially complex.
    """
    def jacobi_dn(u):
        return mp.ellipfun('dn', u, m)

    K = mp.ellipk(m)  # Real period
    Kprime = mp.ellipk(1 - m)  # Imaginary period

    i = 0 
    for attempt in range(max_attempts):
        np.random.seed(i)
        i += 1
        real_guess = np.random.uniform(0, 2 * K)
        imag_guess = np.random.uniform(0, 4 * Kprime)
        initial_guess = complex(real_guess, imag_guess)

        try:
            u_solution = mp.findroot(lambda u: jacobi_dn(u) - x, initial_guess, tol=tol)

            real_part = u_solution.real % (2 * K)
            imag_part = u_solution.imag % (4 * Kprime)
            u_solution_normalized = real_part + 1j * imag_part

            return u_solution_normalized
        except ValueError:
            continue

    raise ValueError("Could not find root within the given number of attempts.")

def conditional_round(number, tolerance=1e-5):
    """
    Rounds the number only if it is within a specified tolerance of its rounded value.

    Parameters:
        number (float): The number to conditionally round.
        tolerance (float): The tolerance within which the number will be rounded. Default is 1e-5.

    Returns:
        float: The rounded number if close enough, otherwise the original number.
    """
    rounded_number = round(number)
    
    # Check if the number is close to its rounded value within the tolerance
    if np.isclose(number, rounded_number, atol=tolerance):
        return rounded_number
    else:
        return number

def calculate_elliptic_params(M):
    """
    Calculate elliptic parameters based on M.

    Parameters:
        M (float): Input parameter M.

    Returns:
        tuple: m, K, and Kprime, where K is the elliptic integral and Kprime is the complementary elliptic integral.
    """
    if M < 1:
        m = 1 - M
    else:
        m = (M - 1) / M

    K = mp.ellipk(m)
    Kprime = mp.ellipk(1 - m)
    return m, K, Kprime

def calculate_t(f, M, m):
    """
    Calculate the value of t using inverse Jacobi elliptic functions.

    Parameters:
        f (float): A parameter in the calculation.
        M (float): Input parameter.
        m (float): Elliptic modulus.

    Returns:
        complex: The value of t.
    """
    if M > 1:
        return inverse_jacobi_dn_complex(1 / mp.sqrt(f), m)
    else:
        return inverse_jacobi_dn_complex(mp.sqrt(f), m)

def handle_complex_values(t, K, Kprime, re_condition):
    """
    Adjust and handle the complex values using elliptic integrals.

    Parameters:
        t (complex): The result of the elliptic function.
        K (float): Elliptic integral.
        Kprime (float): Complementary elliptic integral.
        re_condition (float): Condition for real part of t.

    Returns:
        complex: The adjusted complex value.
    """
    real_part = float(mp.re(t)) / float(K)
    imag_part = float(mp.im(t)) / float(Kprime)

    if np.isclose(conditional_round(real_part) % 2, re_condition % 2, atol=1e-5):
        return complex(re_condition * K, adjust_imaginary(imag_part) * Kprime)
    else:
        return complex(real_part * K, adjust_imaginary(imag_part) * Kprime)

def adjust_imaginary(imag_part):
    """
    Adjust the imaginary part of a complex number based on its value.

    Parameters:
        imag_part (float): The imaginary part of the number.

    Returns:
        float: The adjusted imaginary part.
    """
    if imag_part < 0:
        imag_part += 4  # Shift negative values into the range [0, 4)

    if 2 <= imag_part < 4:
        return 4 - imag_part  # Reflect values in the range [2, 4)

    return imag_part  # Return values already in the range [0, 2)

def t_values(M, alpha, beta, gamma, delta, r, s, f):
    """
    Calculate t-values based on input parameters.

    Parameters:
        M (float): Input parameter.
        alpha, beta, gamma, delta (float): Angles in degrees.
        r, s, f (float): Other parameters used in the calculation.

    Returns:
        complex: Calculated t-value.
    """
    sigma = (alpha + beta + gamma + delta) / 2
    m, K, Kprime = calculate_elliptic_params(M)

    if sigma < 180:
        if r > 1 and s > 1:
            t = calculate_t(f, M, m)
            return handle_complex_values(t, K, Kprime, 0)
        if (r > 1 and s < 1) or (r < 1 and s > 1):
            t = calculate_t(f, M, m)
            return handle_complex_values(t, K, Kprime, 1)
        if r < 1 and s < 1:
            t = calculate_t(f, M, m)
            return handle_complex_values(t, K, Kprime, 2)
    elif sigma > 180:
        if r > 1 and s > 1:
            t = calculate_t(f, M, m)
            return handle_complex_values(t, K, Kprime, 2)
        if (r > 1 and s < 1) or (r < 1 and s > 1):
            t = calculate_t(f, M, m)
            return handle_complex_values(t, K, Kprime, 3)
        if r < 1 and s < 1:
            t = calculate_t(f, M, m)
            return handle_complex_values(t, K, Kprime, 0)

    return None

def read_values_from_file(file_path):
    """
    Read the values of Alphas, Betas, Gammas, Deltas, and x_values from the specified text file.

    Parameters:
        file_path (str): Path to the input file.

    Returns:
        dict: Dictionary containing lists of Alphas, Betas, Gammas, Deltas, and x_values.
    """
    data = {}
    
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('Switch'):
                data['Switch'] = int(line.split('=')[1].strip())
            elif line.startswith('x_values'):
                data['x_values'] = ast.literal_eval(line.split('=')[1].strip())
            elif line.startswith('Alphas'):
                data['Alphas'] = ast.literal_eval(line.split('=')[1].strip())
            elif line.startswith('Betas'):
                data['Betas'] = ast.literal_eval(line.split('=')[1].strip())
            elif line.startswith('Gammas'):
                data['Gammas'] = ast.literal_eval(line.split('=')[1].strip())
            elif line.startswith('Deltas'):
                data['Deltas'] = ast.literal_eval(line.split('=')[1].strip())
            elif line.startswith('DihedralAngles'):
                data['DihedralAngles'] = ast.literal_eval(line.split('=')[1].strip())

    return data

def sign(x):
    """
    Return '+' if x is positive, otherwise return '-'.
    
    Parameters:
        x (int or float): A number to check the sign.

    Returns:
        str: '+' if x > 0, '-' otherwise.
    """
    return "+" if x > 0 else "-"

def check_combinations(t_results, modulus_real, modulus_imag):
    """
    Check combinations of signs for t1, t2, t3, t4 and print the matching ones.
    
    Parameters:
        t_results (list): List of t-values (complex numbers).
        modulus_real (float): Modulus for real part check.
        modulus_imag (float): Modulus for imaginary part check.
    """
    # Iterate through all combinations of -1 and 1 for i, j, k
    for i in [-1, 1]:
        for j in [-1, 1]:
            for k in [-1, 1]:
                # Calculate the sum for real and imaginary parts
                real_sum = conditional_round(float(mp.re(t_results[0]) + i * mp.re(t_results[1]) + 
                            j * mp.re(t_results[2]) + k * mp.re(t_results[3]))) % modulus_real
                imag_sum = conditional_round(float(mp.im(t_results[0]) + i * mp.im(t_results[1]) + 
                            j * mp.im(t_results[2]) + k * mp.im(t_results[3]))) % modulus_imag
                
                # Check if both the real and imaginary sums are close to zero
                if np.isclose(real_sum, 0, atol=1e-5) and np.isclose(imag_sum, 0, atol=1e-5):
                    print(f"t1 {sign(i)} t2 {sign(j)} t3 {sign(k)} t4")

def double_check(Alphas, Betas, Gammas, Deltas, tolerance=1e-5):
    """
    Perform consistency checks on Alphas, Betas, Gammas, and Deltas.

    This function calculates the sigma values, then computes the sine ratios (A, B, C, D) for Alphas, Betas, Gammas, and Deltas.
    Using these ratios, it computes the r, s, and M values. The function checks if these values are consistent within a 
    specified tolerance. If all conditions are met, it returns True, otherwise False.

    Parameters:
        Alphas, Betas, Gammas, Deltas (list): Lists of angle values in degrees.
        tolerance (float): The tolerance within which the consistency is checked. Default is 1e-5.

    Returns:
        bool: True if all calculated values meet the consistency conditions, False otherwise.
    """
    
    # Step 1: Calculate sigma values as (alpha + beta + gamma + delta) / 2 for each corresponding set
    Sigmas = [(a + b + g + d) / 2 for a, b, g, d in zip(Alphas, Betas, Gammas, Deltas)]
    
    # Step 2: Calculate the sine ratios for Alphas, Betas, Gammas, Deltas
    As = [sin(v) / sin(sigma - v) for v, sigma in zip(Alphas, Sigmas)]
    Bs = [sin(v) / sin(sigma - v) for v, sigma in zip(Betas, Sigmas)]
    Cs = [sin(v) / sin(sigma - v) for v, sigma in zip(Gammas, Sigmas)]
    Ds = [sin(v) / sin(sigma - v) for v, sigma in zip(Deltas, Sigmas)]

    # Step 3: Compute r, s, and M values using the ratios
    r_values = [a * d for a, d in zip(As, Ds)]  # r1, r2, r3, r4
    s_values = [c * d for c, d in zip(Cs, Ds)]  # s1, s2, s3, s4
    M_values = [a * b * c * d for a, b, c, d in zip(As, Bs, Cs, Ds)]  # M1, M2, M3, M4

    # Step 4: Check if r, s, and M values are consistent within the specified tolerance
    consistent_r = np.isclose(r_values[0], r_values[1], atol=tolerance) and np.isclose(r_values[2], r_values[3], atol=tolerance)
    consistent_s = np.isclose(s_values[0], s_values[3], atol=tolerance) and np.isclose(s_values[1], s_values[2], atol=tolerance)
    consistent_M = np.isclose(np.min(M_values), np.max(M_values), atol=tolerance)

    # Return True if all conditions are met, otherwise False
    return consistent_r and consistent_s and consistent_M

def main():
    """
    Main function to read input data, compute t-values using elliptic parameters, and print the results.

    This function reads input values from a specified file, calculates intermediate parameters such as
    r, s, f, and M, computes elliptic integrals (K and K'), and then calculates t-values using these
    parameters. The results are printed along with a decomposition in terms of K and K'. Finally, 
    it checks the computed values against specific conditions and prints whether the conditions are satisfied.
    """
    
    # Read input values from the specified file (replace file_path with the actual file location)
    file_path = 'angles_2024_09_21_18_08_21.txt'  # Example file path, adjust as needed
    values = read_values_from_file(file_path)  # Function to read and return values from file

    # Unpack the values from the file
    x_values = values['x_values']  # The list of x-values
    Alphas = values['Alphas']  # List of Alpha angles
    Betas = values['Betas']  # List of Beta angles
    Gammas = values['Gammas']  # List of Gamma angles
    Deltas = values['Deltas']  # List of Delta angles
    DihedralAngles = values['DihedralAngles']  # List of Dihedral angles

    # Unpack individual x_values into variables for clarity
    x1, x2, x3, x4, x5, x6, x7, x8, x9 = x_values

    # Calculate the r, s, and f values based on x_values
    r_values = [1 + 1 / x1, 1 + 1 / x2]  # r1 = 1 + 1/x1, r2 = 1 + 1/x2
    s_values = [1 + 1 / x3, 1 + 1 / x4]  # s1 = 1 + 1/x3, s2 = 1 + 1/x4
    f_values = [1 + 1 / x5, 1 + 1 / x6, 1 + 1 / x7, 1 + 1 / x8]  # f1, f2, f3, f4 from x5 to x8
    M = 1 - x9  # Compute M as 1 - x9

    # Calculate the elliptic integrals K and Kprime based on M
    m, K, Kprime = calculate_elliptic_params(M)  # Function to calculate elliptic parameters

    # Compute t-values for each set of angles (Alpha, Beta, Gamma, Delta) using a loop
    t_results = [
        t_values(
            M, 
            Alphas[i], 
            Betas[i], 
            Gammas[i], 
            Deltas[i], 
            r_values[i // 2],  # Use r1 for i=0,1 and r2 for i=2,3
            s_values[0 if i in {0, 3} else 1],  # Use s1 for i=0,3 and s2 for i=1,2
            f_values[i]  # Corresponding f-value for each i
        ) 
        for i in range(4)  # Loop for the 4 sets of angles
    ]

    # Print the elliptic integrals K and Kprime
    print(f"K = {K}")
    print(f"K' = {Kprime}")

    # Print each t-value along with its decomposition in terms of K and K'
    for i, t in enumerate(t_results, 1):  # Start enumeration from 1 for t1, t2, etc.
        print(f"t{i} = {t} = {mp.re(t) / K}*K + {mp.im(t) / Kprime}*K'")  # Decompose t-values using K and K'

    # Based on the value of M, check and print the results
    if M > 1:
        # For M > 1, check the combinations with 2*K and 2*Kprime
        check_combinations(t_results, 2 * float(K), 2 * float(Kprime))
    else:
        # For M <= 1, check combinations with 4*K and 2*Kprime
        check_combinations(t_results, 4 * float(K), 2 * float(Kprime))
    
    # Perform a double-check to ensure calculated angles meet tolerance conditions
    if double_check(Alphas, Betas, Gammas, Deltas):
        print("Conditions satisfied: The calculated values meet the required tolerance.")
    else:
        print("Conditions NOT satisfied: The calculated values do not meet the required tolerance.")

if __name__ == "__main__":
    main()