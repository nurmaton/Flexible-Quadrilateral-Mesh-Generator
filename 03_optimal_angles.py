import numpy as np
import time

def cos(x):
    """Convert angle from degrees to radians and return its cosine."""
    return np.cos(np.deg2rad(x, dtype = 'float64'))

def sin(x):
    """Convert angle from degrees to radians and return its sine."""
    return np.sin(np.deg2rad(x, dtype = 'float64'))

def euclidean_distance(list_a, list_b):
    """Calculate the Euclidean distance between two lists of equal length."""
    return np.sqrt(np.sum((np.array(list_a) - np.array(list_b)) ** 2))

def find_closest_list(lists, target_list):
    """Find and return the list in 'lists' that is closest to the 'target_list' based on Euclidean distance."""
    distances = [euclidean_distance(lst, target_list) for lst in lists]
    closest_index = np.argmin(distances)
    return lists[closest_index]

def adjust_angle(angle):
    """Ensure the angle is within the [0, 180] degree range."""
    return angle if angle > 0 else 180 + angle

def compute_angles(g, sigma):
    """
    Compute the angles (alpha, beta, gamma, delta) based on 'g' values and 'sigma'.
    
    Args:
        g (list): List of four g-values.
        sigma (float): Angle in degrees.

    Returns:
        Tuple (alpha, beta, gamma, delta) in degrees.
    """
    alpha = adjust_angle(np.degrees(np.arctan(g[0] * sin(sigma) / (1 + g[0] * cos(sigma)))))
    beta = adjust_angle(np.degrees(np.arctan(g[1] * sin(sigma) / (1 + g[1] * cos(sigma)))))
    gamma = adjust_angle(np.degrees(np.arctan(g[2] * sin(sigma) / (1 + g[2] * cos(sigma)))))
    delta = adjust_angle(np.degrees(np.arctan(g[3] * sin(sigma) / (1 + g[3] * cos(sigma)))))
    return alpha, beta, gamma, delta

def process_sigma(cosSigma, g, deltaS):
    """
    Process and validate sigma to compute angles (alpha, beta, gamma, delta), returning valid angle sets.

    Args:
        cosSigma (float): Cosine of the sigma angle.
        g (list): List of four g-values.
        deltaS (float): Target delta angle for validation.

    Returns:
        Bag (list): List of valid angle sets [alpha, beta, gamma, delta].
    """
    Bag = []
    for sigma in [np.degrees(np.arccos(cosSigma)), np.degrees(-np.arccos(cosSigma) + 2 * np.pi)]:
        alpha, beta, gamma, delta = compute_angles(g, sigma)
        if np.isclose((alpha + beta + gamma + delta) / 2, sigma, atol=1e-5) and np.isclose(delta, deltaS, atol=1e-5):
            Bag.append([alpha, beta, gamma, delta])
    return Bag

def compute_g_values(r, f, s, M):
    """
    Compute the g values for the given r, f, s, and M.
    
    Args:
        r, f, s, M (float): Parameters for g-value computation.

    Returns:
        g_values (list): List of four g-values.
    """
    g1 = np.sqrt(r * f / s)
    g2 = M / (s * g1)
    g3 = f / g1
    g4 = (s / f) * g1
    return [g1, g2, g3, g4]

def compute_cosSigma(M, f, r, s, g1):
    """
    Compute the cosSigma value based on the provided inputs.
    
    Args:
        M, f, r, s, g1 (float): Parameters for cosSigma computation.

    Returns:
        cosSigma (float): The cosine of the sigma angle.
    """
    return (M * (f + r + s - 1) + r * f * s - r * f - r * s - f * s) / (2 * g1 * s * (1 - M))

def process_g_values_and_sigma(r, f, s, M, deltaS):
    """
    Compute g values and process sigma for both positive and negative square root cases.

    Args:
        r, f, s, M (float): Parameters for g-value computation.
        deltaS (float): Target delta angle for validation.

    Returns:
        Bag (list): List of valid angle sets [alpha, beta, gamma, delta].
    """
    Bag = []
    for sign in [1, -1]:
        g_values = compute_g_values(r, f, s, M)
        g_values = [sign * g for g in g_values]  # Apply positive and negative sqrt cases
        cosSigma = compute_cosSigma(M, f, r, s, g_values[0])
        Bag.extend(process_sigma(cosSigma, g_values, deltaS))
    return Bag

def read_x_values_from_file(file_path):
    """
    Read the first set of x1 to x9 values from the specified text file.

    Args:
        file_path (str): Path to the input file containing x-values.

    Returns:
        x_values (list): List of 9 float values for x1 through x9.
        switch (float): The switch value found in the file.
    """
    x_values = []
    switch = None  # Initialize switch to ensure it always has a value

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('S'):  # Check for the 'switch' line
                switch = int(line.split('=')[1].strip())
            if line.startswith('x'):  # Read lines starting with x1 to x9
                value = float(line.split('=')[1].strip())
                x_values.append(value)
                if len(x_values) == 9:  # Stop after reading the first 9 x-values
                    break

    # Ensure switch is not None (in case 's' line is missing)
    if switch is None:
        raise ValueError("No 'Switch' line found in the file.")

    return x_values, switch

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
    Main function to read input values, process them, and print the best angle sets.

    This function loads predefined sample angles (Alphas, Betas, Gammas, Deltas, and DihedralAngles),
    reads the x_values from a file, computes r, s, f, and M values, processes the data, and 
    determines the best angle sets that match the provided inputs.

    Results are printed and saved to a timestamped file.
    """

    # Load sample values for Alphas, Betas, Gammas, Deltas, and DihedralAngles from the file
    AlphasS, BetasS, GammasS, DeltasS, DihedralAnglesS = read_sample_angles_from_file("sample_angles.txt", 1)

    # Load x_values and switch from a separate file
    file_path = 'roots_2024_09_21_01_26_18.txt'  # Update this path with the actual file path if needed
    x_values, switch = read_x_values_from_file(file_path)

    # Compute r, s, f values and M based on the x_values loaded from the file
    r_values = [1 + 1 / x_values[0], 1 + 1 / x_values[1]]  # r_values[0] and r_values[1] from x_values
    s_values = [1 + 1 / x_values[2], 1 + 1 / x_values[3]]  # s_values[0] and s_values[1] from x_values
    f_values = [1 + 1 / x_values[4], 1 + 1 / x_values[5], 1 + 1 / x_values[6], 1 + 1 / x_values[7]]  # f_values from x_values
    M = 1 - x_values[8]  # Compute M as 1 - x_values[8]

    # Initialize a list to hold the best angle sets found
    bestBag = []
    
    # Process all cases, matching each set of angles with computed r, s, f, and delta values
    for i, (r, f, s, alphaS, betaS, gammaS, deltaS) in enumerate(zip([r_values[0], r_values[0], r_values[1], r_values[1]], 
                                                                      f_values, 
                                                                      [s_values[0], s_values[1], s_values[1], s_values[0]], 
                                                                      AlphasS, BetasS, GammasS, DeltasS), 1):
        # Process g values and sigma for each set, using the r, f, s, and M values
        Bag = process_g_values_and_sigma(r, f, s, M, deltaS)
        
        # If valid angle sets are found, find the closest set to the provided sample angles
        if Bag:
            bestABGD = find_closest_list(Bag, [alphaS, betaS, gammaS, deltaS])
            bestBag.append(bestABGD)
    
    # Unpack the best angles from the processed data into separate lists
    bestA, bestB, bestG, bestD = zip(*bestBag)
    
    # Convert the tuples into lists for easier manipulation
    bestA, bestB, bestG, bestD = list(bestA), list(bestB), list(bestG), list(bestD)

    # Print the best angles found for Alphas, Betas, Gammas, and Deltas
    print(f"Best Alphas: {bestA}\nBest Betas: {bestB}\nBest Gammas: {bestG}\nBest Deltas: {bestD}")
    
    # Get the current time to create a timestamped file name for saving the results
    now = time.localtime()
    current_time = time.strftime("%Y_%m_%d_%H_%M_%S", now)
    
    # Save the best angles and other data to a file with a timestamped name
    with open(f"optimal_angles_{current_time}.txt", "a") as optimalAnglesFile:
        optimalAnglesFile.write(f"Switch = {switch}\n\nx_values = {x_values}\n\nAlphas = {bestA}\n\nBetas = {bestB}\n\nGammas = {bestG}\n\nDeltas = {bestD}\n\nDihedralAngles = {DihedralAnglesS}")

if __name__ == "__main__":
    main()