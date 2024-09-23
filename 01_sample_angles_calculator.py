import numpy as np

def cos(x):
    """Convert angle from degrees to radians and return its cosine."""
    return np.cos(np.deg2rad(x, dtype='float64'))

def sin(x):
    """Convert angle from degrees to radians and return its sine."""
    return np.sin(np.deg2rad(x, dtype='float64'))

def parse_file(file_path):
    """
    Parses the file to extract normal vectors and returns them as numpy arrays.

    The file is expected to contain lines in the format:
    NC = (x, y, z)
    N1 = (x, y, z)
    ...
    The method reads the normal vectors and converts them into numpy arrays.

    Parameters:
    file_path (str): Path to the file containing the normal vectors.

    Returns:
    tuple: A tuple containing the normal vectors as numpy arrays in the following order:
    (NC, N1, N2, N3, N4, NC1, N21, N41, NC2, N22, N42, NC3, N23, N43, NC4, N24, N44)
    """
    normals = []
    with open(file_path, 'r') as file:
        for line in file:
            # Only process lines that contain '=' and ignore others
            if '=' in line:
                # Extract the part after the '=' and remove unwanted characters
                vector_str = line.split('=')[1].strip().strip('()')
                # Convert the string to a list of floats
                vector = [float(num) for num in vector_str.split(',')]
                normals.append(np.array(vector, dtype='float64'))
    
    return tuple(normals)

def normalize_vectors(*vectors):
    """
    Normalize a list of vectors.

    Parameters:
    vectors (tuple): A tuple of vectors to be normalized.

    Returns:
    list: A list of normalized vectors.
    """
    return [vec / np.linalg.norm(vec) for vec in vectors]

def calculate_dihedral_angles(NC, N1, N2, N3, N4, NC1, N21, N41, NC2, N22, N42, NC3, N23, N43, NC4, N24, N44):
    """
    Calculate the dihedral angles between NC and other normal vectors.

    Parameters:
    NC (numpy array): Normal vector NC.
    N1, N2, N3, N4 (numpy arrays): Normal vectors N1, N2, N3, N4.
    NC1, N21, N41, NC2, N22, N42, NC3, N23, N43, NC4, N24, N44 (numpy arrays): 
    Normal vectors for various faces in different coordinate systems.

    Returns:
    list: A list of calculated dihedral angles:
    [phi, psi2, theta, psi1, psi11, psi21, psi12, psi22, psi13, psi23, psi14, psi24]
    """
    # Calculate primary dihedral angles
    phi = np.degrees(np.arccos(np.dot(N1, NC)))
    psi2 = np.degrees(np.arccos(np.dot(N2, NC)))
    theta = np.degrees(np.arccos(np.dot(N3, NC)))
    psi1 = np.degrees(np.arccos(np.dot(N4, NC)))
    
    # Calculate the remaining dihedral angles using other coordinate systems
    psi11 = np.degrees(np.arccos(np.dot(N41, NC1)))
    psi21 = np.degrees(np.arccos(np.dot(N21, NC1)))
    psi12 = np.degrees(np.arccos(np.dot(N42, NC2)))
    psi22 = np.degrees(np.arccos(np.dot(N22, NC2)))
    psi13 = np.degrees(np.arccos(np.dot(N43, NC3)))
    psi23 = np.degrees(np.arccos(np.dot(N23, NC3)))
    psi14 = np.degrees(np.arccos(np.dot(N44, NC4)))
    psi24 = np.degrees(np.arccos(np.dot(N24, NC4)))
    
    DihedralAngles = [phi, psi2, theta, psi1, psi11, psi21, psi12, psi22, psi13, psi23, psi14, psi24]
    
    # Check if any angle is 0 or 180 (which is not allowed)
    if 0 in DihedralAngles[0:4] or 180 in DihedralAngles[0:4]:
        raise ValueError("Dihedral angle should not be 0 or 180 degrees.")
    
    return DihedralAngles

def calculate_angles(NC, N1, N2, N3, N4, NC1, N21, N41, NC2, N22, N42, NC3, N23, N43, NC4, N24, N44, phi, psi2, theta, psi1, psi11, psi21, psi12, psi22, psi13, psi23, psi14, psi24):
    """
    Calculate the delta, alpha, beta, and gamma angles for the system based on the provided normal vectors and dihedral angles.

    Parameters:
    NC, N1, N2, N3, N4, NC1, N21, N41, NC2, N22, N42, NC3, N23, N43, NC4, N24, N44 (numpy arrays):
        The normal vectors for the corresponding planes.

    phi, psi1, psi2, theta, psi11, psi21, psi12, psi22, psi13, psi23, psi14, psi24 (float):
        The dihedral angles between the respective normal vectors.

    Returns:
    list: A list of calculated angles in the format [Alphas, Betas, Gammas, Deltas],
          where Alphas, Betas, Gammas, and Deltas are the corresponding angles for each system.
    """
    
    # Initialize lists to store the alpha, beta, gamma, and delta angles.
    Alphas = []
    Betas = []
    Gammas = []
    Deltas = []
    output = []
    
    # Calculate delta1 based on N4 and psi1.
    delta1 = np.degrees(np.arccos(N4[1] / sin(psi1)))
    if not np.isclose(N4[0], sin(delta1) * sin(psi1), atol=1e-5):
        raise ValueError("Inconsistency detected in N4 calculations.")
    
    # Calculate delta2 based on N2 and psi2.
    delta2 = np.degrees(np.arccos(N2[1] / sin(psi2)))
    if not np.isclose(N2[0], -sin(delta2) * sin(psi2), atol=1e-5):
        raise ValueError("Inconsistency detected in N2 calculations.")
    
    # Calculate delta23 and determine if further correction is needed for N3.
    delta23 = np.degrees(np.arccos(-N3[1] / sin(theta)))
    parameter23 = False
    if not np.isclose(N3[0], sin(delta23) * sin(theta), atol=1e-5):
        parameter23 = True
    
    # Adjust delta3 and delta4 based on delta23 and delta2.
    if delta23 > delta2:
        delta3 = delta23 - delta2
        delta4 = 360 - delta1 - delta2 - delta3
        Deltas.append([delta1, delta2, delta3, delta4])
        delta23 = np.degrees(-np.arccos(-N3[1] / sin(theta)) + 2 * np.pi)
        if 0 < delta23 - delta2 < 180:
            delta3 = delta23 - delta2
            delta4 = 360 - delta1 - delta2 - delta3
            Deltas.append([delta1, delta2, delta3, delta4])
    else:
        delta23 = np.degrees(-np.arccos(-N3[1] / sin(theta)) + 2 * np.pi)
        if 0 < delta23 - delta2 < 180:
            delta3 = delta23 - delta2
            delta4 = 360 - delta1 - delta2 - delta3
            Deltas.append([delta1, delta2, delta3, delta4])
        if not np.isclose(N3[0], sin(delta23) * sin(theta), atol=1e-5):
            if parameter23:
                raise ValueError("Inconsistency detected in N3 calculations.")
    
    # Calculate alpha1 based on N21 and psi21.
    alpha1 = np.degrees(np.arccos(N21[1] / sin(psi21)))
    if not np.isclose(N21[0], -sin(alpha1) * sin(psi21), atol=1e-5):
        raise ValueError("Inconsistency detected in N21 calculations.")
    
    # Calculate alpha2 based on N41 and psi11.
    alpha2 = np.degrees(np.arccos(N41[1] / sin(psi11)))
    if not np.isclose(N41[0], sin(alpha2) * sin(psi11), atol=1e-5):
        raise ValueError("Inconsistency detected in N41 calculations.")
    
    # Calculate alpha3 based on N23 and psi23.
    alpha3 = np.degrees(np.arccos(N23[1] / sin(psi23)))
    if not np.isclose(N23[0], -sin(alpha3) * sin(psi23), atol=1e-5):
        raise ValueError("Inconsistency detected in N23 calculations.")
    
    # Calculate alpha4 based on N43 and psi13.
    alpha4 = np.degrees(np.arccos(N43[1] / sin(psi13)))
    if not np.isclose(N43[0], sin(alpha4) * sin(psi13), atol=1e-5):
        raise ValueError("Inconsistency detected in N43 calculations.")
    
    # Store calculated alpha angles.
    Alphas.extend([alpha1, alpha2, alpha3, alpha4])
    
    # Calculate gamma2 based on N22 and psi22.
    gamma2 = np.degrees(np.arccos(N22[1] / sin(psi22)))
    if not np.isclose(N22[0], -sin(gamma2) * sin(psi22), atol=1e-5):
        raise ValueError("Inconsistency detected in N22 calculations.")
    
    # Calculate gamma3 based on N42 and psi12.
    gamma3 = np.degrees(np.arccos(N42[1] / sin(psi12)))
    if not np.isclose(N42[0], sin(gamma3) * sin(psi12), atol=1e-5):
        raise ValueError("Inconsistency detected in N42 calculations.")
    
    # Calculate gamma1 based on N44 and psi14.
    gamma1 = np.degrees(np.arccos(N44[1] / sin(psi14)))
    if not np.isclose(N44[0], sin(gamma1) * sin(psi14), atol=1e-5):
        raise ValueError("Inconsistency detected in N44 calculations.")
    
    # Calculate gamma4 based on N24 and psi24.
    gamma4 = np.degrees(np.arccos(N24[1] / sin(psi24)))
    if not np.isclose(N24[0], -sin(gamma4) * sin(psi24), atol=1e-5):
        raise ValueError("Inconsistency detected in N24 calculations.")
    
    # Store calculated gamma angles.
    Gammas.extend([gamma1, gamma2, gamma3, gamma4])
    
    # Calculate beta1, beta2, beta3, and beta4 using the provided formulas (details in documentation).
    for deltas in Deltas:
        delta1, delta2, delta3, delta4 = deltas
        beta1 = np.degrees(np.arccos(
            cos(alpha1) * cos(gamma1) * cos(delta1) +
            cos(psi1) * cos(alpha1) * sin(gamma1) * sin(delta1) +
            cos(phi) * sin(alpha1) * cos(gamma1) * sin(delta1) -
            cos(phi) * cos(psi1) * sin(alpha1) * sin(gamma1) * cos(delta1) +
            sin(phi) * sin(psi1) * sin(alpha1) * sin(gamma1)
        ))

        beta2 = np.degrees(np.arccos(
            cos(alpha2) * cos(gamma2) * cos(delta2) +
            cos(psi2) * cos(alpha2) * sin(gamma2) * sin(delta2) +
            cos(phi) * sin(alpha2) * cos(gamma2) * sin(delta2) -
            cos(phi) * cos(psi2) * sin(alpha2) * sin(gamma2) * cos(delta2) +
            sin(phi) * sin(psi2) * sin(alpha2) * sin(gamma2)
        ))

        beta3 = np.degrees(np.arccos(
            cos(alpha3) * cos(gamma3) * cos(delta3) +
            cos(psi2) * cos(alpha3) * sin(gamma3) * sin(delta3) +
            cos(theta) * sin(alpha3) * cos(gamma3) * sin(delta3) -
            cos(theta) * cos(psi2) * sin(alpha3) * sin(gamma3) * cos(delta3) +
            sin(theta) * sin(psi2) * sin(alpha3) * sin(gamma3)
        ))
        
        beta4 = np.degrees(np.arccos(
            cos(alpha4) * cos(gamma4) * cos(delta4) +
            cos(psi1) * cos(alpha4) * sin(gamma4) * sin(delta4) +
            cos(theta) * sin(alpha4) * cos(gamma4) * sin(delta4) -
            cos(theta) * cos(psi1) * sin(alpha4) * sin(gamma4) * cos(delta4) +
            sin(theta) * sin(psi1) * sin(alpha4) * sin(gamma4)
        ))

        # Store calculated beta angles.
        Betas.append([beta1, beta2, beta3, beta4])
    
    # Combine calculated angles and store in the output list.
    for i in range(len(Deltas)):
        output.append([Alphas, Betas[i], Gammas, Deltas[i]])
    
    return output

def main():
    """
    Main function to load normal vectors from a file, normalize them, check for consistency,
    calculate dihedral and delta angles, and ensure they follow the system's geometric constraints.

    Steps:
    1. Load normal vectors from the specified file.
    2. Normalize the loaded normal vectors to ensure unit length.
    3. Verify that the normal vectors for the faces (NC, NC1, NC2, NC3, NC4) align with their corresponding z-axis.
    4. Calculate the dihedral angles between the faces.
    5. Calculate the delta angles based on the dihedral angles.
    6. Print and save the calculated angles to a file.

    Output:
    - The calculated options (Alphas, Betas, Gammas, Deltas, Dihedral Angles) are saved to 'sample_angles.txt'.
    """
    
    # Step 1: Load the normal vectors from a file
    file_path = 'normals.txt'  # Update this path to the actual file location.
    NC, N1, N2, N3, N4, NC1, N21, N41, NC2, N22, N42, NC3, N23, N43, NC4, N24, N44 = parse_file(file_path)

    # Step 2: Normalize the normal vectors to ensure unit length
    NC, N1, N2, N3, N4, NC1, N21, N41, NC2, N22, N42, NC3, N23, N43, NC4, N24, N44 = normalize_vectors(
        NC, N1, N2, N3, N4, NC1, N21, N41, NC2, N22, N42, NC3, N23, N43, NC4, N24, N44
    )

    # Step 3: Verify that the face normals (NC, NC1, NC2, NC3, NC4) align with the z-axis
    for normal, label in zip([NC, NC1, NC2, NC3, NC4], ['NC', 'NC1', 'NC2', 'NC3', 'NC4']):
        if not np.allclose(normal, [0, 0, 1], atol=1e-6):
            raise ValueError(f"Inconsistency detected in {label}: {normal}")

    # Step 4: Calculate dihedral angles between the faces
    DihedralAngles = calculate_dihedral_angles(
        NC, N1, N2, N3, N4, NC1, N21, N41, NC2, N22, N42, NC3, N23, N43, NC4, N24, N44
    )
    
    # Step 5: Extract the dihedral angles needed for delta calculations
    phi, psi2, theta, psi1, psi11, psi21, psi12, psi22, psi13, psi23, psi14, psi24 = DihedralAngles

    # Step 6: Calculate delta angles based on the dihedral angles and normal vectors
    Angles = calculate_angles(
        NC, N1, N2, N3, N4, NC1, N21, N41, NC2, N22, N42, NC3, N23, N43, NC4, N24, N44,
        phi, psi2, theta, psi1, psi11, psi21, psi12, psi22, psi13, psi23, psi14, psi24
    )

    # Step 7: Prepare output text for each angle configuration
    text = ""
    for i, element in enumerate(Angles, start=1):
        text += (f"Option {i}\n\n"
                 f"Alphas = {element[0]}\n"
                 f"Betas = {element[1]}\n"
                 f"Gammas = {element[2]}\n"
                 f"Deltas = {element[3]}\n"
                 f"DihedralAngles = {DihedralAngles[0:4]}\n\n\n")

    # Step 8: Write the results to a file
    with open("sample_angles.txt", "w") as sampleAnglesFile:
        sampleAnglesFile.write(text)

if __name__ == "__main__":
    main()
