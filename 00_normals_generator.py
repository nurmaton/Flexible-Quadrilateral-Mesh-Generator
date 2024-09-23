import numpy as np

def cos(x):
    """Convert angle from degrees to radians and return its cosine."""
    return np.cos(np.deg2rad(x, dtype = 'float64'))

def sin(x):
    """Convert angle from degrees to radians and return its sine."""
    return np.sin(np.deg2rad(x, dtype = 'float64'))

def calculate_betas():
    """
    Calculate the values of beta1, beta2, beta3, and beta4 using the provided angles.
    
    Returns:
    list: A list containing beta1, beta2, beta3, and beta4.
    """
    Betas = []
    
    # Calculate beta1
    beta1 = np.degrees(np.arccos(
        cos(alpha1) * cos(gamma1) * cos(delta1) +
        cos(psi1) * cos(alpha1) * sin(gamma1) * sin(delta1) +
        cos(phi) * sin(alpha1) * cos(gamma1) * sin(delta1) -
        cos(phi) * cos(psi1) * sin(alpha1) * sin(gamma1) * cos(delta1) +
        sin(phi) * sin(psi1) * sin(alpha1) * sin(gamma1)
    ))

    # Calculate beta2
    beta2 = np.degrees(np.arccos(
        cos(alpha2) * cos(gamma2) * cos(delta2) +
        cos(psi2) * cos(alpha2) * sin(gamma2) * sin(delta2) +
        cos(phi) * sin(alpha2) * cos(gamma2) * sin(delta2) -
        cos(phi) * cos(psi2) * sin(alpha2) * sin(gamma2) * cos(delta2) +
        sin(phi) * sin(psi2) * sin(alpha2) * sin(gamma2)
    ))

    # Calculate beta3
    beta3 = np.degrees(np.arccos(
        cos(alpha3) * cos(gamma3) * cos(delta3) +
        cos(psi2) * cos(alpha3) * sin(gamma3) * sin(delta3) +
        cos(theta) * sin(alpha3) * cos(gamma3) * sin(delta3) -
        cos(theta) * cos(psi2) * sin(alpha3) * sin(gamma3) * cos(delta3) +
        sin(theta) * sin(psi2) * sin(alpha3) * sin(gamma3)
    ))
    
    # Calculate beta4
    beta4 = np.degrees(np.arccos(
        cos(alpha4) * cos(gamma4) * cos(delta4) +
        cos(psi1) * cos(alpha4) * sin(gamma4) * sin(delta4) +
        cos(theta) * sin(alpha4) * cos(gamma4) * sin(delta4) -
        cos(theta) * cos(psi1) * sin(alpha4) * sin(gamma4) * cos(delta4) +
        sin(theta) * sin(psi1) * sin(alpha4) * sin(gamma4)
    ))

    Betas.extend([beta1, beta2, beta3, beta4])
    return Betas

def main():
    """
    Main function to compute normal vectors based on provided angles.

    This function calculates normal vectors at various points using the provided
    angles (Alphas, Betas, Gammas, Deltas, and DihedralAngles). The results are
    printed and saved to a file named 'normals.txt'.

    Angle Definitions:
        - alpha1, alpha2, alpha3, alpha4: Angles at points A1, A2, A3, and A4 respectively.
            - alpha1 is angle A2A1B1
            - alpha2 is angle A1A2B2
            - alpha3 is angle A4A3B3
            - alpha4 is angle A3A4B4

        - beta1, beta2, beta3, beta4: Angles at points B1, B2, B3, and B4 respectively.
            - beta1 is angle B1A1C1
            - beta2 is angle B2A2C2
            - beta3 is angle B3A3C3
            - beta4 is angle B4A4C4

        - gamma1, gamma2, gamma3, gamma4: Angles opposite to alphas at respective points.
            - gamma1 is angle A4A1C1
            - gamma2 is angle A3A2C2
            - gamma3 is angle A2A3C3
            - gamma4 is angle A1A4C4

        - delta1, delta2, delta3, delta4: Internal angles of quadrilateral A1A2A3A4.
            - delta1 is angle A4A1A2
            - delta2 is angle A1A2A3
            - delta3 is angle A2A3A4
            - delta4 is angle A3A4A1

    Dihedral Angle Definitions:
        - phi: Angle between faces A1A2A3A4 (face C, 11, 12, 13, or 14) and A1A2B2B1 (face 1 or C1).
        - psi2: Angle between faces A1A2A3A4 (face C, 11, 12, 13, or 14) and A2A3C3C2 (face 2 or C2).
        - theta: Angle between faces A1A2A3A4 (face C, 11, 12, 13, or 14) and A3A4B4B3 (face 3 or C3).
        - psi1: Angle between faces A1A2A3A4 (face C, 11, 12, 13, or 14) and A1A4C4C1 (face 4 or C4).

    Additional Angles:
        - psi11: Angle between faces C1 and A2B2C2 (face 41 or 22).
        - psi21: Angle between faces C1 and A1B1C1 (face 21 or 44).
        - psi12: Angle between faces C2 and A3B3C3 (face 42 or 23).
        - psi22: Angle between faces C2 and A2B2C2 (face 22 or 41).
        - psi13: Angle between faces C3 and A4B4C4 (face 43 or 24).
        - psi23: Angle between faces C3 and A3B3C3 (face 23 or 42).
        - psi14: Angle between faces C4 and A1B1C1 (face 44 or 21).
        - psi24: Angle between faces C4 and A4B4C4 (face 24 or 43).

    Note:
        The faces are numbered for clarity:
        - Face C corresponds to faces 11, 12, 13, or 14.
        - Faces 1 (or C1), 2 (or C2), 3 (or C3), and 4 (or C4) are adjacent to face C.
        - The numbers in parentheses indicate alternative face labels.
    """
    # Global values for angles
    global alpha1, alpha2, alpha3, alpha4
    global gamma1, gamma2, gamma3, gamma4
    global delta1, delta2, delta3, delta4
    global phi, psi2, theta, psi1
    
    # Set up values for Alphas, Betas, Gammas, Deltas, and DihedralAngles
    alpha1, alpha2, alpha3, alpha4 = [28.734040598418567, 151.26595940207199, 141.82734464226124, 38.17265535784382]
    gamma1, gamma2, gamma3, gamma4 = [159.76424407072537, 60.45797263708481, 96.65466960386112, 156.74078364249067]
    delta1, delta2, delta3, delta4 = [50, 130, 100, 80]
    phi, psi2, theta, psi1 = [125, 135, 140, 145]
    
    # Calculate the beta angles
    beta1, beta2, beta3, beta4 = calculate_betas()
    
    # Coordinate system A2xyz:
    # - Quadrilateral A1A2A3A4 lies in the xy-plane.
    # - x-axis is along A2A1.
    # - y-axis forms an acute angle with A2A4.

    # Normal vector to plane C in A2xyz coordinate system.
    NC = (0, 0, 1)

    # Normal vectors to faces 1, 2, 3, and 4 in A2xyz coordinate system.
    N1 = (0, -sin(phi), cos(phi))  # Face 1
    N2 = (-sin(delta2) * sin(psi2), cos(delta2) * sin(psi2), cos(psi2))  # Face 2
    N3 = (sin(delta2 + delta3) * sin(theta),
          -cos(delta2 + delta3) * sin(theta), cos(theta))  # Face 3
    N4 = (sin(delta1) * sin(psi1), cos(delta1) * sin(psi1), cos(psi1))  # Face 4

    # Coordinate system A1x1y1z1:
    # - Quadrilateral A1B1B2A2 lies in the x1y1-plane.
    # - x1-axis is along A1A2.
    # - y1-axis forms an acute angle with A1B1.

    # Compute dihedral angles psi11 and psi21.
    # psi11 is angle between faces C1 and A2B2C2 (face 41 or 22).
    # psi21 is angle between faces C1 and A1B1C1 (face 21 or 44).
    psi11_numerator = (
        sin(alpha2) * cos(gamma2) +
        sin(gamma2) * cos(alpha2) * cos(phi) * cos(psi2)
    ) * cos(delta2) - (
        sin(delta2) * cos(gamma2) * cos(phi) +
        sin(gamma2) * sin(phi) * sin(psi2)
    ) * cos(alpha2) + sin(alpha2) * sin(delta2) * sin(gamma2) * cos(psi2)
    psi11_denominator = sin(beta2)
    psi11 = np.degrees(np.arccos(psi11_numerator / psi11_denominator))

    psi21_numerator = (
        sin(alpha1) * cos(gamma1) +
        sin(gamma1) * cos(alpha1) * cos(phi) * cos(psi1)
    ) * cos(delta1) - (
        sin(delta1) * cos(gamma1) * cos(phi) +
        sin(gamma1) * sin(phi) * sin(psi1)
    ) * cos(alpha1) + sin(alpha1) * sin(delta1) * sin(gamma1) * cos(psi1)
    psi21_denominator = sin(beta1)
    psi21 = np.degrees(np.arccos(psi21_numerator / psi21_denominator))

    # Normal vector to plane C1 in A1x1y1z1 coordinate system.
    NC1 = (0, 0, 1)

    # Normal vectors to faces 21 and 41 in A1x1y1z1 coordinate system.
    N21 = (-sin(alpha1) * sin(psi21), cos(alpha1) * sin(psi21), cos(psi21))  # Face 21
    N41 = (sin(alpha2) * sin(psi11), cos(alpha2) * sin(psi11), cos(psi11))   # Face 41

    # Coordinate system A2x2y2z2:
    # - Quadrilateral A2C2C3A3 lies in the x2y2-plane.
    # - x2-axis is along A2A3.
    # - y2-axis forms an acute angle with A2C2.

    # Compute dihedral angles psi12 and psi22.
    # psi12 is angle between faces C2 and A3B3C3 (face 42 or 23).
    # psi22 is angle between faces C2 and A2B2C2 (face 22 or 41).
    psi12_numerator = (
        sin(delta3) * sin(gamma3) * cos(theta) -
        sin(psi2) * sin(theta) * cos(gamma3)
    ) * sin(alpha3) + (
        sin(alpha3) * cos(gamma3) * cos(psi2) * cos(theta) +
        sin(gamma3) * cos(alpha3)
    ) * cos(delta3) - sin(delta3) * cos(alpha3) * cos(gamma3) * cos(psi2)
    psi12_denominator = sin(beta3)
    psi12 = np.degrees(np.arccos(psi12_numerator / psi12_denominator))

    psi22_numerator = (
        sin(delta2) * sin(gamma2) * cos(phi) -
        sin(psi2) * sin(phi) * cos(gamma2)
    ) * sin(alpha2) + (
        sin(alpha2) * cos(gamma2) * cos(psi2) * cos(phi) +
        sin(gamma2) * cos(alpha2)
    ) * cos(delta2) - sin(delta2) * cos(alpha2) * cos(gamma2) * cos(psi2)
    psi22_denominator = sin(beta2)
    psi22 = np.degrees(np.arccos(psi22_numerator / psi22_denominator))

    # Normal vector to plane C2 in A2x2y2z2 coordinate system.
    NC2 = (0, 0, 1)

    # Normal vectors to faces 22 and 42 in A2x2y2z2 coordinate system.
    N22 = (-sin(gamma2) * sin(psi22), cos(gamma2) * sin(psi22), cos(psi22))  # Face 22
    N42 = (sin(gamma3) * sin(psi12), cos(gamma3) * sin(psi12), cos(psi12))   # Face 42

    # Coordinate system A3x3y3z3:
    # - Quadrilateral A3B3B4A4 lies in the x3y3-plane.
    # - x3-axis is along A3A4.
    # - y3-axis forms an acute angle with A3B3.

    # Compute dihedral angles psi13 and psi23.
    # psi13 is angle between faces C3 and A4B4C4 (face 43 or 24).
    # psi23 is angle between faces C3 and A3B3C3 (face 23 or 42).
    psi13_numerator = (
        sin(alpha4) * cos(gamma4) +
        sin(gamma4) * cos(alpha4) * cos(psi1) * cos(theta)
    ) * cos(delta4) - (
        sin(delta4) * cos(gamma4) * cos(theta) +
        sin(gamma4) * sin(psi1) * sin(theta)
    ) * cos(alpha4) + sin(alpha4) * sin(delta4) * sin(gamma4) * cos(psi1)
    psi13_denominator = sin(beta4)
    psi13 = np.degrees(np.arccos(psi13_numerator / psi13_denominator))

    psi23_numerator = (
        sin(alpha3) * cos(gamma3) +
        sin(gamma3) * cos(alpha3) * cos(psi2) * cos(theta)
    ) * cos(delta3) - (
        sin(delta3) * cos(gamma3) * cos(theta) +
        sin(gamma3) * sin(psi2) * sin(theta)
    ) * cos(alpha3) + sin(alpha3) * sin(delta3) * sin(gamma3) * cos(psi2)
    psi23_denominator = sin(beta3)
    psi23 = np.degrees(np.arccos(psi23_numerator / psi23_denominator))

    # Normal vector to plane C3 in A3x3y3z3 coordinate system.
    NC3 = (0, 0, 1)

    # Normal vectors to faces 23 and 43 in A3x3y3z3 coordinate system.
    N23 = (-sin(alpha3) * sin(psi23), cos(alpha3) * sin(psi23), cos(psi23))  # Face 23
    N43 = (sin(alpha4) * sin(psi13), cos(alpha4) * sin(psi13), cos(psi13))   # Face 43

    # Coordinate system A4x4y4z4:
    # - Quadrilateral A4C4C1A1 lies in the x4y4-plane.
    # - x4-axis is along A4A1.
    # - y4-axis forms an acute angle with A4C4.

    # Compute dihedral angles psi14 and psi24.
    # psi14 is angle between faces C4 and A1B1C1 (face 44 or 21).
    # psi24 is angle between faces C4 and A4B4C4 (face 24 or 43).
    psi14_numerator = (
        sin(delta1) * sin(gamma1) * cos(phi) -
        sin(psi1) * sin(phi) * cos(gamma1)
    ) * sin(alpha1) + (
        sin(alpha1) * cos(gamma1) * cos(psi1) * cos(phi) +
        sin(gamma1) * cos(alpha1)
    ) * cos(delta1) - sin(delta1) * cos(alpha1) * cos(gamma1) * cos(psi1)
    psi14_denominator = sin(beta1)
    psi14 = np.degrees(np.arccos(psi14_numerator / psi14_denominator))

    psi24_numerator = (
        sin(delta4) * sin(gamma4) * cos(theta) -
        sin(psi1) * sin(theta) * cos(gamma4)
    ) * sin(alpha4) + (
        sin(alpha4) * cos(gamma4) * cos(psi1) * cos(theta) +
        sin(gamma4) * cos(alpha4)
    ) * cos(delta4) - sin(delta4) * cos(alpha4) * cos(gamma4) * cos(psi1)
    psi24_denominator = sin(beta4)
    psi24 = np.degrees(np.arccos(psi24_numerator / psi24_denominator))

    # Normal vector to plane C4 in A4x4y4z4 coordinate system.
    NC4 = (0, 0, 1)

    # Normal vectors to faces 24 and 44 in A4x4y4z4 coordinate system.
    N24 = (-sin(gamma4) * sin(psi24), cos(gamma4) * sin(psi24), cos(psi24))  # Face 24
    N44 = (sin(gamma1) * sin(psi14), cos(gamma1) * sin(psi14), cos(psi14))   # Face 44

    # Prepare the output text.
    text = (
        f"A2xyz\n"
        f"NC = {NC}\n"
        f"N1 = {N1}\n"
        f"N2 = {N2}\n"
        f"N3 = {N3}\n"
        f"N4 = {N4}\n\n"
        f"A1x1y1z1\n"
        f"NC1 = {NC1}\n"
        f"N21 = {N21}\n"
        f"N41 = {N41}\n\n"
        f"A2x2y2z2\n"
        f"NC2 = {NC2}\n"
        f"N22 = {N22}\n"
        f"N42 = {N42}\n\n"
        f"A3x3y3z3\n"
        f"NC3 = {NC3}\n"
        f"N23 = {N23}\n"
        f"N43 = {N43}\n\n"
        f"A4x4y4z4\n"
        f"NC4 = {NC4}\n"
        f"N24 = {N24}\n"
        f"N44 = {N44}"
    )

    # Print the results to the console.
    print(text)

    # Write results to a file.
    with open("normals.txt", "w") as normalsFile:
        normalsFile.write(text)

if __name__ == "__main__":
    main()