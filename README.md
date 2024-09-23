# Flexible Quadrilateral Mesh Generator (3x3 Case)

<p align="center">
  <img src="./pics/01.svg" alt="Example of 3x3 mesh" width="800"/>
</p>

<p align="center"><strong>Figure 1: Example of 3x3 mesh.</strong></p>


A robust Python toolset for generating flexible quadrilateral meshes of equimodular elliptic type. This project is dedicated to solving optimization problems related to mesh generation and transformation, ensuring the creation of valid and non-self-intersecting quadrilateral meshes. The tools are designed to support complex mathematical computations, such as calculating normals, angles, and optimizing for specific mesh properties.

## Table of Contents
- [Overview](#overview)
- [Key Features](#key-features)
- [Installation](#installation)
- [Usage](#usage)
  - [Step 1: Normals Generation](#step-1-normals-generation)
  - [Step 2: Sample Angles Calculation](#step-2-sample-angles-calculation)
  - [Step 3: Optimal Solution](#step-3-optimal-solution)
  - [Step 4: Optimal Angle Calculation](#step-4-optimal-angle-calculation)
  - [Step 5: t-Values Calculation and Double Check](#step-5-t-values-calculation-and-double-check)
  - [Step 6: Transforming Angles to Mesh Vertices](#step-6-transforming-angles-to-mesh-vertices)
- [Acknowledgments](#acknowledgments)
- [License](#license)

## Overview
This project solves the problem of generating flexible quadrilateral meshes of equimodular elliptic type by utilizing complex mathematical models and optimization techniques. The workflow is divided into several key steps: generating the normals for mesh faces, calculating sample angles, solving optimization problems, validating results through t-values, and transforming angles into vertices for 3D visualization.

The core of the project lies in transforming given angles and constraints into a flexible mesh structure with specific properties (non-self-intersecting, equal moduli at common vertices, etc.). It ensures the quadrilateral mesh generated meets the required geometric conditions.

## Key Features
- **Polynomial System Generation:** Ability to generate 4 and 9 polynomial equations for mesh flexibility based on the given deltas and dihedral angles.
- **Optimization Solver:** Uses the Sequential Least SQuares Programming (SLSQP) method to solve minimization problems.
- **t-Values Computation:** Ensures that the quadrilateral mesh satisfies the periodicity conditions.
- **3D Mesh Visualization:** Transform computed angles into mesh vertices, outputting 3D OBJ files for external visualization.
- **Handles Self-intersections:** Avoids mesh self-intersections during vertex calculation.

## Installation
1. Clone the repository:

   ```bash
   git clone https://github.com/nurmaton/Flexible-Quadrilateral-Mesh-Generator.git
   cd Flexible-Quadrilateral-Mesh-Generator
   ```
2. Install the required dependencies
   ```bash
   pip install -r requirements.txt
   ```
   
## Usage

### Step 1: Normals Generation
The first step in the workflow is generating the normals to the faces of the quadrilateral mesh. This is done using the `00_normals_generator.py file`. This script generates normal vectors for mesh faces, which will later be used for further computations and validations.

### Step 2: Sample Angles Calculation
Once the normals are generated, the next step is to calculate the sample angles based on the provided normals. This is handled by `01_sample_angles_calculator.py`. This step calculates sample angles for the quadrilateral mesh based on given normals in different coordinate systems. These angles serve as an initial set for further optimization. 

### Step 3: Optimal Solution
In this step, the optimization problem is defined and solved using the Sequential Least SQuares Programming (SLSQP) method. The optimization problem minimizes the differences between theoretical angles and the actual angles of the quadrilateral mesh. The SLSQP algorithm is an iterative method for constrained nonlinear optimization, used to find optimal solutions for the given constraints and variables. This step involves 2000 iterations with random initial guesses to obtain the best solution.

### Step 4: Optimal Angle Calculation
Once the optimal solution is found, the `03_optimal_angles.py` script is run to calculate the best set of angles based on the solution of the optimization problem. This step finalizes the calculated angles for the quadrilateral mesh, ensuring that they meet the constraints set by the optimization problem.

### Step 5: t-Values Calculation and Double Check
To ensure the validity of the generated quadrilateral mesh, `04_optimal_ts_and_double_check.py` calculates the necessary t-values and ensures that the mesh satisfies all required constraints, such as equal moduli and correct periodicity. The script computes the t-values and verifies that the sum of the shifts is a period. It also checks if the quadrilateral mesh has equal moduli and amplitudes at common vertices.

### Step 6: Transforming Angles to Mesh Vertices
Finally, `05_angles_to_vertices.py` transforms the given angles into 3D coordinates of the quadrilateral mesh. It also avoids self-intersection cases and visualizes the mesh. This script takes in angles and a single edge length, computes the vertex coordinates, and exports the resulting mesh as an OBJ file for external visualization.


## Overview

This project is designed to generate flexible quadrilateral meshes of elliptic type with an emphasis on minimizing differences between theoretically computed and given angles. The project uses various optimization techniques to compute optimal angle sets that match the given geometry constraints.

### Steps Overview:

1. **Step 1: Generating Equations from Given Delta Angles**  
   For a quadrilateral mesh, given four delta angles (`delta1`, `delta2`, `delta3`, `delta4`), a system of 4 polynomial equations with 4 unknowns is formulated. These equations generate flexible mesh examples of elliptic type where we can choose the deltas.

2. **Step 2: Adding Dihedral Angles**  
   With the delta angles and dihedral angles (`phi`, `psi2`, `theta`, `psi1`), a system of 9 polynomial equations with 9 unknowns is developed. This step allows us to generate a specific type of flexible mesh where both the deltas and dihedral angles are chosen.

3. **Step 3: Minimizing the Differences Between Given and Computed Angles**  
   In this step, the goal is to minimize the differences between given angles (`alpha`, `beta`, `gamma`) and their corresponding computed values while respecting the constraints from the equations generated in step 2.

4. **Step 4: Handling Normals for Angle Computation**  
   Using computed normals to the faces of the quadrilateral mesh, we can calculate "sample" angles for the mesh. This is useful when the normals are given but the angles are not directly provided.

### File Descriptions:

- **00_normals_generator.py**  
  This script generates the normals to the faces of the mesh given a set of angles. It is useful for creating sample meshes to test the optimization process.

- **01_sample_angles_calculator.py**  
  Given normals in different coordinate systems, this script calculates angles for a 3x3 equimodular elliptic type flexible mesh. These angles are optimized to create a mesh that closely matches the given object while maintaining flexibility.

- **02_optimal_solver.py**  
  This script defines and solves the minimization problem using the `SLSQP` (Sequential Least Squares Quadratic Programming) optimization method. It uses random initial guesses and 2000 iterations to converge on an optimal solution that minimizes the differences between theoretical and given angles.

  **What is SLSQP?**  
  `SLSQP` is an optimization algorithm for solving nonlinear constrained optimization problems. It minimizes a scalar objective function subject to equality and inequality constraints. The algorithm ensures that the solution adheres to geometric constraints while finding the best angles for the mesh.

### Dependencies

- Python 3.x
- `numpy`
- `scipy`

- **03_optimal_angles.py**  
  After solving the minimization problem in the previous step, this script calculates the optimal angles for the quadrilateral mesh. The angles are calculated based on the solution obtained from the minimization process and are constrained by the equations defined in earlier steps.
  
  **Key Highlights:**  
  - The script ensures that the calculated angles (`alpha`, `beta`, `gamma`, and `delta`) match the required constraints from the minimization problem.
  - These optimal angles are the final output of the system and are printed to the console or saved to a file for further analysis.
  - The code also includes a mechanism to handle boundary cases or edge cases where the angles might result in invalid mesh configurations.

- **04_optimal_ts_and_double_check.py**  
  In this step, the script calculates the optimal values of the complex variables `t`, which are essential for the mesh construction, and ensures that the sum of the shifts is a period. Additionally, it performs a double-check to validate that the generated quadrilateral mesh has equal moduli and that the amplitudes at common vertices are equal.

  **Key Highlights:**
  - **t-values calculation:** These are derived from the solution and checked to ensure they meet the periodicity condition.
  - **Double-check mechanism:** The code ensures that the quadrilateral has equal moduli at all vertices, meaning the mesh satisfies the desired symmetry and flexibility conditions.
  - The script verifies that all amplitudes at common vertices are equal, guaranteeing geometric consistency in the quadrilateral mesh.
  
  This step is crucial as it provides the final validation of the constructed flexible mesh, ensuring it adheres to the requirements of the equimodular elliptic type mesh.

- **05_angles_to_vertices.py**  
  This script transforms the calculated angles and the length of one edge into the 3D coordinates of the flexible quadrilateral mesh. It also avoids self-intersection cases during this transformation, ensuring that the resulting mesh is valid and can be easily visualized.

  **Key Highlights:**
  - **Coordinate transformation:** The provided angles and edge length are used to compute the precise vertex coordinates for the mesh. This transformation ensures the quadrilateral's flexibility while avoiding self-intersections.
  - **Mesh visualization:** The script includes functions to visualize the generated flexible mesh, allowing for a better understanding of its structure.
  - **OBJ file creation:** It outputs the generated mesh as an OBJ file, which can be loaded into 3D modeling software for further inspection or manipulation.
  
  This step is crucial as it enables the final visualization and export of the flexible quadrilateral mesh, ensuring that the structure can be reviewed and analyzed in different formats.



Steps Overview
Step 1: Polynomial System for Deltas
In the first step, given four angles delta1, delta2, delta3, and delta4, a system of four polynomial equations with four unknowns is generated. This system is used to create examples of quadrilateral flexible meshes of the elliptic type. At this stage, we can only choose the delta angles to define the mesh structure.

Step 2: Polynomial System for Deltas and Dihedral Angles
Building upon Step 1, in Step 2, both the delta angles (delta1, delta2, delta3, delta4) and the dihedral angles (phi, psi2, theta, psi1) are used. This generates a system of nine polynomial equations with nine unknowns, allowing for the generation of flexible quadrilateral meshes of the elliptic type with more control. Here, we can choose both the delta and dihedral angles, giving us the ability to create the type of mesh we desire.

Step 3: Optimization Problem for Angle Matching (Pending)
In this step, the optimization problem begins. Given a full set of angles—delta, alpha, beta, gamma—as well as dihedral angles, the objective is to minimize the difference between the given angles and those calculated theoretically, while respecting the constraints imposed by the system of nine equations from Step 2.

Step 4: Normal-Based Sample Angle Generation (Pending)
In this final step, formulas are derived to compute unit normals to the faces of the quadrilateral mesh in various coordinate systems. Using these normals, sample angles of the mesh can be obtained. This is useful when only normals are provided and not the angles directly.

Code Structure
The project contains several Python scripts, each addressing different parts of the quadrilateral mesh generation and optimization process:

00_normals_generator.py
This script generates the normals to the faces of the mesh based on the given angles.
It serves as a tool for quickly generating examples and testing other parts of the code.
01_sample_angles_calculator.py
This script calculates the angles of the object to be optimized.
The goal is to find a 3x3 equimodular elliptic-type flexible mesh that closely matches a given object, based on the provided normals in different coordinate systems.



## 02_optimal_solver.py

This script defines and solves an optimization problem that seeks to minimize the difference between theoretical and given angles (`alpha`, `beta`, `gamma`) subject to geometric constraints of the flexible quadrilateral mesh.

### Major Steps:
1. **Trigonometric Calculations**: Helper functions calculate the cosine, sine, and tangent of angles in degrees, which is essential as all input angles are specified in degrees.
   
2. **Coefficient Calculation**: 
   - The function `coefficients()` computes specific coefficients based on input angles (`delta`) and parameters (`x`, `y`).
   - These coefficients represent the geometric properties and constraints of the mesh, contributing to the formulation of the optimization problem.

3. **Constraint Generation**:
   - Constraints are generated to ensure that the solution remains geometrically valid and that the mesh retains flexibility. These constraints are non-linear and are crucial for ensuring the quadrilateral mesh meets the requirements.

4. **Optimization Using SLSQP**:
   - The script uses the `SLSQP` method (`Sequential Least Squares Quadratic Programming`) from the `scipy.optimize` library to solve the constrained non-linear optimization problem.
   - `SLSQP` is particularly well-suited for this problem because it handles both equality and inequality constraints, ensuring the solution respects geometric boundaries.

5. **Random Initial Guess**:
   - The optimization process starts with a random guess for the solution, which helps avoid local minima and increases the chances of finding an optimal solution.
   - The algorithm is set to run for up to 2000 iterations, maximizing the chances of convergence.

6. **Minimization of the Objective Function**:
   - The goal is to minimize the differences between the given angles and the computed angles based on theoretical models.
   - The solution is determined by balancing these differences while satisfying all the geometric constraints.

### About SLSQP:
`SLSQP` (Sequential Least Squares Quadratic Programming) is a powerful optimization algorithm for solving non-linear constrained optimization problems. It minimizes a scalar objective function while adhering to both equality and inequality constraints. In this case, it ensures that the mesh generation process results in a flexible quadrilateral mesh with the desired properties.

### Key Details:
- **Method**: `SLSQP`
- **Iterations**: 2000 (random initial guesses are used)
- **Purpose**: Solve for angles and parameters that minimize discrepancies in mesh geometry.



## License

This project is licensed under the Apache License 2.0. See the [LICENSE](./LICENSE) file for more details.

