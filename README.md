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
- [License](#license)

## Overview
This project solves the problem of generating flexible quadrilateral meshes of equimodular elliptic type by utilizing complex mathematical models and optimization techniques. The workflow is divided into several key steps: generating the normals for mesh faces, calculating sample angles, solving optimization problems, validating results through t-values, and transforming angles into vertices for 3D visualization.

The core of the project lies in transforming given angles and constraints into a flexible mesh structure with specific properties (non-self-intersecting, equal moduli, equal amplitudes at common vertices and sum of shifts to be a period). It ensures the quadrilateral mesh generated meets the required geometric conditions.

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
The first step in the workflow is generating the normals to the faces of the quadrilateral mesh. This is done using the `00_normals_generator.py` file. This script generates the normals to the faces of the mesh given a set of angles. If you already have a set of angles, you can skip this step and step 2, starting directly from step 3. However, this step is useful for creating sample normals to test step 2, which generates sample angles based on the given normals.

### Step 2: Sample Angles Calculation
Once the normals are given, the next step is to calculate the sample angles based on the provided normals. This is handled by `01_sample_angles_calculator.py`. This step calculates sample angles for the quadrilateral mesh based on given normals in different coordinate systems. These angles serve as an initial set for further optimization. 

### Step 3: Optimal Solution
In this step, the optimization problem is defined and solved using the Sequential Least SQuares Programming (SLSQP) method. The optimization problem minimizes the differences between theoretical angles and the actual angles of the quadrilateral mesh. The SLSQP algorithm is an iterative method for constrained nonlinear optimization, used to find optimal solutions for the given constraints and variables. It minimizes a scalar objective function subject to equality and inequality constraints. The algorithm ensures that the solution adheres to geometric constraints while finding the best angles for the mesh.

### Step 4: Optimal Angle Calculation
Once the optimal solution is found, the `03_optimal_angles.py` script is run to calculate the best set of angles based on the solution of the optimization problem. This step finalizes the calculated angles for the quadrilateral mesh, ensuring that they meet the constraints set by the optimization problem.

### Step 5: t-Values Calculation and Double Check
To ensure the validity of the generated quadrilateral mesh, `04_optimal_ts_and_double_check.py` calculates the necessary t-values and ensures that the mesh satisfies all required constraints: The script computes the t-values and verifies that the sum of the shifts is a period. It also checks if the quadrilateral mesh has equal moduli and amplitudes at common vertices.

### Step 6: Transforming Angles to Mesh Vertices
Finally, `05_angles_to_vertices.py` transforms the given angles into 3D coordinates of the quadrilateral mesh. It also avoids self-intersection cases and visualizes the mesh. This script takes in angles and a single edge length, computes the vertex coordinates, and exports the resulting mesh as an OBJ file for external visualization.


## License

This project is licensed under the Apache License 2.0. See the [LICENSE](./LICENSE) file for more details.

