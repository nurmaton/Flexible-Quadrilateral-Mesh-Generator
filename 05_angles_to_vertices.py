import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from mpl_toolkits.mplot3d import Axes3D, proj3d
import pylab, ast

def cos(x):
    """Convert angle from degrees to radians and return its cosine."""
    return np.cos(np.deg2rad(x, dtype = 'float64'))

def sin(x):
    """Convert angle from degrees to radians and return its sine."""
    return np.sin(np.deg2rad(x, dtype = 'float64'))

def generate_edge(A1, angle1, angle2, angle3, booster=1):
    """
    Generates the edge length A2 based on the provided angles and A1, applying conditions
    from a quadrilateral configuration.

    The method uses trigonometric relationships between angles to calculate A2, ensuring
    that the quadrilateral satisfies certain geometric constraints derived from a3 > 0 and a4 > 0.

    Parameters:
    A1 (float): Length of the first edge.
    angle1 (float): The angle between A1 and another edge (in degrees).
    angle2 (float): The second angle between the edges (in degrees).
    angle3 (float): The third angle between the edges (in degrees).
    booster (float, optional): A scaling factor used to boost the value of A2. Default is 1.

    Returns:
    float: The calculated edge length A2.
    """
    # Conditions derived from a3 > 0 and a4 > 0 (see reference PDF for details).
    # For the quadrilateral to be valid with angle1, angle2, angle3, A1, A2 are taken 
    # in the same order as delta1, delta2, delta3 (a1=A1A2, a2=A2A3) corresponding 
    # to the quadrilateral A1A2A3A4.
    
    if angle2 + angle3 >= 180:
        # Case when the sum of angle2 and angle3 is larger than or equal to 180 degrees
        if angle1 + angle2 >= 180:
            # If angle1 + angle2 is also larger than or equal to 180 degrees
            A2 = booster * A1
            # Checking and debugging output when booster is 1000
            if booster == 1000:
                print("A2 = A1")
        else:
            # Calculate A2 based on the sine rule when angle1 + angle2 is less than 180
            A2 = booster * 0.8 * A1 * sin(angle1) / sin(angle1 + angle2)
            # Checking and debugging output when booster is 1000
            if booster == 1000:
                print(0.8)
    else:
        # Case when angle2 + angle3 is less than 180 degrees
        if angle1 + angle2 >= 180:
            # Calculate A2 based on the sine rule when angle1 + angle2 is greater than or equal to 180
            A2 = booster * 1.2 * A1 * sin(angle2 + angle3) / sin(angle3)
            # Checking and debugging output when booster is 1000
            if booster == 1000:
                print(1.2)
        else:
            # Let's compare sin(angle1) / sin(angle1 + angle2) and sin(angle2 + angle3) / sin(angle3)
            leftScale = sin(angle2 + angle3) / sin(angle3)
            rightScale = sin(angle1) / sin(angle1 + angle2)
            
            # Determine scaleRatio based on the comparison of the two sine ratios
            if leftScale > rightScale:
                scaleRatio = leftScale // rightScale + 1
            elif rightScale > leftScale:
                scaleRatio = rightScale // leftScale + 1
            else:
                scaleRatio = 2

            # Checking and debugging output when booster is 1000
            if booster == 1000:
                print("scaleRatio:", scaleRatio)
            
            # Compute A2 using a weighted average of the two scales based on booster and scaleRatio
            A2 = A1 * ((1 - booster/scaleRatio) * leftScale + booster/scaleRatio * rightScale)
    
    return A2

def generate_tilde_angle(angle1, angle2, booster=1):
    """
    Generates a tilde angle based on the input angles (angle1, angle2) and a booster value.

    The method computes the tilde angle using the formula:
    angle3 = k * (360 - angle1 - angle2)
    angle4 = (1 - k) * (360 - angle1 - angle2)
    The function ensures that 0 < angle3 < 180 and 0 < angle4 < 180.

    Parameters:
    angle1 (float): The first angle (in degrees).
    angle2 (float): The second angle (in degrees).
    booster (float, optional): A scaling factor used to adjust the tilde angle. Default is 1.

    Returns:
    float: The generated tilde angle.
    """
    # Check for invalid input where the sum of angles exceeds 360 degrees
    if angle1 + angle2 >= 360:
        print("Angles Generating Error: The sum of two angles exceeds 360 degrees.")
        return None

    # If the sum of angle1 and angle2 is 180 degrees or more
    if angle1 + angle2 >= 180:
        # Calculate the factor k when the sum of angles is greater than or equal to 180
        k = min(booster * 0.5 * 180 / (360 - angle1 - angle2), 0.5)
        # Debugging output when booster is set to 1000
        if booster == 1000:
            print("Case 1: Sum of angles >= 180")
            print("k =", k)
    else:
        # Calculate k for when the sum of angles is less than 180
        k = ((1 - 0.6 * booster) * (180 - angle1 - angle2) + 0.6 * booster * 180) / (360 - angle1 - angle2)
        # Debugging output when booster is set to 1000
        if booster == 1000:
            print("Case 2: Sum of angles < 180")
            print("k =", k)

    return (360 - angle1 - angle2) * k

def convert_angles_to_vertices(Alphas, Gammas, Deltas, DihedralAngles, Edge):
    """
    Converts angles and edges into vertices for a quadrilateral in 3D space.

    The quadrilateral A1A2A3A4 lies in the xy-plane, with A2A1 along the x-axis, and the y-axis
    forms an acute angle with A2A3. The orientation of the faces is clockwise for faces 1, 2, 3, and 4, 
    while face C is counterclockwise. The vertices are calculated based on the given angles and edges 
    using trigonometric relationships (see reference PDF for details).

    Parameters:
    Alphas (list): A list of four alpha angles [alpha1, alpha2, alpha3, alpha4] (in degrees).
    Gammas (list): A list of four gamma angles [gamma1, gamma2, gamma3, gamma4] (in degrees).
    Deltas (list): A list of four delta angles [delta1, delta2, delta3, delta4] (in degrees).
    DihedralAngles (list): A list of four dihedral angles [phi, psi2, theta, psi1] (in degrees).
    Edge (float): The length of the edge A1A2.

    Returns:
    dict: A dictionary of vertex coordinates labeled as A1, A2, A3, A4, B1, B2, B3, B4, C1, C2, C3, C4.
    """
    # Initialize the 12 angles (we don't need betas to construct vertices)
    alpha1, alpha2, alpha3, alpha4 = map(np.float64, Alphas)
    gamma1, gamma2, gamma3, gamma4 = map(np.float64, Gammas)
    delta1, delta2, delta3, delta4 = map(np.float64, Deltas)

    # Initialize the four dihedral angles
    phi, psi2, theta, psi1 = map(np.float64, DihedralAngles)

    # Initialize the edge A1A2
    A1A2 = np.float64(Edge)

    # Initialize the dictionary to store the vertices
    vertices = {}

    # Plane C: Calculate the coordinates of vertices A1, A2, A3, A4
    a1 = A1A2
    a2 = generate_edge(a1, delta1, delta2, delta3)  # Edge A2A3
    xA1, yA1, zA1 = a1, 0.0, 0.0  # A1 lies on the x-axis
    xA2, yA2, zA2 = 0.0, 0.0, 0.0  # A2 is at the origin

    # Calculate vertex A3
    xA3 = xA2 + a2 * cos(delta2)
    yA3 = yA2 + a2 * sin(delta2)
    zA3 = zA2

    # Calculate vertex A4
    A1A4 = (a2 * sin(delta3) - a1 * sin(delta2 + delta3)) / sin(delta4)
    a3 = A1A4
    xA4 = xA1 - a3 * cos(delta1)
    yA4 = yA1 + a3 * sin(delta1)
    zA4 = zA1

    # Plane 1: Calculate the coordinates of vertices B1, B2
    alpha3Tilde = generate_tilde_angle(alpha1, alpha2)  # Angle A2B2B1
    alpha4Tilde = 360 - alpha1 - alpha2 - alpha3Tilde # Angle A1B1B2
    A2B2 = generate_edge(a1, alpha1, alpha2, alpha3Tilde, booster=0.4)
    b2 = A2B2
    A1B1 = (b2 * sin(alpha3Tilde) - a1 * sin(alpha2 + alpha3Tilde)) / sin(alpha4Tilde)
    b1 = A1B1

    # Calculate vertex B1
    xB1 = xA1 - b1 * cos(alpha1)
    yB1 = yA1 + b1 * sin(alpha1) * cos(phi)
    zB1 = zA1 + b1 * sin(alpha1) * sin(phi)

    # Calculate vertex B2
    xB2 = xA2 + b2 * cos(alpha2)
    yB2 = yA2 + b2 * sin(alpha2) * cos(phi)
    zB2 = zA2 + b2 * sin(alpha2) * sin(phi)

    # Plane 2: Calculate the coordinates of vertices C2, C3
    gamma1Tilde = generate_tilde_angle(gamma3, gamma2) # Angle A2C2C3
    gamma4Tilde = 360 - gamma2 - gamma3 - gamma1Tilde # Angle A3C3C2
    A2C2 = generate_edge(a2, gamma3, gamma2, gamma1Tilde, booster=50)
    c2 = A2C2
    A3C3 = (c2 * sin(gamma1Tilde) - a2 * sin(gamma2 + gamma1Tilde)) / sin(gamma4Tilde)
    c3 = A3C3

    # Calculate vertex C2
    xC2 = xA2 + c2 * (cos(gamma2) * cos(delta2) + sin(gamma2) * cos(psi2) * sin(delta2))
    yC2 = yA2 + c2 * (cos(gamma2) * sin(delta2) - sin(gamma2) * cos(psi2) * cos(delta2))
    zC2 = zA2 + c2 * sin(gamma2) * sin(psi2)

    # Calculate vertex C3
    xC3 = xA3 - c3 * (cos(gamma3) * cos(delta2) - sin(gamma3) * cos(psi2) * sin(delta2))
    yC3 = yA3 - c3 * (cos(gamma3) * sin(delta2) + sin(gamma3) * cos(psi2) * cos(delta2))
    zC3 = zA3 + c3 * sin(gamma3) * sin(psi2)

    # Plane 3: Calculate the coordinates of vertices B3, B4
    alpha1Tilde = generate_tilde_angle(alpha4, alpha3) # Angle A3B3B4
    alpha2Tilde = 360 - alpha3 - alpha4 - alpha1Tilde # Angle A4B4B3
    A3A4 = (a1 * sin(delta1) - a2 * sin(delta1 + delta2)) / sin(delta4)
    a4 = A3A4
    A3B3 = generate_edge(a4, alpha4, alpha3, alpha1Tilde, booster=0.4)
    b3 = A3B3
    A4B4 = (b3 * sin(alpha1Tilde) - a4 * sin(alpha3 + alpha1Tilde)) / sin(alpha2Tilde)
    b4 = A4B4

    # Calculate vertex B3
    xB3 = xA3 - b3 * (cos(alpha3) * cos(delta2 + delta3) + sin(alpha3) * cos(theta) * sin(delta2 + delta3))
    yB3 = yA3 - b3 * (cos(alpha3) * sin(delta2 + delta3) - sin(alpha3) * cos(theta) * cos(delta2 + delta3))
    zB3 = zA3 + b3 * sin(alpha3) * sin(theta)

    # Calculate vertex B4
    xB4 = xA4 + b4 * (cos(alpha4) * cos(delta2 + delta3) - sin(alpha4) * cos(theta) * sin(delta2 + delta3))
    yB4 = yA4 + b4 * (cos(alpha4) * sin(delta2 + delta3) + sin(alpha4) * cos(theta) * cos(delta2 + delta3))
    zB4 = zA4 + b4 * sin(alpha4) * sin(theta)

    # Plane 4: Calculate the coordinates of vertices C1, C4
    gamma2Tilde = generate_tilde_angle(gamma1, gamma4) # Angle A4C4C1
    gamma3Tilde = 360 - gamma1 - gamma4 - gamma2Tilde # Angle A1C1C4
    A4C4 = generate_edge(a3, gamma1, gamma4, gamma2Tilde, booster=25)
    c4 = A4C4
    A1C1 = (c4 * sin(gamma2Tilde) - a3 * sin(gamma4 + gamma2Tilde)) / sin(gamma3Tilde)
    c1 = A1C1

    # Calculate vertex C4
    xC4 = xA4 + c4 * (cos(gamma4) * cos(delta1) - sin(gamma4) * cos(psi1) * sin(delta1))
    yC4 = yA4 - c4 * (cos(gamma4) * sin(delta1) + sin(gamma4) * cos(psi1) * cos(delta1))
    zC4 = zA4 + c4 * sin(gamma4) * sin(psi1)

    # Calculate vertex C1
    xC1 = xA1 - c1 * (cos(gamma1) * cos(delta1) + sin(gamma1) * cos(psi1) * sin(delta1))
    yC1 = yA1 + c1 * (cos(gamma1) * sin(delta1) - sin(gamma1) * cos(psi1) * cos(delta1))
    zC1 = zA1 + c1 * sin(gamma1) * sin(psi1)

    # Store the vertices in the dictionary
    vertices[r"$A_1$"] = (xA1, yA1, zA1)
    vertices[r"$A_2$"] = (xA2, yA2, zA2)
    vertices[r"$A_3$"] = (xA3, yA3, zA3)
    vertices[r"$A_4$"] = (xA4, yA4, zA4)
    vertices[r"$B_1$"] = (xB1, yB1, zB1)
    vertices[r"$B_2$"] = (xB2, yB2, zB2)
    vertices[r"$B_3$"] = (xB3, yB3, zB3)
    vertices[r"$B_4$"] = (xB4, yB4, zB4)
    vertices[r"$C_1$"] = (xC1, yC1, zC1)
    vertices[r"$C_2$"] = (xC2, yC2, zC2)
    vertices[r"$C_3$"] = (xC3, yC3, zC3)
    vertices[r"$C_4$"] = (xC4, yC4, zC4)

    return vertices

def create_obj_file_of_vertices_and_faces(Vertices):
    """
    Creates an .obj file with the given vertices and faces.
    
    This function generates a 3D object file (OBJ format) by first converting 
    the vertex coordinates into the required format and then adding face definitions.
    
    Parameters:
    Vertices (dict): A dictionary of vertices where each key is a vertex name, 
                     and each value is a tuple of coordinates (x, y, z).
    
    Output:
    A file named "VerticesFaces.obj" containing the vertices and faces in OBJ format.
    """
    # Generate vertex text in OBJ format
    vText = ""
    for val in Vertices.values():
        vText += "v"  # Start each line with 'v' indicating a vertex in OBJ format
        for el in val:
            vText += f" {el}"  # Append each coordinate to the line
        vText += "\n"  # Move to the next line for the next vertex

    # Define the faces (polygon surfaces) for the object
    fText = (
        "f 1 2 3 4\n"  # Face connecting vertices 1, 2, 3, 4
        "f 1 5 6 2\n"  # Face connecting vertices 1, 5, 6, 2
        "f 2 6 10\n"  # Face connecting vertices 2, 6, 10
        "f 2 10 11 3\n"  # Face connecting vertices 2, 10, 11, 3
        "f 3 11 7\n"  # Face connecting vertices 3, 11, 7
        "f 3 7 8 4\n"  # Face connecting vertices 3, 7, 8, 4
        "f 4 8 12\n"  # Face connecting vertices 4, 8, 12
        "f 4 12 9 1\n"  # Face connecting vertices 4, 12, 9, 1
        "f 1 9 5"  # Face connecting vertices 1, 9, 5
    )

    # Combine the vertex and face information into a single OBJ string
    objText = vText + fText

    # Write the text to an .obj file
    with open("VerticesFaces.obj", "w") as objFile:
        objFile.write(objText)

def for_geogebra(Vertices):
    """
    Creates a text file formatted for use with GeoGebra, including the vertices 
    and polygons that form the faces of a 3D shape.
    
    Parameters:
    Vertices (dict): A dictionary of vertices where each key is a vertex name, 
                     and each value is a tuple of coordinates (x, y, z).
    
    Output:
    A file named "forGeoGebra.txt" containing the vertices and polygons in GeoGebra format.
    """
    # Generate the vertex declarations for GeoGebra
    verticesText = ""
    for val in Vertices.keys():
        verticesText += f"{val}=" + f"{Vertices[val]}" + "\n"

    # Define the polygons (faces) for GeoGebra using the defined vertices
    polygonText = (
        "\nPolygon($A_1$,$A_2$,$A_3$,$A_4$)\n"  # Polygon face A1A2A3A4
        "Polygon($A_1$,$B_1$,$B_2$,$A_2$)\n"  # Polygon face A1B1B2A2
        "Polygon($A_2$,$B_2$,$C_2$)\n"  # Polygon face A2B2C2
        "Polygon($A_2$,$C_2$,$C_3$,$A_3$)\n"  # Polygon face A2C2C3A3
        "Polygon($A_3$,$C_3$,$B_3$)\n"  # Polygon face A3C3B3
        "Polygon($A_3$,$B_3$,$B_4$,$A_4$)\n"  # Polygon face A3B3B4A4
        "Polygon($A_4$,$B_4$,$C_4$)\n"  # Polygon face A4B4C4
        "Polygon($A_4$,$C_4$,$C_1$,$A_1$)\n"  # Polygon face A4C4C1A1
        "Polygon($A_1$,$C_1$,$B_1$)"  # Polygon face A1C1B1
    )

    # Combine the vertex declarations and polygon definitions
    geoGebraText = verticesText + polygonText

    # Write the text to a file
    with open("forGeoGebra.txt", "w") as geoGebraFile:
        geoGebraFile.write(geoGebraText)

def visualize_here(Vertices):
    """
    Visualizes 3D vertices and faces using Matplotlib.

    This function creates a 3D plot of the given vertices and connects them to form faces. 
    Each vertex is labeled, and the 3D model can be rotated with mouse events. 
    Additionally, key events can toggle mouse rotation.

    Parameters:
    Vertices (dict): A dictionary where the keys are vertex labels and the values are tuples of vertex coordinates (x, y, z).
    
    Interaction:
    - Press 'z' to toggle mouse rotation.
    """
    global verticesList, ax, fig, L, labelsList, x2DList, y2DList, rotation_enabled

    # Prepare vertices and labels
    verticesList = [list(val) for val in Vertices.values()]
    verticesLabels = list(Vertices.keys())

    # Define faces (polygon surfaces) using vertex indices
    Faces = [[0, 1, 2, 3], [0, 4, 5, 1], [1, 9, 10, 2], [2, 6, 7, 3], [3, 11, 8, 0]]
    AllFaces = [[0, 1, 2, 3], [0, 4, 5, 1], [1, 5, 9], [1, 9, 10, 2], [2, 10, 6], 
                [2, 6, 7, 3], [3, 7, 11], [3, 11, 8, 0], [0, 8, 4]]

    # Map vertex indices to their coordinates
    VerticesToDraw = [[verticesList[i] for i in face] for face in Faces]
    AllVerticesToDraw = [[verticesList[i] for i in face] for face in AllFaces]

    # Define face colors
    Colors = ['slategrey', 'crimson', 'slateblue', 'turquoise', 'orangered']
    AllColors = ['slategrey', 'crimson', 'slateblue', 'turquoise', 'orangered', 
                 'seagreen', 'purple', 'royalblue', 'teal']

    # Initialize the 3D figure
    fig = plt.figure(figsize=(10, 8))
    ax = plt.axes(projection='3d')
    ax.set_axis_off()  # Turn off the axes
    ax.set(xlim=(-6, 10), ylim=(-6, 10), zlim=(-6, 10))  # Set axis limits

    # Create the 3D polygon collection for the faces
    # collection = Poly3DCollection(VerticesToDraw, color = Colors, linewidths = 1.1, alpha = 0.05) 
    # collection = Poly3DCollection(AllVerticesToDraw, color = AllColors, linewidths = 1, alpha = 0.05)
    collection = Poly3DCollection(AllVerticesToDraw, edgecolor='darkblue', linewidths=1.5, alpha=0.5)
    ax.add_collection3d(collection)  # Add the collection to the plot

    # Store the number of vertices and initialize labels and positions for annotation
    L = len(verticesList)
    labelsList = [None] * L
    x2DList, y2DList = [0] * L, [0] * L

    # Adding labels for each vertex
    for i in range(L):
        x2DList[i], y2DList[i], _ = proj3d.proj_transform(verticesList[i][0], verticesList[i][1], verticesList[i][2], ax.get_proj())
        labelsList[i] = pylab.annotate(f"{verticesLabels[i]}", xy=(x2DList[i], y2DList[i]), xytext=(0, 0),
                                       color='darkblue', fontsize=16, textcoords='offset points', ha='right', va='bottom')
        labelsList[i].draggable()  # Make labels draggable
        ax.scatter(verticesList[i][0], verticesList[i][1], verticesList[i][2], color='darkblue', s=30)  # Plot vertices

    # Initial state of mouse rotation (enabled)
    rotation_enabled = True
    fig.canvas.mpl_connect('key_press_event', onkey)  # Connect key press event
    fig.canvas.mpl_connect('button_release_event', update_position)  # Connect mouse release event to update labels
    plt.show()

def update_position(event):
    """
    Updates the position of vertex labels in 2D space based on the 3D view.

    This function is called on mouse release to adjust the position of vertex labels 
    after the 3D plot has been rotated or moved.

    Parameters:
    event (Event): The Matplotlib event triggered by a mouse release.
    """
    for i in range(L):
        x2DList[i], y2DList[i], _ = proj3d.proj_transform(verticesList[i][0], verticesList[i][1], verticesList[i][2], ax.get_proj())
        labelsList[i].xy = x2DList[i], y2DList[i]
        labelsList[i].update_positions(fig.canvas.renderer)
    fig.canvas.draw()

def onkey(event):
    """
    Toggles mouse rotation for the 3D plot on and off.

    This function is triggered when the 'z' key is pressed and disables or enables 
    mouse rotation of the 3D plot.

    Parameters:
    event (Event): The Matplotlib event triggered by a key press.
    """
    global rotation_enabled
    if event.key == 'z':  # Press 'z' to toggle rotation
        if rotation_enabled:
            ax.disable_mouse_rotation()  # Disable rotation
        else:
            ax.mouse_init()  # Enable rotation
        rotation_enabled = not rotation_enabled  # Toggle the flag
        plt.draw()
        
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
        
def main():
    """
    Main function to set up the parameters, compute the vertices from angles, 
    and generate outputs for visualization and file creation.
    
    The function initializes the angles, edge lengths, and dihedral angles, 
    then calls necessary functions to:
    1. Convert the angles into 3D vertices.
    2. Export the vertices for GeoGebra.
    3. Create an .obj file for 3D visualization.
    4. Visualize the 3D shape in a Matplotlib window.
    """

    # SET-UP: Define angles and edge lengths for the shape
    # Read input values from file
    file_path = 'optimal_angles_2024_09_23_04_32_08.txt'  # Replace with the actual file path
    values = read_values_from_file(file_path)
    
    # Unpack parameters
    Alphas = values['Alphas']
    Betas = values['Betas']
    Gammas = values['Gammas']
    Deltas = values['Deltas']
    DihedralAngles = values['DihedralAngles']

    # A1A2: length of the edge A1A2
    Edge = 4

    # RESULT: Convert angles into 3D vertices and visualize the result
    # Convert angles and edge length to vertices
    Vertices = convert_angles_to_vertices(Alphas, Gammas, Deltas, DihedralAngles, Edge)

    # Generate a file for GeoGebra containing the vertices and polygons
    # for_geogebra(Vertices)

    # Create an .obj file for use in 3D modeling software
    # create_obj_file_of_vertices_and_faces(Vertices)

    # Visualize the 3D vertices and faces in a Matplotlib plot
    visualize_here(Vertices)


if __name__ == "__main__":
    main()