import numpy as np


def generate_alpha_water(file_path, resx=100, resy=100, diameter_fraction=0.5, center_x=0.5, center_y=0.5):
    # Calculate the radius of the circle
    radius_x = (diameter_fraction * resx) / 2
    radius_y = radius_x
    
    # Create the grid
    grid = np.ones((resx, resy))
    
    # Populate the grid with zeros within the circle
    for i in range(resx):
        for j in range(resy):
            if ((i - center_x) / radius_x)**2 + ((j - center_y) / radius_y)**2 <= 1:
                grid[i, j] = 0.0
    
    # Write the grid to the file in OpenFOAM format
    with open(file_path, 'w') as file:
        file.write("/*--------------------------------*- C++ -*----------------------------------*\\\n")
        file.write("| =========                 |                                                 |\n")
        file.write("| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n")
        file.write("|  \\    /   O peration     | Version:  7                                     |\n")
        file.write("|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |\n")
        file.write("|    \\/     M anipulation  |                                                 |\n")
        file.write("\\*---------------------------------------------------------------------------*/\n")
        file.write("/* File: alpha.water\n")
        file.write(" * Description: Generated scalar field file\n")
        file.write(" * ---------------------------------------------------------------------------\n")
        file.write(" */\n")
        file.write("\n")
        file.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n")
        file.write("\n")
        file.write("FoamFile\n")
        file.write("{\n")
        file.write("    version     2.0;\n")
        file.write("    format      ascii;\n")
        file.write("    class       volScalarField;\n")
        file.write("    location    \"0\";\n")
        file.write("    object      alpha.water;\n")
        file.write("}\n")
        file.write("\n")
        file.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n")
        file.write("\n")
        file.write("dimensions      [0 0 0 0 0 0 0];\n")
        file.write(f"internalField   nonuniform List<scalar>\n{resx * resy}\n(\n")
        
        # Write the grid values
        for value in grid.flatten():
            file.write(f"{value}\n")
        
        file.write(");\n")
        file.write("boundaryField\n")
        file.write("{\n")
        file.write("    front\n")
        file.write("    {\n")
        file.write("        type            zeroGradient;\n")
        file.write("    }\n")
        file.write("    back\n")
        file.write("    {\n")
        file.write("        type            zeroGradient;\n")
        file.write("    }\n")
        file.write("    left\n")
        file.write("    {\n")
        file.write("        type            zeroGradient;\n")
        file.write("    }\n")
        file.write("    right\n")
        file.write("    {\n")
        file.write("        type            zeroGradient;\n")
        file.write("    }\n")
        file.write("    topAndBottom\n")
        file.write("    {\n")
        file.write("        type            empty;\n")
        file.write("    }\n")
        file.write("}\n")
        file.write("\n")
        file.write("// ************************************************************************* //\n")


# Example usage
generate_alpha_water('0.orig/alpha.water', resx=200, resy=120, diameter_fraction=0.1, center_x=20, center_y=20)

