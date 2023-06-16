#!/usr/bin/env python3

import os
import sys
import argparse
import meshio
import pyvista as pv

def convert_file(input_file, output_file):
    # Load the mesh from the input file
    mesh = meshio.read(input_file)

    # Write the mesh to the output file in the desired format
    meshio.write(output_file, mesh, 'gmsh')

def display_mesh(output_file):
    # Read the mesh file
    mesh = pv.read(output_file)

    plotter = pv.Plotter()
    plotter.add_mesh(mesh, color="white", show_edges=True)
    plotter.set_background("white")
    plotter.show_axes()
    plotter.show()

if __name__ == '__main__':
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Convert an STL file to a mesh file for CFD simulation.')
    parser.add_argument('input_file', help=r'C:\Users\Roshan George\Documents\Meshing\test1\ship.stl')
    parser.add_argument('output_file', help=r'C:\Users\Roshan George\Documents\Meshing\test1\shipOut.msh')
    args = parser.parse_args()

    # Check if the input file exists
    if not os.path.isfile(args.input_file):
        print(f'Error: Input file "{args.input_file}" not found.')
        sys.exit(1)

    # Convert the file
    try:
        convert_file(args.input_file, args.output_file)
        print(f'Successfully converted "{args.input_file}" to "{args.output_file}" in Gmsh format for CFD simulation.')
        display_mesh(args.output_file)
    except Exception as e:
        print(f'Error: {str(e)}')
        sys.exit(1)

