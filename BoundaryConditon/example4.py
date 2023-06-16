import sys
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
import numpy as np


def set_boundary_conditions(surface_vessel):
    if surface_vessel:
        # Set boundary conditions for surface vessel simulation
        inflow_conditions = input("Inflow conditions: ")
        outflow_conditions = input("Outflow conditions: ")
        free_surface_conditions = input("Free surface conditions: ")
        hull_boundary_conditions = input("Hull boundary conditions: ")

        boundary_conditions = {
            'inflow': inflow_conditions,
            'outflow': outflow_conditions,
            'free_surface': free_surface_conditions,
            'hull_boundary': hull_boundary_conditions
        }
    else:
        # Set boundary conditions for underwater vessel simulation
        inflow_conditions = input("Inflow conditions: ")
        submerged_boundary_conditions = input("Submerged boundary conditions: ")
        outflow_conditions = input("Outflow conditions: ")

        boundary_conditions = {
            'inflow': inflow_conditions,
            'submerged_boundary': submerged_boundary_conditions,
            'outflow': outflow_conditions
        }

    return boundary_conditions


def create_simulation_environment(mesh_file, surface_vessel):
    # Load the mesh file
    mesh = load_mesh(mesh_file)

    # Create an instance of the simulation environment
    simulation = Simulation(mesh)

    # Set the boundary conditions based on the type of vessel
    boundary_conditions = set_boundary_conditions(surface_vessel)

    # Apply the boundary conditions to the simulation environment
    simulation.apply_boundary_conditions(boundary_conditions)

    # Define the blue strip of water and waves
    strip_height = 0.2  # Height of the blue strip
    wave_amplitude = 0.05  # Amplitude of the waves
    wave_length = 0.5  # Length of the waves

    # Generate the mesh points for the blue strip
    strip_points = generate_strip_points(mesh)

    # Generate the wave displacement field
    wave_displacement = generate_wave_displacement(mesh, strip_points, strip_height, wave_amplitude, wave_length)

    # Apply the wave displacement field to the simulation environment
    simulation.apply_wave_displacement(wave_displacement)

    # Run the simulation
    simulation.run()

    # Additional steps for post-processing and analysis


def load_mesh(mesh_file):
    # Load the mesh file
    mesh = ParsedParameterFile(mesh_file)
    return mesh


def generate_strip_points(mesh):
    # Generate the mesh points for the blue strip
    x_coordinates = mesh['x']  # Extract the x-coordinates from the mesh
    y_coordinates = mesh['y']  # Extract the y-coordinates from the mesh
    strip_points = np.column_stack((x_coordinates, y_coordinates))

    return strip_points


def generate_wave_displacement(mesh, strip_points, strip_height, wave_amplitude, wave_length):
    # Generate the wave displacement field
    num_points = len(strip_points)
    wave_displacement = np.zeros((num_points, 3))

    for i, point in enumerate(strip_points):
        x = point[0]
        y = point[1]
        z = strip_height * np.sin(2 * np.pi * x / wave_length) * wave_amplitude
        wave_displacement[i] = [x, y, z]

    return wave_displacement


class Simulation:
    def __init__(self, mesh):
        self.mesh = mesh

    def apply_boundary_conditions(self, boundary_conditions):
        # Apply the boundary conditions to the simulation environment
        for boundary, condition in boundary_conditions.items():
            # Apply each boundary condition based on its type
            self.mesh.boundaryField()[boundary] = condition

    def apply_wave_displacement(self, wave_displacement):
        # Apply the wave displacement field to the simulation environment
        points = self.mesh.points()
        for i, point in enumerate(points):
            point[2] += wave_displacement[i][2]

    def run(self):
        # Run the simulation
        self.mesh.write()


if __name__ == "__main__":
    try:
        # Check if the mesh file path is provided as a command line argument
        if len(sys.argv) < 2:
            print("Please provide the path to the mesh file as a command line argument.")
            sys.exit(1)

        # Read the mesh file path from the command line argument
        mesh_file = sys.argv[1]

        # Prompt the user to input the vessel type
        surface_vessel = input("Is the vehicle a surface vessel? (Y/N): ").lower() == 'y'

        # Prompt the user to set the boundary conditions
        boundary_conditions = set_boundary_conditions(surface_vessel)

        # Create the simulation environment with the user-defined boundary conditions
        create_simulation_environment(mesh_file, surface_vessel)

        # Output a message indicating the successful execution
        print("Simulation completed successfully.")

    except Exception as e:
        # Output an error message if an exception occurs
        print("An error occurred during simulation:", str(e))
