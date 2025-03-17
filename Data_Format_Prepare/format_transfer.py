import pyvista as pv

import os
import numpy as np

def mesh2off(vertices, faces, filename):
    # Check if the number of vertices matches the length of the x, y, z coordinates.
    num_vertices = len(vertices)
    num_faces = len(faces)


    # Open the .off file for writing.
    with open(filename, "w") as file:
        # Write the OFF file header.
        file.write("OFF\n")
        file.write(f"{num_vertices } {num_faces} 0\n")  # num_vertices // 3 faces, 0 edges, 0 other entries

        # Write vertex coordinates.
        for vertex in vertices:
            #print(vertex)
            file.write(f"{vertex[0]} {vertex[1]} {vertex[2]}\n")

        for face in faces:
            file.write(f"{face[0]} {face[1]} {face[2]} {face[3]}\n")

        file.close()

def pts2off(vertices, filename):
    # Check if the number of vertices matches the length of the x, y, z coordinates.
    num_vertices = len(vertices)

    # Open the .off file for writing.
    with open(filename, "w") as file:
        # Write the OFF file header.
        file.write("OFF\n")
        file.write(f"{num_vertices } 0 0\n")  # num_vertices // 3 faces, 0 edges, 0 other entries

        # Write vertex coordinates.
        for vertex in vertices:
            #print(vertex)
            file.write(f"{vertex[0]} {vertex[1]} {vertex[2]}\n")

        file.close()

def vol2off(vertices, cells, filename):
    # Check if the number of vertices matches the length of the x, y, z coordinates.
    num_vertices = len(vertices)
    num_cells = len(cells)

    # Open the .off file for writing.
    with open(filename, "w") as file:
        # Write the OFF file header.
        file.write("OFF\n")
        file.write(f"{num_vertices } {num_cells} 0\n")  # num_vertices // 3 faces, 0 edges, 0 other entries

        # Write vertex coordinates.
        for vertex in vertices:
            #print(vertex)
            file.write(f"{vertex[0]} {vertex[1]} {vertex[2]}\n")

        for cell in cells:
            #print(vertex)
            file.write(f"{cell[0]} {cell[1]} {cell[2]} {cell[3]}  {cell[4]}\n")

        file.close()


if __name__ == '__main__':

    init_file = "/media/yzx/yzx_store1/Dataset/Kelly_all_dataset/P5/Flat/volume.vtk"

    # vtk2off
    vol = pv.read(init_file)
    points = vol.points
    cells_dict = vol.cells_dict
    tet_cell = cells_dict[10]  # v1 v2 v3 v4
    pad_cell = np.ones([len(tet_cell), 5]) * 4  # 5 v1 v2 v3 v4
    pad_cell[:, 1:] = tet_cell
    vol2off(points, pad_cell, init_file.replace(".vtk", ".off"))

    # mesh2off
    mesh = vol.extract_surface()
    pts = mesh.points
    faces = mesh.faces
    faces = faces.reshape(len(faces)//4, 4) # 3 v1 v2 v3
    mesh2off(pts, faces, init_file.replace(".vtk","_surf.off"))

    #pts2off
    pts2off(pts, init_file.replace(".vtk","_pts.off"))










