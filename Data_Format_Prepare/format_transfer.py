import pyvista as pv

import os
import numpy as np


import numpy as np

def mesh2off(vertices, faces, filename):
    """
    Save a triangular mesh to an OFF file.

    Parameters
    ----------
    vertices : (N, 3) array-like
        Vertex coordinates (float).
    faces : (M, 3) array-like
        Triangle vertex indices (integers, 0-based).
    filename : str
        Output OFF file path.
    """
    vertices = np.asarray(vertices, dtype=float)
    faces = np.asarray(faces, dtype=int)

    num_vertices = vertices.shape[0]
    num_faces = faces.shape[0]

    with open(filename, "w") as file:
        # Header
        file.write("OFF\n")
        file.write(f"{num_vertices} {num_faces} 0\n")

        # Write vertices
        for v in vertices:
            file.write(f"{v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")

        # Write faces (ensure integer indices)
        for f in faces:
            file.write(f"{int(f[0])} {int(f[1])} {int(f[2])} {int(f[3])}\n")

    print(f"✅ OFF file successfully saved to: {filename}")




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
            file.write(f"{int(cell[0])} {int(cell[1])} {int(cell[2])} {int(cell[3])}  {int(cell[4])}\n")

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

    # volume volume_surf tgt_pts









