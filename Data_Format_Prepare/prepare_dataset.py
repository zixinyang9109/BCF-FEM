
from stl import mesh
import trimesh
import tetgen
import open3d as o3d
import pyvista as pv
import numpy as np
import os
from format_transfer import vol2off, mesh2off, pts2off

def extract_surface_from_tets(tet_vs, tet_elements):
    """
    Extract surface mesh (triangles) from a tetrahedral mesh.

    Parameters:
    -----------
    tet_vs : (N, 3) array
        Vertex coordinates of the tetrahedral mesh.
    tet_elements : (M, 4) array
        Tetrahedral connectivity (each row contains 4 vertex indices).

    Returns:
    --------
    surface : pv.PolyData
        The extracted surface mesh.
    """
    # Convert to numpy arrays
    tet_vs = np.asarray(tet_vs)
    tet_elements = np.asarray(tet_elements, dtype=np.int32)

    # Create an UnstructuredGrid (PyVista’s volumetric representation)
    n_tets = tet_elements.shape[0]
    cells = np.hstack([
        np.full((n_tets, 1), 4),  # number of points per cell (always 4 for a tet)
        tet_elements
    ]).flatten()

    celltypes = np.full(n_tets, pv.CellType.TETRA, dtype=np.uint8)

    grid = pv.UnstructuredGrid(cells, celltypes, tet_vs)

    # Extract the outer surface as a PolyData (triangular mesh)
    surface = grid.extract_surface()

    return surface

def filter_point_cloud_by_distance(target_points, reference_points, distance_threshold):
    """
    Remove points from target point cloud that are farther than `distance_threshold`
    from the reference point cloud.

    Parameters
    ----------
    target_points : (N,3) np.ndarray
        Point cloud to be cleaned.
    reference_points : (M,3) np.ndarray
        Reference point cloud.
    distance_threshold : float
        Maximum allowed distance to the reference point cloud.

    Returns
    -------
    filtered_points : (K,3) np.ndarray
        Cleaned point cloud.
    """
    # Convert to Open3D point clouds
    pcd_target = o3d.geometry.PointCloud()
    pcd_target.points = o3d.utility.Vector3dVector(target_points)

    pcd_ref = o3d.geometry.PointCloud()
    pcd_ref.points = o3d.utility.Vector3dVector(reference_points)

    # Build KDTree on reference points
    pcd_ref_tree = o3d.geometry.KDTreeFlann(pcd_ref)

    # Keep points that have a nearest neighbor within the threshold
    keep_mask = []
    for pt in target_points:
        [_, idx, dist2] = pcd_ref_tree.search_knn_vector_3d(pt, 1)
        if np.sqrt(dist2[0]) <= distance_threshold:
            keep_mask.append(True)
        else:
            keep_mask.append(False)

    filtered_points = target_points[keep_mask, :]
    return filtered_points

def compare_filtered_unfiltered(unfiltered_points, filtered_points):
    """
    Visualize unfiltered and filtered point clouds with labels/captions.

    Parameters
    ----------
    unfiltered_points : (N,3) np.ndarray
        Original point cloud (red).
    filtered_points : (M,3) np.ndarray
        Filtered point cloud (blue).
    """
    pl = pv.Plotter()

    # Unfiltered: red
    pl.add_points(unfiltered_points, color='red', point_size=5, render_points_as_spheres=True, opacity=0.6)
    # Add a text label
    pl.add_text("Unfiltered (Red)", position='upper_left', font_size=14, color='red')

    # Filtered: blue
    pl.add_points(filtered_points, color='blue', point_size=5, render_points_as_spheres=True, opacity=0.8)
    pl.add_text("Filtered (Blue)", position='upper_right', font_size=14, color='blue')

    pl.add_axes()
    pl.show_grid()
    pl.show()

def simplify_mesh(model_vs, model_fs, target_reduction=0.5):
    """
    Simplify a mesh using Open3D's quadric decimation.

    Parameters
    ----------
    model_vs : (N, 3) np.ndarray
        Input vertices.
    model_fs : (M, 3) np.ndarray
        Input triangle faces.
    target_reduction : float
        Fraction of triangles to remove (0.5 = reduce by 50%).

    Returns
    -------
    new_vs : np.ndarray
        Simplified vertices.
    new_fs : np.ndarray
        Simplified faces.
    mesh_out : o3d.geometry.TriangleMesh
        Open3D simplified mesh object.
    """
    # Create Open3D mesh
    mesh = o3d.geometry.TriangleMesh(
        vertices=o3d.utility.Vector3dVector(model_vs),
        triangles=o3d.utility.Vector3iVector(model_fs)
    )
    mesh.compute_vertex_normals()

    # Target triangle count
    target_triangles = max(4, int(len(model_fs) * (1 - target_reduction)))

    # Simplify
    mesh_simpl = mesh.simplify_quadric_decimation(target_number_of_triangles=target_triangles)

    # Clean up mesh
    mesh_simpl.remove_degenerate_triangles()
    mesh_simpl.remove_duplicated_triangles()
    mesh_simpl.remove_unreferenced_vertices()
    mesh_simpl.remove_non_manifold_edges()
    mesh_simpl.compute_vertex_normals()

    # Convert back to numpy arrays
    new_vs = np.asarray(mesh_simpl.vertices)
    new_fs = np.asarray(mesh_simpl.triangles)

    return new_vs, new_fs, mesh_simpl

def read_ply(file_path):
    """
    Reads a PLY file (mesh or point cloud) and returns vertices and faces if available.

    Parameters
    ----------
    file_path : str
        Path to the PLY file.

    Returns
    -------
    vertices : (N, 3) np.ndarray
        Array of vertex coordinates.
    faces : (M, 3) np.ndarray or None
        Array of triangular faces (if mesh). None if the PLY only has point data.
    """
    mesh = o3d.io.read_triangle_mesh(file_path)
    if len(mesh.vertices) == 0:  # Might be a point cloud instead
        pcd = o3d.io.read_point_cloud(file_path)
        return np.asarray(pcd.points), None
    else:
        vertices = np.asarray(mesh.vertices)
        faces = np.asarray(mesh.triangles)
        return vertices, faces

def read_stl(file_path):
    """
    Reads an STL file and returns vertices and faces as numpy arrays.

    Parameters
    ----------
    file_path : str
        Path to the STL file.

    Returns
    -------
    vertices : (N, 3) np.ndarray
        Array of vertex coordinates.
    faces : (M, 3) np.ndarray
        Array of vertex indices defining each triangular face.
    """
    stl_mesh = mesh.Mesh.from_file(file_path)
    # STL files store triangles as three vertices per face
    vertices = np.unique(stl_mesh.vectors.reshape(-1, 3), axis=0)

    # Map vertices to indices to form faces
    vertex_to_index = {tuple(v): i for i, v in enumerate(vertices)}
    faces = np.array([[vertex_to_index[tuple(v)] for v in tri] for tri in stl_mesh.vectors])

    return vertices, faces

def vox_pts(points, voxel_size_normalized=0.02):
    """
    Normalize the point cloud, voxelize it in normalized space,
    and return the unnormalized (original-scale) downsampled point cloud as NumPy arrays.

    Parameters
    ----------
    points : (N, 3) np.ndarray
        Input point cloud coordinates.
    voxel_size_normalized : float, optional
        Voxel size in normalized coordinates (default 0.02).

    Returns
    -------
    downsampled_points : (M, 3) np.ndarray
        Downsampled (unnormalized) point cloud.
    normalized_points : (N, 3) np.ndarray
        The normalized version of the original point cloud.
    norm_params : dict
        Dictionary with normalization info: {'center', 'scale'}.
    """
    if not isinstance(points, np.ndarray):
        points = np.asarray(points)

    # Compute normalization parameters
    center = points.mean(axis=0)
    scale = np.max(np.linalg.norm(points - center, axis=1))

    # Normalize to unit sphere
    points_norm = (points - center) / scale

    # Convert to Open3D format for voxelization
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(points_norm)

    # Voxelize in normalized space
    pcd_down = pcd.voxel_down_sample(voxel_size=voxel_size_normalized)
    down_points_norm = np.asarray(pcd_down.points)

    # Unnormalize back to original space
    down_points = down_points_norm * scale + center

    return down_points

def generate_tet_from_mesh(vertices, faces, max_volume=None):
    """
    Generate a tetrahedral mesh from a watertight surface mesh using TetGen via PyVista.
    """
    # convert faces for PyVista: prepend '3' per triangle
    faces_flat = np.hstack([np.full((faces.shape[0], 1), 3), faces]).flatten()

    surf_mesh = pv.PolyData(vertices, faces_flat)
    tgen_obj = tetgen.TetGen(surf_mesh)

    # optional switches
    if max_volume is not None:
        nodes, elems = tgen_obj.tetrahedralize(switches=f"pq1.4/10Aaa{max_volume}")
    else:
        nodes, elems = tgen_obj.tetrahedralize()

    # nodes is an (N,3) array; elems is an (M,4) array of indices
    return nodes, elems

def simplify_and_repair_mesh(model_vs, model_fs, target_reduction=0.5):
    """
    Simplify a mesh using Open3D and repair it if not watertight using Trimesh.

    Parameters
    ----------
    model_vs : (N, 3) np.ndarray
        Input vertices.
    model_fs : (M, 3) np.ndarray
        Input triangle faces.
    target_reduction : float
        Fraction of triangles to remove (0.5 = reduce 50%).

    Returns
    -------
    new_vs : np.ndarray
        Simplified (and possibly repaired) vertices.
    new_fs : np.ndarray
        Simplified (and possibly repaired) faces.
    is_watertight : bool
        True if the final mesh is watertight.
    """
    # --- Simplify with Open3D ---
    mesh_o3d = o3d.geometry.TriangleMesh(
        vertices=o3d.utility.Vector3dVector(model_vs),
        triangles=o3d.utility.Vector3iVector(model_fs)
    )
    mesh_o3d.compute_vertex_normals()

    target_triangles = max(4, int(len(model_fs) * (1 - target_reduction)))
    mesh_simpl = mesh_o3d.simplify_quadric_decimation(target_number_of_triangles=target_triangles)
    mesh_simpl.remove_degenerate_triangles()
    mesh_simpl.remove_duplicated_triangles()
    mesh_simpl.remove_unreferenced_vertices()
    mesh_simpl.remove_non_manifold_edges()

    vs_s = np.asarray(mesh_simpl.vertices)
    fs_s = np.asarray(mesh_simpl.triangles)

    # --- Convert to Trimesh for repair ---
    mesh_tri = trimesh.Trimesh(vertices=vs_s, faces=fs_s, process=False)

    # Check watertightness
    is_watertight = mesh_tri.is_watertight

    # Repair if not watertight
    if not is_watertight:
        mesh_tri.fill_holes()
        mesh_tri.remove_degenerate_faces()
        mesh_tri.remove_duplicate_faces()
        mesh_tri.remove_infinite_values()
        is_watertight = mesh_tri.is_watertight

    # Return final vertices/faces
    return mesh_tri.vertices, mesh_tri.faces, is_watertight

def visualize_all_pyvista(points, mesh_vertices, mesh_faces, tet_vertices, tet_elements):
    """
    Visualize downsampled point cloud, surface mesh, and tetrahedral mesh in one PyVista window.

    Parameters
    ----------
    points : (N,3) np.ndarray
        Downsampled point cloud.
    mesh_vertices : (M,3) np.ndarray
        Surface mesh vertices.
    mesh_faces : (K,3) np.ndarray
        Surface mesh triangles.
    tet_vertices : (P,3) np.ndarray
        Tetrahedral mesh vertices.
    tet_elements : (Q,4) np.ndarray
        Tetrahedral mesh connectivity.
    """
    pl = pv.Plotter()

    # --- Point cloud ---
    pl.add_points(points, color='red', point_size=5, render_points_as_spheres=True)

    # --- Surface mesh ---
    # PyVista faces array: prepend 3 to each triangle
    faces_flat = np.hstack([np.full((mesh_faces.shape[0], 1), 3), mesh_faces]).flatten()
    surface_mesh = pv.PolyData(mesh_vertices, faces_flat)
    pl.add_mesh(surface_mesh, color='green', opacity=0.5, show_edges=True)

    # --- Tetrahedral mesh ---
    visualize_tet_edges(tet_vertices, tet_elements)

    # --- Final visualization ---
    pl.add_axes()
    pl.show_grid()
    pl.show()

def visualize_tet_edges(tet_vertices, tet_elements):
    """
    Visualize all tetrahedral edges (including interior) as lines.
    """
    # Create a PyVista PolyData to store lines
    lines = []

    # Tet indices: each row has 4 vertices
    for tet in tet_elements:
        # Each tet has 6 edges: (0-1,0-2,0-3,1-2,1-3,2-3)
        edges = [
            (tet[0], tet[1]),
            (tet[0], tet[2]),
            (tet[0], tet[3]),
            (tet[1], tet[2]),
            (tet[1], tet[3]),
            (tet[2], tet[3])
        ]
        lines.extend(edges)

    # Flatten for PyVista: [N, 3] => [N*3] with number of points per line first
    n_lines = len(lines)
    lines_flat = np.hstack([[2, *edge] for edge in lines])

    # Create PolyData
    tet_lines = pv.PolyData()
    tet_lines.points = tet_vertices
    tet_lines.lines = lines_flat

    # Plot
    pl = pv.Plotter()
    pl.add_mesh(tet_lines, color='gray', line_width=1.0)
    pl.add_axes()
    pl.show_grid()
    pl.show()

if __name__ == '__main__':

    root_path = r"C:\Users\Yang\Desktop\Win_projects\BCF-FEM-main\Patients"
    vis = False

    index = 5
    model_path = os.path.join(root_path, str(index), "3Dmodel.stl")
    pts_path = os.path.join(root_path, str(index), "pts.ply")

    model_vs, model_fs = read_stl(model_path)
    pts, _ = read_ply(pts_path)
    pts_down = vox_pts(pts, voxel_size_normalized=0.02)

    # Simplify mesh to 30% of original triangles
    vs_s, fs_s, mesh_simpl = simplify_and_repair_mesh(model_vs, model_fs, target_reduction=0.9)

    pts_down_filtered= filter_point_cloud_by_distance(pts_down, vs_s, 10)

    print(f"Original faces: {len(model_fs)}, Simplified faces: {len(fs_s)}")
    print(f"Watertight: {mesh_simpl}")

    # Suppose vs_s, fs_s are your simplified (and repaired) mesh
    tet_vs, tet_elements = generate_tet_from_mesh(vs_s, fs_s, max_volume=200.0) # max volume the larger the better

    print("Number of nodes:", tet_vs.shape[0])
    print("Number of tetrahedra:", tet_elements.shape[0])

    if vis:
        compare_filtered_unfiltered(pts_down, pts_down_filtered)
        visualize_all_pyvista(
            points=pts_down,
            mesh_vertices=vs_s,
            mesh_faces=fs_s,
            tet_vertices=tet_vs,
            tet_elements=tet_elements
        )

        visualize_all_pyvista(
            points=pts_down_filtered,
            mesh_vertices=vs_s,
            mesh_faces=fs_s,
            tet_vertices=tet_vs,
            tet_elements=tet_elements
        )

    # save results
    save_vol_file = os.path.join(root_path, str(index), "volume.off")
    save_surf_file = os.path.join(root_path, str(index), "volume_surf.off")
    save_pts_file = os.path.join(root_path, str(index), "tgt_pts.off")

    # save vol
    points = tet_vs
    pad_cell = np.ones([len(tet_elements), 5]) * 4  # 5 v1 v2 v3 v4
    pad_cell[:, 1:] = tet_elements
    vol2off(points, pad_cell, save_vol_file)

    # save surface mesh
    surf_mesh = extract_surface_from_tets(tet_vs, tet_elements)

    pts = surf_mesh.points
    faces = surf_mesh.faces
    faces = faces.reshape(len(faces)//4, 4) #np.ones([len(faces), 4]) * 3  # 4 v1 v2 v3 v4
    #pad_faces[:, 1:] = faces
    mesh2off(pts, faces, save_surf_file)

    # save pts
    pts2off(pts_down_filtered, save_pts_file)

    # #visualize_tet_edges(tet_vs, tet_elements)



