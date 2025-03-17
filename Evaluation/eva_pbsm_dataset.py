import numpy as np
import pyvista as pv


def cal_error(gt, pred, print_error=False):
    diff = np.linalg.norm(gt - pred, axis=1)
    diff_mean = np.mean(diff)
    diff_std = np.std(diff)
    diff_max = np.max(diff)

    diff = pred - gt
    RE = np.sqrt(np.sum(diff * diff) / len(diff))
    if print_error:
        print("mean error: %0.2f, max: %0.2f, std: %0.2f, RE: %0.2f" % (diff_mean, diff_max, diff_std, RE))
    return diff_mean, diff_max, diff_std, RE

def read_off_file(filename):
    vertices = []
    element = []

    with open(filename, 'r') as file:
        lines = file.readlines()

        if not lines[0].strip() == 'OFF':
            raise ValueError("Invalid OFF file format.")

        num_vertices, num_faces, _ = map(int, lines[1].split())

        for line in lines[2:2 + num_vertices]:
            vertex = list(map(float, line.split()))
            vertices.append(vertex)

        for line in lines[2 + num_vertices:]:
            face = list(map(int, line.split()[1:]))
            element.append(face)

    return np.asarray(vertices), element

if __name__ == '__main__':

    index = 4
    root_path = "D:/YZX/BCF_FEM_release/Dataset/PBSM/"
    src_file = root_path + str(index) + "/src_vol_low.off"
    src_def_file = root_path + str(index) + "/Result_ct/deformed_liver_volumetric.off" # change to the right folder, No.4 has two intaoperative surfaces from CT and Stereo
    src_marker_file = root_path + str(index) + "/src_marker.off"
    tgt_marker_file = root_path + str(index) + "/tgt_marker.off"

    src, _ = read_off_file(src_file)
    src_def, _ = read_off_file(src_def_file)
    src_marker, _ = read_off_file(src_marker_file)
    tgt_marker, _ = read_off_file(tgt_marker_file)

    flow = src_def - src
    src_flow_poly = pv.PolyData(src)
    src_flow_poly['x'] = flow[:, 0]
    src_flow_poly['y'] = flow[:, 1]
    src_flow_poly['z'] = flow[:, 2]

    src_marker_poly = pv.PolyData(src_marker)
    src_marker_poly = src_marker_poly.interpolate(src_flow_poly, radius=10, n_points=4)

    inter_flow = np.zeros([len(src_marker_poly.points), 3])
    inter_flow[:, 0] = src_marker_poly['x']
    inter_flow[:, 1] = src_marker_poly['y']
    inter_flow[:, 2] = src_marker_poly['z']

    src_marker_warp = src_marker + inter_flow

    print("Before registration")
    diff_mean_before, _, _, _ = cal_error(tgt_marker, src_marker, print_error=True)
    print("After registration")
    diff_mean, _, _, _ = cal_error(tgt_marker, src_marker_warp, print_error=True)




