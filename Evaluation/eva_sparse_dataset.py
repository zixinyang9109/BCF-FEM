import pyvista as pv
import numpy as np
import os

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

def eva_off_txt(src_file, src_def_file, src_marker_file, tgt_marker_file):

    src, _ = read_off_file(src_file)
    src_def, _ = read_off_file(src_def_file)
    #src_marker = pv.read(src_marker_file).points
    #tgt_marker = pv.read(tgt_marker_file).points

    src_marker = np.loadtxt(src_marker_file)
    tgt_marker = np.loadtxt(tgt_marker_file)

    flow = src_def - src
    src_flow_poly = pv.PolyData(src)
    src_flow_poly['x'] = flow[:, 0]
    src_flow_poly['y'] = flow[:, 1]
    src_flow_poly['z'] = flow[:, 2]

    src_marker_poly = pv.PolyData(src_marker)
    # src_mesh = point_cloud_src_marker.interpolate(point_cloud_src, radius=150)
    src_marker_poly = src_marker_poly.interpolate(src_flow_poly, radius=10, n_points=4)
    #print(len(src_marker_poly.points))

    inter_flow = np.zeros([len(src_marker_poly.points), 3])
    inter_flow[:, 0] = src_marker_poly['x']
    inter_flow[:, 1] = src_marker_poly['y']
    inter_flow[:, 2] = src_marker_poly['z']

    src_marker_warp = src_marker + inter_flow

    diff_mean, diff_max, diff_std, RE = cal_error(tgt_marker, src_marker_warp, print_error=True)

    diff_mean_before, _, _, _ = cal_error(tgt_marker, src_marker, print_error=True)
    print("\n")
    return diff_mean, diff_max, diff_std, RE, diff_mean_before


def get_results(root_path, result_folder, list_file_path):


    src_marker_list = np.loadtxt(list_file_path, dtype=str)[:, 2]
    tgt_marker_list = np.loadtxt(list_file_path, dtype=str)[:, 3]
    # out_stat_path = os.path.join(root_path, method, "stat.txt")

    mean_list = []
    mean_list_before = []

    # with open(out_stat_path, "w") as file:
    for index in np.arange(len(tgt_marker_list)):
        src_file = os.path.join(root_path, "liver_model", "src_vol.off")
        src_def_file = os.path.join(root_path, result_folder, str(index) + "_deformed_liver_volumetric.off")
        src_marker_file = os.path.join(root_path, src_marker_list[index])
        tgt_marker_file = os.path.join(root_path, tgt_marker_list[index])
        diff_mean, diff_max, diff_std, RE, diff_mean_before = eva_off_txt(src_file, src_def_file, src_marker_file, tgt_marker_file)
        mean_list.append(diff_mean)
        mean_list_before.append(diff_mean_before)
    print(result_folder)

            # pred_vs, _ = read_off_file(result_file)
            # diff = np.linalg.norm(gt_vs - pred_vs, axis=1)
            # diff_mean = np.mean(diff)
            # diff_std = np.std(diff)
            # diff_max = np.max(diff)
            #print("mean error: %0.2f, std: %0.2f, max: %0.2f" % (diff_mean, diff_std, diff_max))
            # file.write(" %0.2f $\pm$ %0.2f (%0.2f) \n" % (diff_mean, diff_std, diff_max))
            #file.write(" %0.2f,  %0.2f,  %0.2f \n" % (diff_mean, diff_std, diff_max))

    return mean_list, mean_list_before

def get_noise_vis_stat(mean_list, datasetInfo_fname="D:/YZX/Dataset_all/Sparse_full/datasetInfo.dat"):

    with open(datasetInfo_fname, 'r') as f:
        datasetInfo = [line.split() for line in f]

    #  datasetInfo = [[int(col) for col in row] for row in datasetInfo] # cast strings to int
    #    datasetInfo is Nx3, N = number of target points; [datasetID, extentGroup, deformationGroup]
    #    datasetID is an important column: assumed to be the list of expected results set IDs
    #    extentGroup is assumed to be 1-indexed (i.e. 1: 20-28%, 2: 28-36%, etc.)
    #    deformationGroup is assumed to be 1-indexed (i.e. 1: def A, 2: def B, etc.)

    datasetID = [a[0] for a in datasetInfo]
    extentGroup = np.array([int(a[1]) for a in datasetInfo])

    mean_list = np.array(mean_list)
    group_1 = mean_list[np.where(extentGroup == 1)[0]]
    group_2 = mean_list[np.where(extentGroup == 2)[0]]
    group_3 = mean_list[np.where(extentGroup == 3)[0]]

    group_noise = mean_list[:84]
    group_free = mean_list[84:]

    print("{0:16}   {1:8}   {2:7}   {3}\n".format("Surface Coverage", "Average", "Stdev", "Median"))
    print("{0:16}   {1:8}   {2:7}   {3}\n".format("20-28%", "%.2f" % np.mean(group_1), "%.2f" % np.std(group_1),
                                                  "%.2f" % np.median(group_1)))
    print("{0:16}   {1:8}   {2:7}   {3}\n".format("28-36%", "%.2f" % np.mean(group_2), "%.2f" % np.std(group_2),
                                                  "%.2f" % np.median(group_2)))
    print("{0:16}   {1:8}   {2:7}   {3}\n".format("36-44%", "%.2f" % np.mean(group_3), "%.2f" % np.std(group_3),
                                                  "%.2f" % np.median(group_3)))
    print("{0:16}   {1:8}   {2:7}   {3}\n".format("All Data Sets", "%.2f" % np.mean(mean_list),
                                                  "%.2f" % np.std(mean_list), "%.2f" % np.median(mean_list)))

    print("{0:16}   {1:8}   {2:7}   {3}\n".format("Noise = 2", "%.2f" % np.mean(group_noise),
                                                  "%.2f" % np.std(group_noise), "%.2f" % np.median(group_noise)))

    print("{0:16}   {1:8}   {2:7}   {3}\n".format("Noise = 0", "%.2f" % np.mean(group_free),
                                                  "%.2f" % np.std(group_free), "%.2f" % np.median(group_free)))


if __name__ == '__main__':

    root_path = "D:/YZX/BCF_FEM_release/Dataset/Sparse_full"                           # data root path
    datasetInfo_fname = "D:/YZX/Dataset_all/Sparse_full/datasetInfo.dat"          # dataset info from the organizer
    result_folder = "Result_ICPRigidTransform"                                     # result foler under the root path
    list_file_path = os.path.join(root_path, "ICPRigidTransform_" + "list.txt")   # gt folder under the root path
    mean_list, mean_list_before = get_results(root_path, result_folder, list_file_path)  # print results in terms of visibility ratios, and the entire dataset
    get_noise_vis_stat(mean_list, datasetInfo_fname)                              # print results in terms of noise split



    # get_noise_vis_stat(mean_list)


