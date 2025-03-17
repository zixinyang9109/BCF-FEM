#include "run_registration.h"
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib> // for std::stoi and std::stof

using namespace std;


int main() {
  
    // Parse required arguments
    string root_path = "D:\\YZX\\BCF_FEM_release\\Dataset\\Demo\\";
    string pathOut = root_path + "Result\\";
    string fileName_src_vol = root_path + "volume.off";
    string fileName_force = root_path + + "volume_surf.off";
    string fileName_vis = root_path + +"volume_surf.off";
    string fileName_intra = root_path + +"tgt_pts.off";

    // Optional arguments with default values
    int nIters = 200;
    float tau =  0.01f;
    float beta =  0.00f;
    float poissonRatio = 0.49f;

    cerr << "Take seconds to initialize.\n";

    // Perform the registration
    test_data(pathOut,
        fileName_src_vol, fileName_force,
        fileName_vis, fileName_intra,
        nIters, tau, beta, poissonRatio);

    return 0;
}
