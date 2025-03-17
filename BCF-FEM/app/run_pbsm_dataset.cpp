#include "run_registration.h"
using namespace std;

int main() {

	// Replace with the actual file path
	string id = "4";
	string root_path = "D:\\YZX\\BCF_FEM_release\\Dataset\\PBSM\\";

	string fileName_src_vol = root_path + id + "\\src_vol_low.off";
	cout << fileName_src_vol << endl;

	string fileName_force = root_path + id + "\\src_vol_low_surf.off";

	string fileName_vis = root_path + id + "\\src_vol_low_surf.off";

	string fileName_intra = root_path + id + "\\tgt_pts_intra.off";

	//where to save the results
	string pathOut = root_path + id + "\\Result_intra\\";

	test_phantom(pathOut, fileName_src_vol, fileName_force, fileName_vis, fileName_intra);

}