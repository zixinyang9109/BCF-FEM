#include "run_registration.h"
using namespace std;

int main() {

	// Replace with the actual file path
	string root_path =  "D:\\YZX\\BCF_FEM_release\\Dataset\\Sparse_full\\";  //"C:\\Users\\yzxar\\code\\Project\\Non_rigid_registration\\Dataset_all\\Sparse_full\\";
	int nIters = 200;
	float tau = 0.01f; // kss
	float beta = 0.00f; // no regularization, if wwant to use, suggested to set 0.05. It is designed for high noise situation. The detail can be seen paperV1.
	float poissonRatio = 0.49;
	
	vector<string> folders = { "ICPRigidTransform"}; //"ICPRigidTransform",  "OptimalRigidTransform", "wICPRigidTransform"

	for (string folder : folders)
	{
		string list_file = root_path + folder + "_list_off.txt";  // Replace with the actual file path

		string save_path = root_path + "Result_" + folder + "\\"; // without regularization beta = 0

		std::ifstream file(list_file); // Replace with the actual file path

		if (!file.is_open()) {
			std::cerr << "Error opening the file." << std::endl;

		}

		std::vector<std::string> source_file_paths;
		std::vector<std::string> target_file_paths;
		std::vector<std::string> source_vol_paths;
		std::string source_path, target_path, source_vol, target_marker;


		while (file >> source_path >> target_path >> source_vol >> target_marker) {

			source_file_paths.push_back(linuxToWindowsPath(source_path));
			target_file_paths.push_back(linuxToWindowsPath(target_path));
			source_vol_paths.push_back(linuxToWindowsPath(source_vol));
			
		}

		file.close();

		int list_len = source_file_paths.size();
	

		for (size_t i = 0; i < source_file_paths.size(); ++i) {


			string fileName_src_vol = root_path + source_vol_paths[i];

			cout << fileName_src_vol << endl;

			string fileName_force = root_path + source_file_paths[i];

			string fileName_vis = fileName_force;

			string fileName_intra = root_path + target_file_paths[i];

			//where to save the results
			string pathOut = save_path + to_string(i) + "_"; 

			test_data(pathOut, fileName_src_vol, fileName_force, fileName_vis, fileName_intra, nIters, tau, beta, poissonRatio);
		}

	}

	return 0;

}




