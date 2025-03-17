#include "run_registration.h"
using namespace std;

int main() {

    
    string root_path = "D:\\YZX\\BCF_FEM_release\\Dataset\\Our_phantoms\\"; // Replace with the actual file path
    string pathSave = root_path + "Result\\";
    std::ifstream file(root_path + "list.txt"); 
   

    if (!file.is_open()) {
        std::cerr << "Error opening the file." << std::endl;

    }

    std::vector<std::string> source_file_paths;
    std::vector<std::string> target_file_paths;
    std::string source_path, target_path, source_marker, target_marker;

    while (file >> source_path >> target_path >> source_marker >> target_marker) {

        source_file_paths.push_back(linuxToWindowsPath(source_path));
        target_file_paths.push_back(linuxToWindowsPath(target_path));
        
    }

    file.close();

    int list_len = source_file_paths.size();
  
    for (size_t i = 0; i < source_file_paths.size(); ++i) {

        string fileName_src_vol = root_path + regex_replace(source_file_paths[i], std::regex("Flat.off"), "volume.off");
        cout << fileName_src_vol << endl;
        string fileName_force = root_path + regex_replace(source_file_paths[i], std::regex("Flat.off"), "volume_surf.off");
        string fileName_vis = root_path + regex_replace(source_file_paths[i], std::regex("Flat.off"), "volume_surf.off");
        string fileName_intra = root_path + target_file_paths[i];
        string pathOut = pathSave + to_string(i) + "_";

        test_phantom(pathOut, fileName_src_vol, fileName_force, fileName_vis, fileName_intra);
    }

    return 0;

}
