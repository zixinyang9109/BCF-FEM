#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <regex>
#include < direct.h >
#include "ForceDrivenRegistration.h"
using namespace std;

std::string linuxToWindowsPath(std::string& linuxPath);

void test_data(string pathOut, string fileName_src_vol,
	string fileName_force, string fileName_vis, string fileName_intra,
	int niters = 200, float tau = 0.01f, float beta = 0.05f, float poissonRatio = 0.49);

void test_phantom(string pathOut, 
	string fileName_src_vol, string fileName_force, 
	string fileName_vis, string fileName_intra);
