#include "run_registration.h"
using namespace std;

std::string linuxToWindowsPath(std::string& linuxPath) {
	std::string windowsPath = linuxPath;
	for (size_t i = 0; i < windowsPath.size(); ++i) {
		if (windowsPath[i] == '/') {
			windowsPath[i] = '\\';
		}
	}
	return windowsPath;
}

void test_data(string pathOut, string fileName_src_vol,
	string fileName_force, string fileName_vis, string fileName_intra,
	int niters, float tau, float beta, float poissonRatio) {

	///////////////////////////////////////////////////////////////////////////////////
	// Save Settings
	FILE* fi;
	string filename_log = pathOut + "log.txt";
	fopen_s(&fi, filename_log.c_str(), "w");
	fprintf_s(fi, " %s  %d \n", "nIters", niters);
	fprintf_s(fi, " %s  %f \n", "tau", tau);
	fprintf_s(fi, " %s  %f \n", "beta", beta);
	fprintf_s(fi, " %s  %f \n", "poissonRatio", poissonRatio);
	fclose(fi);
	///////////////////////////////////////////////////////////////////////////////////



	///////////////////////////////////////////////////////////////////////////////////
	//Get Preop Volumetric Mesh
	//Mesh mesh = get_vol_mesh(fileName_src_vol);
	int nNodes, nElements, junk;
	fopen_s(&fi, fileName_src_vol.c_str(), "r");
	fscanf_s(fi, "OFF\n");
	fscanf_s(fi, "%d %d %d\n", &nNodes, &nElements, &junk);


	float* nodesX = new float[nNodes];
	float* nodesY = new float[nNodes];
	float* nodesZ = new float[nNodes];

	float node1, node2, node3;

	for (int i = 0; i < nNodes; i++)
	{
		fscanf_s(fi, "%f %f %f\n", &node1, &node2, &node3);
		node1 = floor(node1 * 100.0f);
		node1 /= 100.0f;
		node2 = floor(node2 * 100.0f);
		node2 /= 100.0f;
		node3 = floor(node3 * 100.0f);
		node3 /= 100.0f;
		nodesX[i] = node1;
		nodesY[i] = node2;
		nodesZ[i] = node3;
	}

	//fscanf_s(fi, "%d\n", &nElements);

	Element* elements = new Element[nElements];

	int n, v1, v2, v3, v4; // Tetrahedron
	for (int i = 0; i < nElements; i++)
	{
		fscanf_s(fi, "%d %d %d %d %d\n", &n, &v1, &v2, &v3, &v4);
		//change oreder of nodes for the FEM to have positive volume
		elements[i].nodeIDs[0] = v1;
		elements[i].nodeIDs[1] = v2;
		elements[i].nodeIDs[2] = v3;
		elements[i].nodeIDs[3] = v4;
	}

	fclose(fi);

	Mesh mesh;
	mesh.nElements = nElements;
	mesh.nNodes = nNodes;
	mesh.nodesX = nodesX;
	mesh.nodesY = nodesY;
	mesh.nodesZ = nodesZ;
	mesh.elements = elements; //tet
	///////////////////////////////////////////////////////////////////////////////////



	///////////////////////////////////////////////////////////////////////////////////
	// Priors of Force Locations. 
	// Not used, assume forces can be anywhere, and thus the force file is preoperative surface mesh.
	float nodeX, nodeY, nodeZ;
	//int nNodes, nElements, junk;

	fopen_s(&fi, fileName_force.c_str(), "r");

	fscanf_s(fi, "OFF\n");

	fscanf_s(fi, "%d %d %d\n", &nNodes, &nElements, &junk);

	int* forceVertices = new int[nNodes];
	Element* forceElements = new Element[nElements];

	int nForceVertices = nNodes;
	int nForceElements = nElements;


	//read in force vertices and convert to match pre-op mesh node index
	for (int i = 0; i < nForceVertices; i++)
	{
		fscanf_s(fi, "%f %f %f\n", &nodeX, &nodeY, &nodeZ);

		nodeX = floor(nodeX * 100.0f);
		nodeX /= 100.0f;
		nodeY = floor(nodeY * 100.0f);
		nodeY /= 100.0f;
		nodeZ = floor(nodeZ * 100.0f);
		nodeZ /= 100.0f;
		bool bFlag = false;
		//convert node index to mesh node index
		for (int j = 0; j < mesh.nNodes; j++)
		{

			if (nodeX == mesh.nodesX[j] && nodeY == mesh.nodesY[j] && nodeZ == mesh.nodesZ[j])
			{
				forceVertices[i] = j;

				//printf("i : %d  j : %d\n", i, j);
				bFlag = true;
				break;
			}
		}
		if (bFlag == false)
			printf("%d\n", i);

	}


	int vert0, vert1, vert2;

	//get force elements
	for (int i = 0; i < nForceElements; i++)
	{
		fscanf_s(fi, "%d %d %d %d\n", &junk, &vert0, &vert1, &vert2);
		forceElements[i].nodeIDs[0] = forceVertices[vert0];
		forceElements[i].nodeIDs[1] = forceVertices[vert1];
		forceElements[i].nodeIDs[2] = forceVertices[vert2];

		//printf("%d %d %d %d\n", i, forceElements[i].nodeIDs[0], forceElements[i].nodeIDs[1], forceElements[i].nodeIDs[2]);
	}

	fclose(fi);
	///////////////////////////////////////////////////////////////////////////////////



	///////////////////////////////////////////////////////////////////////////////////
	//Prior ZBC, not used.
	// generate array which assigns a force vertice to be a zero boundary vertice
	bool* zeroBoundaryVertices = new bool[nNodes];

	memset(zeroBoundaryVertices, 0, nNodes * sizeof(bool));

	int nZeroBoundaryVertices = 0;

	//// an example to set ZBC
	//for (int i = 0; i < nForceVertices; i++)
	//{

	//	int j = forceVertices[i];

	//	if (mesh.nodesY[j] > 160.0f) ///use for sim1
	//	{
	//		zeroBoundaryVertices[i] = true;
	//		nZeroBoundaryVertices++;
	//	}
	//		
	//	
	//}

	
	///////////////////////////////////////////////////////////////////////////////////



	///////////////////////////////////////////////////////////////////////////////////
	// Priors of visible surface or hard constraints 
	// Not used, and thus set the prior as the entire preoperative surface. 
	fopen_s(&fi, fileName_vis.c_str(), "r");
	fscanf_s(fi, "OFF\n");
	nNodes = 0;
	fscanf_s(fi, "%d %d %d\n", &nNodes, &nElements, &junk);

	Element* MatchingSurfaceElements = new Element[nElements];
	int nMatchingSurfaceElements = nElements;

	int* MatchingSurfaceVertices = new int[nNodes];
	int nMatchingSurfaceVertices = nNodes;

	for (int i = 0; i < nNodes; i++)
	{
		fscanf_s(fi, "%f %f %f\n", &nodeX, &nodeY, &nodeZ);
		nodeX = floor(nodeX * 100.0f);
		nodeX /= 100.0f;
		nodeY = floor(nodeY * 100.0f);
		nodeY /= 100.0f;
		nodeZ = floor(nodeZ * 100.0f);
		nodeZ /= 100.0f;
		//convert node index to mesh node index
		for (int j = 0; j < mesh.nNodes; j++)
		{
			if (nodeX == mesh.nodesX[j] && nodeY == mesh.nodesY[j] && nodeZ == mesh.nodesZ[j])
			{
				MatchingSurfaceVertices[i] = j;
				break;
			}
		}
	}

	//int vert0, vert1, vert2;
	for (int i = 0; i < nElements; i++)
	{
		fscanf_s(fi, "%d %d %d %d\n", &junk, &vert0, &vert1, &vert2);
		MatchingSurfaceElements[i].nodeIDs[0] = MatchingSurfaceVertices[vert0];
		MatchingSurfaceElements[i].nodeIDs[1] = MatchingSurfaceVertices[vert1];
		MatchingSurfaceElements[i].nodeIDs[2] = MatchingSurfaceVertices[vert2];
	}

	fclose(fi);
	///////////////////////////////////////////////////////////////////////////////////



	///////////////////////////////////////////////////////////////////////////////////
	//Get Intraoperative Point Cloud Data : anterior surface of deformed liver
	Pt3D* IntraOpPoints;
	int nIntraOpPoints;

	fopen_s(&fi, fileName_intra.c_str(), "r");
	fscanf_s(fi, "OFF\n");
	nNodes = 0;
	fscanf_s(fi, "%d %d %d\n", &nNodes, &nElements, &junk);

	IntraOpPoints = new Pt3D[nNodes];
	nIntraOpPoints = nNodes;

	for (int i = 0; i < nIntraOpPoints; i++) //to sub sample intraop pts, modify the i++
	{
		fscanf_s(fi, "%f %f %f\n", &nodeX, &nodeY, &nodeZ);

		IntraOpPoints[i].x = nodeX;
		IntraOpPoints[i].y = nodeY;
		IntraOpPoints[i].z = nodeZ;

	}

	fclose(fi);
	///////////////////////////////////////////////////////////////////////////////////



	///////////////////////////////////////////////////////////////////////////////////
	// Preprocessing the point cloud if needed
	//sub sample intraop pts
	float nleafCount = -1;// if you want to yse, set a value.
	bool bFlag = true;
	if (nleafCount != -1.0f)
		downsamplePointCloud(IntraOpPoints, nIntraOpPoints, nleafCount, bFlag);

	//add gaussian noise to point cloud
	srand(time(NULL));
	//srand(0);
	float sigmaZ = -1.0f;
	float sigmaX = -1.0f;
	float sigmaY = -1.0f;
	if (sigmaX > 0.0f || sigmaY > 0.0f || sigmaZ > 0.0f)
		addGaussianNoise(IntraOpPoints, nIntraOpPoints, sigmaX, sigmaY, sigmaZ);

	////apply rigid motion to intraop points
	extrinsicParms parms;
	memset(&parms, 0, sizeof(parms)); // not used here
	rigidMotionTransformation(IntraOpPoints, nIntraOpPoints, parms);
	///////////////////////////////////////////////////////////////////////////////////



	///////////////////////////////////////////////////////////////////////////////////
	// Hard constraints, not used.
	//get groundtruth deformed liver surface from the FEM Simulations of undeformed to deformed
	// and apply the same rigid motion as was applied to intraop surface
	Pt3D* DeformedVolumeGTVertices = NULL;
	int nDeformedVolumeGTVertices = 0;
	///////////////////////////////////////////////////////////////////////////////////



	// Run the Registration.
	for (int i = 0; i < 1; i++)
	{
		ForceDrivenRegistration registration;
		ForceDrivenRegistration::InputData data;

		data.PreOpModel = mesh;

		data.ForceVertices = forceVertices;
		data.nForceVertices = nForceVertices;
		data.ForceElements = forceElements;
		data.nForceElements = nForceElements;

		data.MatchingSurfaceElements = MatchingSurfaceElements;
		data.nMatchingSurfaceElements = nMatchingSurfaceElements;
		data.MatchingSurfaceVertices = MatchingSurfaceVertices;
		data.nMatchingSurfaceVertices = nMatchingSurfaceVertices;

		data.IntraOpPoints = IntraOpPoints;
		data.nIntraOpPoints = nIntraOpPoints;

		data.PreOpFiducials = NULL;
		data.IntraOpFiducials = NULL;
		data.nFiducialPoints = 0;

		data.ZeroBoundaryVertices = zeroBoundaryVertices;
		data.nZeroBoundaryVertices = nZeroBoundaryVertices;
		data.DeformedVolumeGTVertices = DeformedVolumeGTVertices;
		data.nDeformedVolumeGTVertices = nDeformedVolumeGTVertices;
		data.alpha = 0.0f;
		data.beta = beta;
		data.gamma = 0.1f;
		data.tau = tau;

		data.modulus = 1.0; // we do not seek for real forces
		data.poissonRatio = poissonRatio;
		data.bUseNesterov = true;

		registration.setup(data);

		int nIters = niters;
		int nRegularizerIters = 5;
		int nDataIters = 1;

		registration.solve(nIters, nRegularizerIters, nDataIters);

		registration.SaveUpdatedModel_test(fileName_force, pathOut);
	}

	delete[] nodesX;
	delete[] nodesY;
	delete[] nodesZ;
	delete[] elements;
	delete[] forceVertices;
	delete[] forceElements;
	delete[] MatchingSurfaceElements;
	delete[] MatchingSurfaceVertices;
	delete[] zeroBoundaryVertices;
	delete[] IntraOpPoints;
	delete[] DeformedVolumeGTVertices;
}


void test_phantom(string pathOut, 
	string fileName_src_vol, string fileName_force, string fileName_vis, string fileName_intra
) {
	///////////////////////////////////////////////////////////////////////////////////
	//Get Preop Volumetric Mesh
	//Mesh mesh = get_vol_mesh(fileName_src_vol);
	//get preop mesh

	FILE* fi;
	int nNodes, nElements, junk;
	fopen_s(&fi, fileName_src_vol.c_str(), "r");
	fscanf_s(fi, "OFF\n");
	fscanf_s(fi, "%d %d %d\n", &nNodes, &nElements, &junk);


	float* nodesX = new float[nNodes];
	float* nodesY = new float[nNodes];
	float* nodesZ = new float[nNodes];

	float node1, node2, node3;

	for (int i = 0; i < nNodes; i++)
	{
		fscanf_s(fi, "%f %f %f\n", &node1, &node2, &node3);
		node1 = floor(node1 * 100.0f);
		node1 /= 100.0f;
		node2 = floor(node2 * 100.0f);
		node2 /= 100.0f;
		node3 = floor(node3 * 100.0f);
		node3 /= 100.0f;
		nodesX[i] = node1;
		nodesY[i] = node2;
		nodesZ[i] = node3;
	}

	//fscanf_s(fi, "%d\n", &nElements);

	Element* elements = new Element[nElements];

	int n, v1, v2, v3, v4;
	for (int i = 0; i < nElements; i++)
	{
		fscanf_s(fi, "%d %d %d %d %d\n", &n, &v1, &v2, &v3, &v4);
		//change oreder of nodes for the FEM to have positive volume
		elements[i].nodeIDs[0] = v1;
		elements[i].nodeIDs[1] = v2;
		elements[i].nodeIDs[2] = v3;
		elements[i].nodeIDs[3] = v4;
	}

	fclose(fi);

	Mesh mesh;
	mesh.nElements = nElements;
	mesh.nNodes = nNodes;
	mesh.nodesX = nodesX;
	mesh.nodesY = nodesY;
	mesh.nodesZ = nodesZ;
	mesh.elements = elements;
	///////////////////////////////////////////////////////////////////////////////////



	///////////////////////////////////////////////////////////////////////////////////
	// Priors of Force Locations. 
	// Not used, assume forces can be anywhere, and thus the force file is preoperative surface mesh.
	float nodeX, nodeY, nodeZ;
	//int nNodes, nElements, junk;

	fopen_s(&fi, fileName_force.c_str(), "r");

	fscanf_s(fi, "OFF\n");

	fscanf_s(fi, "%d %d %d\n", &nNodes, &nElements, &junk);

	int* forceVertices = new int[nNodes];
	Element* forceElements = new Element[nElements];

	int nForceVertices = nNodes;
	int nForceElements = nElements;


	//read in force vertices and convert to match pre-op mesh node index
	for (int i = 0; i < nForceVertices; i++)
	{
		fscanf_s(fi, "%f %f %f\n", &nodeX, &nodeY, &nodeZ);

		nodeX = floor(nodeX * 100.0f);
		nodeX /= 100.0f;
		nodeY = floor(nodeY * 100.0f);
		nodeY /= 100.0f;
		nodeZ = floor(nodeZ * 100.0f);
		nodeZ /= 100.0f;
		bool bFlag = false;
		//convert node index to mesh node index
		for (int j = 0; j < mesh.nNodes; j++)
		{

			if (nodeX == mesh.nodesX[j] && nodeY == mesh.nodesY[j] && nodeZ == mesh.nodesZ[j])
			{
				forceVertices[i] = j;

				//printf("i : %d  j : %d\n", i, j);
				bFlag = true;
				break;
			}
		}
		if (bFlag == false)
			printf("%d\n", i);

	}


	int vert0, vert1, vert2;

	//get force elements
	for (int i = 0; i < nForceElements; i++)
	{
		fscanf_s(fi, "%d %d %d %d\n", &junk, &vert0, &vert1, &vert2);
		forceElements[i].nodeIDs[0] = forceVertices[vert0];
		forceElements[i].nodeIDs[1] = forceVertices[vert1];
		forceElements[i].nodeIDs[2] = forceVertices[vert2];

		//printf("%d %d %d %d\n", i, forceElements[i].nodeIDs[0], forceElements[i].nodeIDs[1], forceElements[i].nodeIDs[2]);
	}

	fclose(fi);
	///////////////////////////////////////////////////////////////////////////////////



	///////////////////////////////////////////////////////////////////////////////////
	//Prior ZBC, not used.
	//generate array which assigns a force vertice to be a zero boundary vertice
	bool* zeroBoundaryVertices = new bool[nNodes];

	memset(zeroBoundaryVertices, 0, nNodes * sizeof(bool));

	int nZeroBoundaryVertices = 0;

	//// an example to set ZBC
	//for (int i = 0; i < nForceVertices; i++)
	//{

	//	int j = forceVertices[i];

	//	if (mesh.nodesY[j] > 160.0f) ///use for sim1
	//	{
	//		zeroBoundaryVertices[i] = true;
	//		nZeroBoundaryVertices++;
	//	}
	//		
	//	
	//}

	///////////////////////////////////////////////////////////////////////////////////



	///////////////////////////////////////////////////////////////////////////////////
	// Priors of visible preoperative surface.
	// Not used, and thus set the visible surface to be the entire preoperative surface. 
	fopen_s(&fi, fileName_vis.c_str(), "r");

	fscanf_s(fi, "OFF\n");
	nNodes = 0;
	fscanf_s(fi, "%d %d %d\n", &nNodes, &nElements, &junk);

	Element* MatchingSurfaceElements = new Element[nElements];
	int nMatchingSurfaceElements = nElements;

	int* MatchingSurfaceVertices = new int[nNodes];
	int nMatchingSurfaceVertices = nNodes;

	for (int i = 0; i < nNodes; i++) // can be optimized to accelerate
	{
		fscanf_s(fi, "%f %f %f\n", &nodeX, &nodeY, &nodeZ);
		nodeX = floor(nodeX * 100.0f);
		nodeX /= 100.0f;
		nodeY = floor(nodeY * 100.0f);
		nodeY /= 100.0f;
		nodeZ = floor(nodeZ * 100.0f);
		nodeZ /= 100.0f;
		//convert node index to mesh node index
		for (int j = 0; j < mesh.nNodes; j++)
		{
			if (nodeX == mesh.nodesX[j] && nodeY == mesh.nodesY[j] && nodeZ == mesh.nodesZ[j])
			{
				MatchingSurfaceVertices[i] = j;
				break;
			}
		}
	}

	//int vert0, vert1, vert2;

	for (int i = 0; i < nElements; i++)
	{
		fscanf_s(fi, "%d %d %d %d\n", &junk, &vert0, &vert1, &vert2);
		MatchingSurfaceElements[i].nodeIDs[0] = MatchingSurfaceVertices[vert0];
		MatchingSurfaceElements[i].nodeIDs[1] = MatchingSurfaceVertices[vert1];
		MatchingSurfaceElements[i].nodeIDs[2] = MatchingSurfaceVertices[vert2];
	}

	fclose(fi);
	///////////////////////////////////////////////////////////////////////////////////



	///////////////////////////////////////////////////////////////////////////////////
	//Get Intraoperative Point Cloud Data : anterior surface of deformed liver
	Pt3D* IntraOpPoints;
	int nIntraOpPoints;

	fopen_s(&fi, fileName_intra.c_str(), "r");

	fscanf_s(fi, "OFF\n");
	nNodes = 0;
	fscanf_s(fi, "%d %d %d\n", &nNodes, &nElements, &junk);

	IntraOpPoints = new Pt3D[nNodes];
	nIntraOpPoints = nNodes;

	//sub sample intraop pts
	for (int i = 0; i < nIntraOpPoints; i++)
	{
		fscanf_s(fi, "%f %f %f\n", &nodeX, &nodeY, &nodeZ);

		IntraOpPoints[i].x = nodeX;
		IntraOpPoints[i].y = nodeY;
		IntraOpPoints[i].z = nodeZ;

	}

	fclose(fi);
	///////////////////////////////////////////////////////////////////////////////////



	///////////////////////////////////////////////////////////////////////////////////
	// Preprocessing the point cloud if needed
	//sub sample intraop pts
	float nleafCount = -1;// -1.0f;
	bool bFlag = true;
	if (nleafCount != -1.0f)
		downsamplePointCloud(IntraOpPoints, nIntraOpPoints, nleafCount, bFlag);

	//add gaussian noise to point cloud
	srand(time(NULL));
	//srand(0);
	float sigmaZ = -1.0f;
	float sigmaX = -1.0f;
	float sigmaY = -1.0f;
	if (sigmaX > 0.0f || sigmaY > 0.0f || sigmaZ > 0.0f)
		addGaussianNoise(IntraOpPoints, nIntraOpPoints, sigmaX, sigmaY, sigmaZ);

	////apply rigid motion to intraop points
	extrinsicParms parms;
	memset(&parms, 0, sizeof(parms));// not used here
	rigidMotionTransformation(IntraOpPoints, nIntraOpPoints, parms);
	///////////////////////////////////////////////////////////////////////////////////



	///////////////////////////////////////////////////////////////////////////////////
	// Hard constraints, not used.
	//get groundtruth deformed liver surface from the FEM Simulations of undeformed to deformed
	// and apply the same rigid motion as was applied to intraop surface
	Pt3D* DeformedVolumeGTVertices = NULL;
	int nDeformedVolumeGTVertices = 0;
	///////////////////////////////////////////////////////////////////////////////////

	// Run the Registration.
	for (int i = 0; i < 1; i++)
	{

		ForceDrivenRegistration registration;
		ForceDrivenRegistration::InputData data;

		data.PreOpModel = mesh;
		data.ForceVertices = forceVertices;
		data.nForceVertices = nForceVertices;
		data.ForceElements = forceElements;
		data.nForceElements = nForceElements;

		data.MatchingSurfaceElements = MatchingSurfaceElements;
		data.nMatchingSurfaceElements = nMatchingSurfaceElements;
		data.MatchingSurfaceVertices = MatchingSurfaceVertices;
		data.nMatchingSurfaceVertices = nMatchingSurfaceVertices;

		data.IntraOpPoints = IntraOpPoints;
		data.nIntraOpPoints = nIntraOpPoints;

		data.PreOpFiducials = NULL;
		data.IntraOpFiducials = NULL;
		data.nFiducialPoints = 0;

		data.ZeroBoundaryVertices = zeroBoundaryVertices;
		data.nZeroBoundaryVertices = nZeroBoundaryVertices;

		data.DeformedVolumeGTVertices = DeformedVolumeGTVertices;
		data.nDeformedVolumeGTVertices = nDeformedVolumeGTVertices;

		data.alpha = 0.0f;
		data.beta = 0.00f; // 0.05
		data.gamma = 0.1f;
		data.tau = 0.01f;

		data.modulus = 1.0;
		data.poissonRatio = 0.49;
		data.bUseNesterov = true;

		registration.setup(data);

		int nIters = 200;
		int nRegularizerIters = 5;
		int nDataIters = 1;

		registration.solve(nIters, nRegularizerIters, nDataIters);

		registration.SaveUpdatedModel_test(fileName_force, pathOut);


	}

	delete[] nodesX;
	delete[] nodesY;
	delete[] nodesZ;
	delete[] elements;
	delete[] forceVertices;
	delete[] forceElements;
	delete[] MatchingSurfaceElements;
	delete[] MatchingSurfaceVertices;
	delete[] zeroBoundaryVertices;
	delete[] IntraOpPoints;
	delete[] DeformedVolumeGTVertices;
}


