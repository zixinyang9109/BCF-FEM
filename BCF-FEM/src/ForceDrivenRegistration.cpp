#define _USE_MATH_DEFINES
#include "ForceDrivenRegistration.h"


inline float clamp(double x, double upper, double lower)
{
	return (float)min(upper, max(x, lower));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Original implementation : iterative solution
//ClosestPt operator using OpenMP 
//optimization :  Nesterov augmented gradient algorithm
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ForceDrivenRegistration::ForceDrivenRegistration()
{
	m_modulus = 1.0;
	m_poissonRatio = 0.49;
	m_bLinearElasticMaterial = true;

	m_Neighbors = NULL;
	m_nNeighbors = NULL;

	m_ForceVertices = NULL;
	m_MatchingSurfaceElements = NULL;
	m_IntraOpPoints = NULL;

	m_deformedPts = NULL;
	m_undeformedPts = NULL;

	m_DataTermErr = NULL;
	m_RegistrationErr = NULL;

	m_dist = NULL;
    m_baryCentric = NULL;
	m_nodeIDs = NULL;

	m_lambda = 1.0;

	m_sigma2 = 1.0;

	m_alpha = m_beta = m_gamma = 0.0;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ForceDrivenRegistration::~ForceDrivenRegistration()
{
	FreeImage(m_Neighbors);
	delete[] m_nNeighbors;

	delete[] m_deformedPts;
	delete[] m_undeformedPts;

	delete[] m_DataTermErr;
	delete[] m_RegistrationErr;

	delete[] m_dist;
	delete[] m_baryCentric;
	delete[] m_nodeIDs;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ForceDrivenRegistration::setup(InputData data)
{
	//set up data
	m_PreOpModel = data.PreOpModel;
	m_ForceVertices = data.ForceVertices;
	m_ForceElements = data.ForceElements;
	m_ZeroBoundaryVertices = data.ZeroBoundaryVertices;
	m_MatchingSurfaceElements = data.MatchingSurfaceElements;
	m_MatchingSurfaceVertices = data.MatchingSurfaceVertices;
	m_IntraOpPoints = data.IntraOpPoints;
	m_DeformedVolumeGTVertices = data.DeformedVolumeGTVertices;

	m_nForceVertices = data.nForceVertices;
	m_nForceElements = data.nForceElements;
	m_nZeroBoundaryVertices = data.nZeroBoundaryVertices;
	m_nMatchingSurfaceElements = data.nMatchingSurfaceElements;
	m_nMatchingSurfaceVertices = data.nMatchingSurfaceVertices;
	m_nIntraOpPoints = data.nIntraOpPoints;
	m_nDeformedVolumeGTVertices = data.nDeformedVolumeGTVertices;

	m_alpha = data.alpha;
	m_beta = data.beta;
	m_gamma = data.gamma;
	m_tau = data.tau;
	m_modulus = data.modulus; // 1.0;
	m_poissonRatio = data.poissonRatio;//0.49;
	m_bUseNesterov = data.bUseNesterov;

	m_FEM.setup(m_PreOpModel, (float)m_modulus, (float)m_poissonRatio, m_bLinearElasticMaterial);

	size_t nNodes = 3 * size_t(m_PreOpModel.nNodes);

	//allocate displacement array
	m_u.resize(nNodes);
	m_u.setZero();

	//allocate force array
	m_f.resize(nNodes);
	m_f.setZero();

	//allocate intermediate force array
	m_fBar.resize(nNodes);
	m_fBar.setZero();

	//allocate grad(= dDataTerm / dDisplacement = dDataTerm / dClosestPt  * dClosestPt / dDisplacement) array
	m_grad.resize(nNodes);
	m_grad.setZero();

	//allocate dDataTerm / dforce array
	m_dJdf.resize(nNodes);
	m_dJdf.setZero();

	//allocate array for current node positions
	m_updateNodePositions.resize(nNodes);
	m_updateNodePositions.setZero();		

	//calculate initial stiffness matrix
	m_FEM.calculateStiffnessMatrix(m_u, m_K);

	//set zero boundary conditions using penalty method by modifying stiffness matrix
	if (m_nZeroBoundaryVertices > 0)
		SetZeroBoundaryConditions();
	else
		SetInitialConditions();

	//Factorization of stiffness matrix
	m_solver.compute(m_K);

	//generate adjacency matrix for regularizer
	Adjacency_Matrix();

	//allocate arrays for rigid motion correction
	m_deformedPts = new Pt3D[size_t(m_PreOpModel.nNodes)];
	m_undeformedPts = new Pt3D[size_t(m_PreOpModel.nNodes)];

	//allocating arrays for updating alpha
	m_dist = new double[3 * m_nIntraOpPoints];
	memset(m_dist, 0, 3 * sizeof(double));

	m_baryCentric = new float[3 * m_nIntraOpPoints];
	memset(m_baryCentric, 0, 3 * sizeof(float));

	m_nodeIDs = new int[3 * m_nIntraOpPoints];
	memset(m_nodeIDs, 0, 3 * sizeof(int));

}


/////////////////////////////////////////////////////////////////////////////////////////
// //
// /////////////////////////////////////////////////////////////////////////////////////////
void ForceDrivenRegistration::SaveUpdatedModel_dev(string pathIn, string pathOut)
{
	//for (int i = 0; i < m_nForceVertices; i++)
	//{
	//	int j = m_ForceVertices[i];
	//	float y = m_PreOpModel.nodesY[j];
	//	float force = -0.0027f * y + 0.432f;
	//	//if (y >= 150.0f)
	//	//	force = -0.0017f * y + 0.255f;
	//	//else if (y <= 90.0f)
	//	//	force = -0.0037 * y + 0.333f;
	//	//else
	//	//	force = 0.0f;
	//	//force *= 4.0f;
	//	if (force < 0.0f)
	//		force = 0.0f;
	//	m_f(3*j + 2) = force / 2.0f;
	//}

	m_u = m_solver.solve(m_f);

	//Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver0;
	//solver0.compute(m_K0);

	//Eigen::VectorXd f0 = m_K0 * m_u;

	////Eigen::VectorXf f00 = m_f - 0.01f * m_u;
	////Eigen::VectorXf f00 = 0.01f * m_u;
	//for (int i = 0; i < 3 * m_PreOpModel.nNodes; i++)
	//{
	//	if (m_f(i) == 0.0f)
	//	{
	//		f0(i) = 0.0f;
	//	}
	//}

	////Eigen::VectorXf u = m_solver.solve(f00);
	//Eigen::VectorXd u0 = solver0.solve(f0);

	//Eigen::VectorXd u00 = m_u - u0;

	//m_u = u00;
	//FILE* fi0;

	//string fileName0 = pathOut + "forces.txt";
	//fopen_s(&fi0, fileName0.c_str(), "w");
	//for (int i = 0; i < 3 * m_PreOpModel.nNodes; i++)
	//	fprintf_s(fi0, "%f\t%f\t%f\t%f\t%f\n", m_f(i), f0(i), m_u(i),u0(i),u00(i));

	//fclose(fi0);

	//for(int i = 0; i < 3*m_PreOpModel.nNodes; i++)
	//	m_u(i) = u0(i);

	//m_u = solver0.solve(f0);

	//path = "D:\\LiverMeshesDatabase\\Simulations\\KellysPhantomLiver\\Sim1\\";
	string path = pathIn;

	//get undeformed liver surface and generate deformed surface model
	FILE* fi;
	string fileName = path + "undeformed_liver_surface.off";
	fopen_s(&fi, fileName.c_str(), "r");

	fscanf_s(fi, "OFF\n");
	int nNodes, nElements, junk;
	fscanf_s(fi, "%d %d %d\n", &nNodes, &nElements, &junk);

	Element* SurfaceElements = new Element[nElements];
	int nSurfaceElements = nElements;

	int* SurfaceVertices = new int[nNodes];
	int nSurfaceVertices = nNodes;

	//convert vertices to match preop mesh node index
	float nodeX, nodeY, nodeZ;
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
		for (int j = 0; j < nNodes; j++)
		{
			if (nodeX == m_PreOpModel.nodesX[j] && nodeY == m_PreOpModel.nodesY[j] && nodeZ == m_PreOpModel.nodesZ[j])
			{
				SurfaceVertices[i] = j;
				break;
			}
		}
	}

	int vert0, vert1, vert2;

	for (int i = 0; i < nElements; i++)
	{
		fscanf_s(fi, "%d %d %d %d\n", &junk, &vert0, &vert1, &vert2);
		SurfaceElements[i].nodeIDs[0] = vert0;
		SurfaceElements[i].nodeIDs[1] = vert1;
		SurfaceElements[i].nodeIDs[2] = vert2;
	}

	fclose(fi);

	for (size_t i = 0; i < nNodes; i++)
	{
		int j = SurfaceVertices[i];
		m_updateNodePositions[3 * i] = m_PreOpModel.nodesX[j] + (float)m_u[3 * j];
		m_updateNodePositions[3 * i + 1] = m_PreOpModel.nodesY[j] + (float)m_u[3 * j + 1];
		m_updateNodePositions[3 * i + 2] = m_PreOpModel.nodesZ[j] + (float)m_u[3 * j + 2];
	}

	//location to save data 
	//path = "D:\\LiverMeshesDatabase\\Simulations\\KellysPhantomLiver\\Sim1\\Results\\26\\";
	path = pathOut;

	//save deformed liver surface
	string filename = path + "deformed_liver_surface.off";
	//FILE* fi;
	fopen_s(&fi, filename.c_str(), "w");
	fprintf_s(fi, "OFF\n");
	fprintf_s(fi, "%d %d %d\n", nNodes, nElements, 0);

	for (int i = 0; i < nNodes; i++)
	{
		fprintf_s(fi, "%f %f %f\n", m_updateNodePositions[3 * i], m_updateNodePositions[3 * i + 1], m_updateNodePositions[3 * i + 2]);
	}

	for (int i = 0; i < nElements; i++)

		fprintf_s(fi, "%d %d %d %d\n", 3, SurfaceElements[i].nodeIDs[0], SurfaceElements[i].nodeIDs[1],
			SurfaceElements[i].nodeIDs[2]);

	fclose(fi);

	if (m_nDeformedVolumeGTVertices > 0)
	{
		//save ground truth deformed liver surface
		filename = path + "gt_deformed_liver_surface.off";
		//FILE* fi;
		fopen_s(&fi, filename.c_str(), "w");
		fprintf_s(fi, "OFF\n");
		fprintf_s(fi, "%d %d %d\n", nNodes, nElements, 0);

		for (int i = 0; i < nNodes; i++)
		{
			fprintf_s(fi, "%f %f %f\n", m_DeformedVolumeGTVertices[i].x, m_DeformedVolumeGTVertices[i].y, m_DeformedVolumeGTVertices[i].z);
		}

		for (int i = 0; i < nElements; i++)

			fprintf_s(fi, "%d %d %d %d\n", 3, SurfaceElements[i].nodeIDs[0], SurfaceElements[i].nodeIDs[1],
				SurfaceElements[i].nodeIDs[2]);

		fclose(fi);
	}
	//save deformed liver posterior surface
	filename = path + "deformed_liver_posterior_surface.off";
	//FILE* fi;
	fopen_s(&fi, filename.c_str(), "w");
	fprintf_s(fi, "OFF\n");
	fprintf_s(fi, "%d %d %d\n", m_nForceVertices, m_nForceElements, 0);

	for (int j = 0; j < m_nForceVertices; j++)
	{
		int i = m_ForceVertices[j];
		fprintf_s(fi, "%f %f %f\n", m_updateNodePositions[3 * i], m_updateNodePositions[3 * i + 1], m_updateNodePositions[3 * i + 2]);
	}

	fclose(fi);

	filename = path + "deformed_liver_anterior_surface.off";
	//FILE* fi;
	fopen_s(&fi, filename.c_str(), "w");
	fprintf_s(fi, "OFF\n");
	fprintf_s(fi, "%d %d %d\n", m_nMatchingSurfaceVertices, m_nMatchingSurfaceElements, 0);

	for (int j = 0; j < m_nMatchingSurfaceVertices; j++)
	{
		int i = m_MatchingSurfaceVertices[j];
		fprintf_s(fi, "%f %f %f\n", m_updateNodePositions[3 * i], m_updateNodePositions[3 * i + 1], m_updateNodePositions[3 * i + 2]);
	}

	fclose(fi);

	filename = path + "ForceVectors.txt";
	fopen_s(&fi, filename.c_str(), "w");
	fprintf(fi, "%d\n", m_nForceVertices);
	for (int i = 0; i < m_nForceVertices; i++)
	{
		int j = m_ForceVertices[i];
		//fprintf_s(fi, "%d %f %f %f %f %f %f\n", j, m_updateNodePositions[3 * j], m_updateNodePositions[3 * j + 1], m_updateNodePositions[3 * j + 2],
		//	                         m_f(3 * j), m_f(3 * j + 1), m_f(3 * j + 2));
		float mag = (float)sqrt(m_f(3 * j) * m_f(3 * j) + m_f(3 * j + 1) * m_f(3 * j + 1) + m_f(3 * j + 2) * m_f(3 * j + 2));
		if (m_f(3 * j + 2) < 0.0f)
			mag *= -1.0f;
		fprintf_s(fi, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n", m_PreOpModel.nodesX[j], m_PreOpModel.nodesY[j], m_PreOpModel.nodesZ[j],
			m_f(3 * j), m_f(3 * j + 1), m_f(3 * j + 2), mag);
	}
	fclose(fi);

	filename = path + "SurfaceForceVectors.txt";
	fopen_s(&fi, filename.c_str(), "w");
	fprintf(fi, "%d\n", nNodes);
	for (int i = 0; i < nNodes; i++)
	{
		int j = SurfaceVertices[i];
		float mag = (float)sqrt(m_f(3 * j) * m_f(3 * j) + m_f(3 * j + 1) * m_f(3 * j + 1) + m_f(3 * j + 2) * m_f(3 * j + 2));
		if (m_f(3 * j + 2) < 0.0f)
			mag *= -1.0f;
		fprintf_s(fi, "%f\n", mag);
	}
	fclose(fi);

	delete[]SurfaceElements;
	delete[]SurfaceVertices;

	//save deformed liver volume
	for (size_t i = 0; i < m_PreOpModel.nNodes; i++)
	{
		m_updateNodePositions[3 * i] = m_PreOpModel.nodesX[i] + (float)m_u[3 * i];
		m_updateNodePositions[3 * i + 1] = m_PreOpModel.nodesY[i] + (float)m_u[3 * i + 1];
		m_updateNodePositions[3 * i + 2] = m_PreOpModel.nodesZ[i] + (float)m_u[3 * i + 2];
	}

	filename = path + "deformed_liver_volumetric.off";
	//FILE* fi;
	fopen_s(&fi, filename.c_str(), "w");
	fprintf_s(fi, "OFF\n");
	fprintf_s(fi, "%d %d %d\n", m_PreOpModel.nNodes, m_PreOpModel.nElements, 0);

	for (int i = 0; i < m_PreOpModel.nNodes; i++)
	{
		fprintf_s(fi, "%f %f %f\n", m_updateNodePositions[3 * i], m_updateNodePositions[3 * i + 1], m_updateNodePositions[3 * i + 2]);
	}

	for (int i = 0; i < m_PreOpModel.nElements; i++)

		fprintf_s(fi, "%d %d %d %d %d\n", 4, m_PreOpModel.elements[i].nodeIDs[0], m_PreOpModel.elements[i].nodeIDs[1],
			m_PreOpModel.elements[i].nodeIDs[2], m_PreOpModel.elements[i].nodeIDs[3]);

	fclose(fi);

	if (m_nDeformedVolumeGTVertices > 0)
	{
		filename = path + "gt_deformed_liver_volumetric.off";
		//FILE* fi;
		fopen_s(&fi, filename.c_str(), "w");
		fprintf_s(fi, "OFF\n");
		fprintf_s(fi, "%d %d %d\n", m_PreOpModel.nNodes, m_PreOpModel.nElements, 0);

		for (int i = 0; i < m_PreOpModel.nElements; i++)

			fprintf_s(fi, "%d %d %d %d %d\n", 4, m_PreOpModel.elements[i].nodeIDs[0], m_PreOpModel.elements[i].nodeIDs[1],
				m_PreOpModel.elements[i].nodeIDs[2], m_PreOpModel.elements[i].nodeIDs[3]);

		fclose(fi);
	}

	//save intraoperative point cloud
	filename = path + "intraoperative_point_cloud.off";
	//FILE* fi;
	fopen_s(&fi, filename.c_str(), "w");
	fprintf_s(fi, "OFF\n");
	fprintf_s(fi, "%d %d %d\n", m_nIntraOpPoints, 0, 0);
	for (int i = 0; i < m_nIntraOpPoints; i++)
	{
		fprintf_s(fi, "%f %f %f\n", m_IntraOpPoints[i].x, m_IntraOpPoints[i].y, m_IntraOpPoints[i].z);
	}

	fclose(fi);

	//save displacemnt vecotrs for dirchlet vertices
	filename = path + "DirchletVectors.txt";
	fopen_s(&fi, filename.c_str(), "w");
	fprintf(fi, "%d\n", m_nForceVertices);
	for (int i = 0; i < m_nForceVertices; i++)
	{
		int j = m_ForceVertices[i];

		float mag = (float)sqrt(m_u(3 * j) * m_u(3 * j) + m_u(3 * j + 1) * m_u(3 * j + 1) + m_u(3 * j + 2) * m_u(3 * j + 2));
		fprintf_s(fi, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n", m_PreOpModel.nodesX[j], m_PreOpModel.nodesY[j], m_PreOpModel.nodesZ[j],
			m_u(3 * j), m_u(3 * j + 1), m_u(3 * j + 2), mag);
	}
	fclose(fi);

	//save registration err vs iters
	filename = path + "RegErrs_vs_Iters.txt";
	fopen_s(&fi, filename.c_str(), "w");
	fprintf(fi, "%d\n", m_nAlgorithmIters);
	for (int i = 0; i < m_nAlgorithmIters; i++)
	{
		fprintf_s(fi, "%d\t%f\t%f\n", i, m_DataTermErr[i], m_RegistrationErr[i]);
	}
	fclose(fi);

	//save registration error vs node
	if (m_nDeformedVolumeGTVertices > 0)
	{
		float avg = 0.0f;
		float avg2 = 0.0f;
		float max = FLT_MIN;
		filename = path + "RegistationErrors1.txt";
		fopen_s(&fi, filename.c_str(), "w");
		fprintf(fi, "%d\n", m_nDeformedVolumeGTVertices);
		for (int i = 0; i < m_nDeformedVolumeGTVertices; i++)
		{
			float delX = m_PreOpModel.nodesX[i] + (float)m_u[3 * i] - m_DeformedVolumeGTVertices[i].x;
			float delY = m_PreOpModel.nodesY[i] + (float)m_u[3 * i + 1] - m_DeformedVolumeGTVertices[i].y;
			float delZ = m_PreOpModel.nodesZ[i] + (float)m_u[3 * i + 2] - m_DeformedVolumeGTVertices[i].z;
			float err = delX * delX + delY * delY + delZ * delZ;
			err = sqrt(err);
			avg += err;
			avg2 += err * err;
			if (err > max)
				max = err;
			fprintf(fi, "%f\t%f\t%f\t%f\n", delX, delY, delZ, err);
		}

		

		avg /= (float)m_nDeformedVolumeGTVertices;
		avg2 /= (float)m_nDeformedVolumeGTVertices;

		float sd = sqrt(avg2 - avg * avg);

		printf("avg : %f sd : %f max : %f\n", avg, sd, max);
		fprintf(fi, "avg : %f sd : %f max : %f\n", avg, sd, max);

		fclose(fi);
	}

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//This is only good for zero displacement boundary conditions : Penalty method
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ForceDrivenRegistration::GroundTruthError(double &avg, double &sd, double &max, bool bRMS)
{
	double err = 0.0;
	avg = 0.0f;
	double avg2 = 0.0f;
	max = DBL_MIN;
	if (m_nDeformedVolumeGTVertices > 0)
	{
		for (int i = 0; i < m_nDeformedVolumeGTVertices; i++)
		{
			double delX = m_PreOpModel.nodesX[i] + m_u[3 * i] - m_DeformedVolumeGTVertices[i].x;
			double delY = m_PreOpModel.nodesY[i] + m_u[3 * i + 1] - m_DeformedVolumeGTVertices[i].y;
			double delZ = m_PreOpModel.nodesZ[i] + m_u[3 * i + 2] - m_DeformedVolumeGTVertices[i].z;

			err = delX * delX + delY * delY + delZ * delZ;
			if(!bRMS)
				err = sqrt(err);
			avg += err;
			avg2 += err * err;
			if (err > max)
				max = err;
		}

		avg /= (double)m_nDeformedVolumeGTVertices;
		avg2 /= (double)m_nDeformedVolumeGTVertices;

		sd = sqrt(avg2 - avg * avg);

	}

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//This is only good for zero displacement boundary conditions : Penalty method
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ForceDrivenRegistration::SetZeroBoundaryConditions()
{
	int nDirichletConstraints = 3 * m_nZeroBoundaryVertices;
	Constraint* DirichletConstraints = new Constraint[nDirichletConstraints];

	int index = 0;
	for (int i = 0; i < m_nForceVertices; i++)
	{
		if (m_ZeroBoundaryVertices[i] == true)
		{
			DirichletConstraints[index].node = m_ForceVertices[i];
			DirichletConstraints[index].DOF = 0;
			index++;
			DirichletConstraints[index].node = m_ForceVertices[i];
			DirichletConstraints[index].DOF = 1;
			index++;
			DirichletConstraints[index].node = m_ForceVertices[i];
			DirichletConstraints[index].DOF = 2;
			index++;
		}
	}

	//find maximum value of stiffness matrix
	double max = 0.0f;
	
	for (int k = 0; k < m_K.outerSize(); k++)
	{
		for (Eigen::SparseMatrix<double>::InnerIterator it(m_K, k); it; ++it)
		{

			if (it.valueRef() > max)
				max = it.valueRef();

		}
	}

	//add penalty to stiffness matrix
	double penalty = max * 10000.0;
	
	int nDims = 3;

	int nConstraints = 0;
	for (int k = 0; k < m_K.outerSize(); k++)
	{
		for (Eigen::SparseMatrix<double>::InnerIterator it(m_K, k); it; ++it)
		{
			for (int i = 0; i < nDirichletConstraints; i++)
			{

				int index = nDims * DirichletConstraints[i].node + DirichletConstraints[i].DOF;

				if (it.row() == index || it.col() == index)
				{
					if (it.row() == it.col())
					{
						it.valueRef() = penalty;
					}
				}

			}
		}
	}

	delete[] DirichletConstraints;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Set
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ForceDrivenRegistration::SetInitialConditions()
{
	double penalty = m_tau;

	int nConstraints = 3 * m_nForceVertices;
	Constraint* Constraints = new Constraint[nConstraints];

	int index = 0;
	for (int i = 0; i < m_nForceVertices; i++)
	{
		if (m_ZeroBoundaryVertices[i] == false)
		{
			Constraints[index].node = m_ForceVertices[i];
			Constraints[index].DOF = 0;
			index++;
			Constraints[index].node = m_ForceVertices[i];
			Constraints[index].DOF = 1;
			index++;
			Constraints[index].node = m_ForceVertices[i];
			Constraints[index].DOF = 2;
			index++;
		}
	}

	/*int nConstraints = 3 * m_PreOpModel.nNodes;
	Constraint* Constraints = new Constraint[nConstraints];

	int index = 0;
	for (int i = 0; i < m_PreOpModel.nNodes; i++)
	{
		if (m_ZeroBoundaryVertices[i] == false)
		{
			Constraints[index].node = i;
			Constraints[index].DOF = 0;
			index++;
			Constraints[index].node = i;
			Constraints[index].DOF = 1;
			index++;
			Constraints[index].node = i;
			Constraints[index].DOF = 2;
			index++;
		}
	}*/

	nConstraints = index;

	int nDims = 3;

	for (int k = 0; k < m_K.outerSize(); k++)
	{
		for (Eigen::SparseMatrix<double>::InnerIterator it(m_K, k); it; ++it)
		{
			for (int i = 0; i < nConstraints; i++)
			{

				int index = nDims * Constraints[i].node + Constraints[i].DOF;

				if (it.row() == index || it.col() == index)
				{
					if (it.row() == it.col())
					{
						it.valueRef() += penalty;
					}
				}

			}
		}
	}

	delete[] Constraints;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Use this for the force vector smoothness regularizer
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ForceDrivenRegistration::Adjacency_Matrix()
{
	size_t nNodes = (size_t)m_PreOpModel.nNodes;

	//create node adjacency matrix
	m_Neighbors = AllocateImage<int>(10, (int)nNodes);
	m_nNeighbors = new int[nNodes];

	memset(m_nNeighbors, 0, nNodes * sizeof(int));
	memset(m_Neighbors[0], 0, nNodes * 10 * sizeof(int));


	for (int i = 0; i < m_nForceElements; i++)
	{

		for (int j = 0; j < 3; j++)
		{
			int nodej = m_ForceElements[i].nodeIDs[j];

			for (int jj = 0; jj < 3; jj++)
			{
				if (j != jj)
				{
					bool bInList = false;
					int nodejj = m_ForceElements[i].nodeIDs[jj];

					for (int k = 0; k < m_nNeighbors[nodej]; k++)
					{
						int nodek = m_Neighbors[nodej][k];
						if (nodek == nodejj)
							bInList = true;
					}

					if (!bInList)
					{
						m_Neighbors[nodej][m_nNeighbors[nodej]] = nodejj;
						m_nNeighbors[nodej]++;
					}
				}

			}

		}
	} //elements

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ForceDrivenRegistration::solve(int nIters, int nRegularizerIters, int nDataIters)
{
	double DataTerm;
	double gamma = 0.0f;

	double alp = 0.05f; 

	//allocate arrays used for Nesterov accelerated gradient algorithm
	Eigen::VectorXd p;
	Eigen::VectorXd fOld;

	p.resize(3 * size_t(m_PreOpModel.nNodes)); // updated forces with momentum
	p.setZero();

	fOld.resize(3 * size_t(m_PreOpModel.nNodes)); // previous forces
	fOld.setZero();

	m_DataTermErr = new double[nIters];
	m_RegistrationErr = new double[nIters];

	memset(m_DataTermErr, 0, nIters * sizeof(double));
	memset(m_RegistrationErr, 0, nIters * sizeof(double));

	m_nAlgorithmIters = nIters;

	
	double err = 0;
	//bool bRMS = false;
	// double err, errSD, errMax;
	//GroundTruthError(err, errSD, errMax, bRMS); // for dev
	m_RegistrationErr[0] = err;// sqrt(err);

	auto start = std::chrono::high_resolution_clock::now();

	for (int iters = 0; iters < nIters; iters++)
	{
		//GroundTruthError(err, errSD, errMax, bRMS); //for dev

		m_RegistrationErr[iters] = err;// sqrt(err);

		int iters3 = iters + 3;
		gamma = (double)iters / (double)(iters3);

		if(!m_bUseNesterov)
		   gamma = 0.0f; //set gamma to zero to turn off Nesterov acceleration

		p = m_f + gamma * (m_f - fOld);

		fOld = m_f;

		gradientDataTerm(p, m_dJdf, DataTerm); // use forces to update displacements and update correspondences

		updateAlpha(alp); //update alpha 

		//gradientDataTermEM(p, m_dJdf, DataTerm);
		// 
		// update new forces
		for (int i = 0; i < m_nForceVertices; i++)
		{
			if (m_ZeroBoundaryVertices[i] == false)
			{
				size_t vertID = (size_t)m_ForceVertices[i];
				m_f[3 * vertID] = p[3 * vertID] - alp * m_dJdf[3 * vertID];
				m_f[3 * vertID + 1] = p[3 * vertID + 1] - alp * m_dJdf[3 * vertID + 1];
				m_f[3 * vertID + 2] = p[3 * vertID + 2] - alp * m_dJdf[3 * vertID + 2] ;
			}
		}

		m_DataTermErr[iters] = DataTerm;

		Regularizer((float)alp, nRegularizerIters); //not used as beta=0 or nRegularizerIters=0

		printf("Data term : %f  Registration Err %f alpha %f\n", m_DataTermErr[iters], m_RegistrationErr[iters], alp);
	}


	// Stop the clock
	auto stop = std::chrono::high_resolution_clock::now();

	// Calculate the duration
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);

	// Print the duration
	std::cout << "Time taken by function: " << duration.count() << " seconds" << std::endl;


}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ForceDrivenRegistration::gradientDataTerm(Eigen::VectorXd& f, Eigen::VectorXd& dJdf, double& DataTerm)
{
	Pt3D closestPt;
	float barycentric[3];
	float distance, avgDist2, distX, distY, distZ;
	int faceID, cnt;
	size_t vertID;

	// solve displacement. The m_solver has been initialized by m_solver.compute(m_K). K u = f. Set f to solve u.
	m_u = m_solver.solve(f);

	for (size_t i = 0; i < (size_t)m_PreOpModel.nNodes; i++)
	{
		m_updateNodePositions[3 * i] = m_PreOpModel.nodesX[i] + (float)m_u[3 * i];
		m_updateNodePositions[3 * i + 1] = m_PreOpModel.nodesY[i] + (float)m_u[3 * i + 1];
		m_updateNodePositions[3 * i + 2] = m_PreOpModel.nodesZ[i] + (float)m_u[3 * i + 2];
	}

	m_grad.setZero(); // dJdu
	cnt = 0;
	avgDist2 = 0.0f;

	for (int i = 0; i < m_nIntraOpPoints; i++)
	{
		ClosestPointMP(m_IntraOpPoints[i], closestPt, barycentric, distance, faceID); // generate correspondence matrix per line on fly
		//printf("%d %f %f %f %d\n", i, barycentric[0], barycentric[1], barycentric[2], faceID);

		cnt++;
		avgDist2 += distance * distance;

		distX = closestPt.x - m_IntraOpPoints[i].x;
		distY = closestPt.y - m_IntraOpPoints[i].y;
		distZ = closestPt.z - m_IntraOpPoints[i].z;

		vertID = (size_t)m_MatchingSurfaceElements[faceID].nodeIDs[0];
		m_grad[3 * vertID] += distX * barycentric[0];
		m_grad[3 * vertID + 1] += distY * barycentric[0];
		m_grad[3 * vertID + 2] += distZ * barycentric[0];

		vertID = (size_t)m_MatchingSurfaceElements[faceID].nodeIDs[1];
		m_grad[3 * vertID] += distX * barycentric[1];
		m_grad[3 * vertID + 1] += distY * barycentric[1];
		m_grad[3 * vertID + 2] += distZ * barycentric[1];

		vertID = (size_t)m_MatchingSurfaceElements[faceID].nodeIDs[2];
		m_grad[3 * vertID] += distX * barycentric[2];
		m_grad[3 * vertID + 1] += distY * barycentric[2];
		m_grad[3 * vertID + 2] += distZ * barycentric[2];

		//update arrays for calculation of alpha
		m_baryCentric[3*i] = barycentric[0];
		m_baryCentric[3 * i + 1] = barycentric[1];
		m_baryCentric[3 * i + 2] = barycentric[2];

		m_nodeIDs[3 * i] = m_MatchingSurfaceElements[faceID].nodeIDs[0];
		m_nodeIDs[3 * i + 1] = m_MatchingSurfaceElements[faceID].nodeIDs[1];
		m_nodeIDs[3 * i + 2] = m_MatchingSurfaceElements[faceID].nodeIDs[2];

		m_dist[3 * i] = distX;
		m_dist[3 * i + 1] = distY;
		m_dist[3 * i + 2] = distZ;

	}

	DataTerm = avgDist2 / (double)cnt;
	//printf("Data term : %f\n", DataTerm);

	//dDataterm_dForce. 
	dJdf = m_solver.solve(m_grad); //dJdu = K dJdf (dJdu dudf = dJdf, dudf=K^-1)

	dJdf /= (double)cnt;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ForceDrivenRegistration::updateAlpha(double& alp)
{
	double top = 0.0f;
	double bot = 0.0f;

	m_fBar.setZero();

	for (int i = 0; i < m_nForceVertices; i++)
	{
		if (m_ZeroBoundaryVertices[i] == false)
		{
			size_t vertID = (size_t)m_ForceVertices[i];
			m_fBar[3 * vertID] = m_dJdf(3 * vertID);
			m_fBar[3 * vertID + 1] = m_dJdf[3 * vertID + 1];
			m_fBar[3 * vertID + 2] = m_dJdf[3 * vertID + 2];
		}
	}

	//use m_fBar as a temp array
	m_fBar = m_solver.solve(m_fBar);

	for (int i = 0; i < m_nIntraOpPoints; i++)
	{
		double fac1 = m_fBar[3*m_nodeIDs[3 * i]] * (double)m_baryCentric[3 * i] + m_fBar[3*m_nodeIDs[3 * i + 1]] * (double)m_baryCentric[3 * i + 1] +
			                   m_fBar[3*m_nodeIDs[3 * i + 2]] * (double)m_baryCentric[3 * i + 2];

		double fac2 = m_fBar[3*m_nodeIDs[3 * i] + 1] * (double)m_baryCentric[3 * i] + m_fBar[3*m_nodeIDs[3 * i + 1] + 1] * (double)m_baryCentric[3 * i + 1] +
			m_fBar[3*m_nodeIDs[3 * i + 2] + 1] * (double)m_baryCentric[3 * i + 2];

		double fac3 = m_fBar[3*m_nodeIDs[3 * i] + 2] * (double)m_baryCentric[3 * i] + m_fBar[3*m_nodeIDs[3 * i + 1] + 2] * (double)m_baryCentric[3 * i + 1] +
			m_fBar[3*m_nodeIDs[3 * i + 2] + 2] * (double)m_baryCentric[3 * i + 2];

		top += fac1 * m_dist[3 * i] + fac2 * m_dist[3 * i + 1] + fac3 * m_dist[3 * i + 2];
		bot += fac1 * fac1 + fac2 * fac2 + fac3 * fac3;

	}

	alp = top / bot;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void ForceDrivenRegistration::updateAlpha(double& alp)
//{
//	double top = 0.0f;
//	double bot = 0.0f;
//
//	m_fBar.setZero();
//
//	for (int i = 0; i < m_nForceVertices; i++)
//	{
//		if (m_ZeroBoundaryVertices[i] == false)
//		{
//			size_t vertID = (size_t)m_ForceVertices[i];
//			m_fBar[3 * vertID] = m_dJdf(3 * vertID);
//			m_fBar[3 * vertID + 1] = m_dJdf[3 * vertID + 1];
//			m_fBar[3 * vertID + 2] = m_dJdf[3 * vertID + 2];
//		}
//	}
//
//	//use m_fBar as a temp array
//	m_fBar = m_solver.solve(m_fBar);
//
//	for (int i = 0; i < m_nIntraOpPoints; i++)
//	{
//		double fac1 = m_fBar[3 * m_nodeIDs[3 * i]] * (double)m_baryCentric[3 * i] + m_fBar[3 * m_nodeIDs[3 * i + 1]] * (double)m_baryCentric[3 * i + 1] +
//			m_fBar[3 * m_nodeIDs[3 * i + 2]] * (double)m_baryCentric[3 * i + 2];
//
//		double fac2 = m_fBar[3 * m_nodeIDs[3 * i] + 1] * (double)m_baryCentric[3 * i] + m_fBar[3 * m_nodeIDs[3 * i + 1] + 1] * (double)m_baryCentric[3 * i + 1] +
//			m_fBar[3 * m_nodeIDs[3 * i + 2] + 1] * (double)m_baryCentric[3 * i + 2];
//
//		double fac3 = m_fBar[3 * m_nodeIDs[3 * i] + 2] * (double)m_baryCentric[3 * i] + m_fBar[3 * m_nodeIDs[3 * i + 1] + 2] * (double)m_baryCentric[3 * i + 1] +
//			m_fBar[3 * m_nodeIDs[3 * i + 2] + 2] * (double)m_baryCentric[3 * i + 2];
//
//		top += fac1 * m_dist[3 * i] + fac2 * m_dist[3 * i + 1] + fac3 * m_dist[3 * i + 2];
//		bot += fac1 * fac1 + fac2 * fac2 + fac3 * fac3;
//
//	}
//
//	if (m_bUseDisplacement)
//	{
//
//	}
//	alp = top / bot;
//
//
//
//}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ForceDrivenRegistration::Regularizer(float alp, int nIters)
{
	float* nodesX = m_PreOpModel.nodesX;
	float* nodesY = m_PreOpModel.nodesY;
	float* nodesZ = m_PreOpModel.nodesZ;

	float tau = (float)m_beta;

	m_fBar = m_f;

	int j;
	double dist;
	alp = 1.0f;

	//string fileName = "D:\\Source(2019)\\DIGS\\Liver\\test_regulaizer.txt";
	//FILE* fi;
	//fopen_s(&fi, fileName.c_str(), "w");

	for (size_t iters = 0; iters < nIters; iters++)
	{
		for (size_t i = 0; i < m_nForceVertices; i++)
		{
			if (m_ZeroBoundaryVertices[i] == false)
			{
				j = m_ForceVertices[i];
				double topX = 0.0;
				double topY = 0.0;
				double topZ = 0.0;
				double bot = 0.0;
				double delx = 0.0;
				double dely = 0.0;
				double delz = 0.0;
				//topX = m_fBar[3 * j];
				//topY = m_fBar[3 * j + 1];
				//topZ = m_fBar[3 * j + 2];
				bot = 1.0f;
				topX = 0.0f;
				topY = 0.0f;
				topZ = 0.0f;
				bot = 0.0f;
				for (size_t k = 0; k < m_nNeighbors[j]; k++)
				{
					size_t nodek = (size_t)m_Neighbors[j][k];
					delx = nodesX[j] - nodesX[nodek];
					dely = nodesY[j] - nodesY[nodek];
					delz = nodesZ[j] - nodesZ[nodek];
					dist = sqrt(delx * delx + dely * dely + delz * delz);
					topX += m_f[3 * nodek] / dist;
					topY += m_f[3 * nodek + 1] / dist;
					topZ += m_f[3 * nodek + 2] / dist;
					bot += 1.0f / dist;
				}
				topX = m_fBar[3 * j] + tau * alp * topX;
				topY = m_fBar[3 * j + 1] + tau * alp * topY;
				topZ = m_fBar[3 * j + 2] + tau * alp * topZ;
				bot = 1.0f + tau * alp * bot;
				m_f[3 * j] = topX / bot;
				m_f[3 * j + 1] = topY / bot;
				m_f[3 * j + 2] = topZ / bot;

				//fprintf_s(fi, "%f %f %f %f %f %f\n", m_fBar[3 * j], m_f[3 * j], m_fBar[3 * j + 1], m_f[3 * j + 1], m_fBar[3 * j + 2],m_f[3 * j + 2]);
			}
		}
	}

	//fclose(fi);
	//float* nodesX = m_PreOpModel.nodesX;
	//float* nodesY = m_PreOpModel.nodesY;
	//float* nodesZ = m_PreOpModel.nodesZ;

	//m_fBar = m_f;

	//int j;
	//float dist;

	////string fileName = "D:\\Source(2019)\\DIGS\\Liver\\test_regulaizer.txt";
	////FILE* fi;
	////fopen_s(&fi, fileName.c_str(), "w");

	//for (size_t iters = 0; iters < nIters; iters++)
	//{
	//	for (size_t i = 0; i < m_nForceVertices; i++)
	//	{
	//		if (m_ZeroBoundaryVertices[i] == false)
	//		{
	//			j = m_ForceVertices[i];
	//			float topX = 0.0f;
	//			float topY = 0.0f;
	//			float topZ = 0.0f;
	//			float bot = 0.0f;
	//			float delx = 0.0f;
	//			float dely = 0.0f;
	//			float delz = 0.0f;
	//			topX = m_beta * m_fBar[3 * j];
	//			topY = m_beta * m_fBar[3 * j + 1];
	//			topZ = m_beta * m_fBar[3 * j + 2];
	//			bot = m_beta;
	//			for (size_t k = 0; k < m_nNeighbors[j]; k++)
	//			{
	//				size_t nodek = (size_t)m_Neighbors[j][k];
	//				delx = nodesX[j] - nodesX[nodek];
	//				dely = nodesY[j] - nodesY[nodek];
	//				delz = nodesZ[j] - nodesZ[nodek];
	//				dist = sqrt(delx * delx + dely * dely + delz * delz);
	//				topX += m_f[3 * nodek] / dist;
	//				topY += m_f[3 * nodek + 1] / dist;
	//				topZ += m_f[3 * nodek + 2] / dist;
	//				bot += 1.0f / dist;
	//			}

	//			m_f[3 * j] = topX / bot;
	//			m_f[3 * j + 1] = topY / bot;
	//			m_f[3 * j + 2] = topZ / bot;

	//			//fprintf_s(fi, "%f %f %f %f %f %f\n", m_fBar[3 * j], m_f[3 * j], m_fBar[3 * j + 1], m_f[3 * j + 1], m_fBar[3 * j + 2],m_f[3 * j + 2]);
	//		}
	//	}
	//}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//find minimum distance to a face (triangluar element)
// https://www.gamedev.net/forums/topic/552906-closest-point-on-triangle/
// https://github.com/embree/embree/blob/master/tutorials/common/math/closest_point.h
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ForceDrivenRegistration::ClosestPoint(Pt3D point, Pt3D &closestPt, float* barycentric, float&minDistance, int &minFaceID)
{
	//bool ret = false;
	//float *nodesX = m_PreOpModel.nodesX;
	//float* nodesY = m_PreOpModel.nodesY;
	//float* nodesZ = m_PreOpModel.nodesZ;

	Pt3D v0, v1, v2, d0, e0, e1, testPt;
	float a, b, c, d, e, det, s, t, distance;
	size_t nodeID;
	minDistance = FLT_MAX;
	minFaceID = -1;

	s = 0.0f;
	t = 0.0f;

	for (int i = 0; i < m_nMatchingSurfaceElements; i++)
	{
		nodeID = (size_t)m_MatchingSurfaceElements[i].nodeIDs[0];
		v0.x = (float)m_updateNodePositions[3 * nodeID]; 
		v0.y = (float)m_updateNodePositions[3 * nodeID + 1];
		v0.z = (float)m_updateNodePositions[3 * nodeID + 2];

		nodeID = (size_t)m_MatchingSurfaceElements[i].nodeIDs[1];
		v1.x = (float)m_updateNodePositions[3 * nodeID];
		v1.y = (float)m_updateNodePositions[3 * nodeID + 1];
		v1.z = (float)m_updateNodePositions[3 * nodeID + 2];

		nodeID = (size_t)m_MatchingSurfaceElements[i].nodeIDs[2];
		v2.x = (float)m_updateNodePositions[3 * nodeID];
		v2.y = (float)m_updateNodePositions[3 * nodeID + 1];
		v2.z = (float)m_updateNodePositions[3 * nodeID + 2];

		//d0 = v0 - point
		d0.x = v0.x - point.x;  d0.y = v0.y - point.y;  d0.z = v0.z - point.z;

		//e0 = edge between v1 and v0
		e0.x = v1.x - v0.x;  e0.y = v1.y - v0.y;  e0.z = v1.z - v0.z;

		//e1 = edge between v2 and v0
		e1.x = v2.x - v0.x;  e1.y = v2.y - v0.y;  e1.z = v2.z - v0.z;

		//a = dot product between e0 and e0
		a = e0.x * e0.x + e0.y * e0.y + e0.z * e0.z;

		//b = dot product between e0 and e1
		b = e0.x * e1.x + e0.y * e1.y + e0.z * e1.z;

		//c = dot product between e1 and e1
		c = e1.x * e1.x + e1.y * e1.y + e1.z * e1.z;

		//d = dot product between e0 and d0
		d = e0.x * d0.x + e0.y * d0.y + e0.z * d0.z;

		//e = dot product between e1 and d0
		e = e1.x * d0.x + e1.y * d0.y + e1.z * d0.z;

		det = a * c - b * b;
		s = b * e - c * d;
		t = b * d - a * e;

		if (s + t < det)
		{
			if (s < 0.0f)
			{
				if (t < 0.0f)
				{
					if (d < 0.0f)
					{
						s = clamp(-d / a, 0.0f, 1.0f);
						t = 0.0f;
					}
					else
					{
						s = 0.0f;
						t = clamp(-e / c, 0.0f, 1.0f);
					}
				}
				else
				{
					s = 0.0f;
					t = clamp(-e / c, 0.0f, 1.0f);
				}
			}
			else if( t < 0.0f)
			{
				s = clamp(-d/a, 0.0f, 1.0f);
				t = 0.0f;
			}
			else
			{
				s /= det;
				t /= det;
			}
		}
		else
		{
			if (s < 0.0f)
			{
				float tmp0 = b + d;
				float tmp1 = c + e;
				if (tmp1 > tmp0)
				{
					float numer = tmp1 - tmp0;
					float denom = a - 2.0f * b + c;
					s = clamp(numer / denom, 0.0f, 1.0f);
					t = 1 - s;
				}
				else
				{
					t = clamp(-e / c, 0.0f, 1.0f);
					s = 0.0f;
				}
			}
			else if (t < 0.0f)
			{
				if (a + d > b + e)
				{
					float numer = c + e - b - d;
					float denom = a - 2.0f * b + c;
					s = clamp(numer / denom, 0.0f, 1.0f);
					t = 1.0f - s;
				}
				else
				{
					s = clamp(-e / c, 0.0f, 1.0f);
					t = 0.0f;
				}
			}
			else
			{
				float numer = c + e - b - d;
				float denom = a - 2.0f * b + c;
				s = clamp(numer / denom, 0.0f, 1.0f);
				t = 1.0f - s;
			}
		}

	
		testPt.x = v0.x + e0.x * s + e1.x * t;
		testPt.y = v0.y + e0.y * s + e1.y * t;
		testPt.z = v0.z + e0.z * s + e1.z * t;
		distance = (testPt.x - point.x) * (testPt.x - point.x) + (testPt.y - point.y) * (testPt.y - point.y) +
			(testPt.z - point.z) * (testPt.z - point.z);

		if (distance < minDistance)
		{
			minDistance = distance;
			minFaceID = i;
			closestPt.x = testPt.x;
			closestPt.y = testPt.y;
			closestPt.z = testPt.z;

			//barycentric.X * v0 + barycentric.Y * v1  + barycentric.Z * v2
			barycentric[0] = 1.0f - s - t;
			barycentric[1] = s;
			barycentric[2] = t;
		}
		
	}


	minDistance = sqrt(minDistance);

	return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Open MP Closest polint algorithm
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ForceDrivenRegistration::ClosestPointMP(Pt3D point, Pt3D& closestPt, float* barycentric, float& minDistance, int& minFaceID)
{

	minDistance = FLT_MAX;
	minFaceID = -1;

	Pt3D* testPtArr = new Pt3D[m_nMatchingSurfaceElements];
	float* sArr = new float[m_nMatchingSurfaceElements];
	float* tArr = new float[m_nMatchingSurfaceElements];
	float* distanceArr = new float[m_nMatchingSurfaceElements];

#pragma omp parallel for
	for (int i = 0; i < m_nMatchingSurfaceElements; i++)
	{

		ClosestPoint(point, testPtArr, sArr, tArr, distanceArr, i);
	}

	for (int i = 0; i < m_nMatchingSurfaceElements; i++)
	{
		if (distanceArr[i] < minDistance)
		{
			minFaceID = i;
			minDistance = distanceArr[i];
		}
	}

	barycentric[0] = 1.0f - sArr[minFaceID] - tArr[minFaceID];
	barycentric[1] = sArr[minFaceID];
	barycentric[2] = tArr[minFaceID];
	closestPt.x = testPtArr[minFaceID].x;
	closestPt.y = testPtArr[minFaceID].y;
	closestPt.z = testPtArr[minFaceID].z;

	minDistance = sqrt(minDistance);

	delete[] sArr;
	delete[] tArr;
	delete[] distanceArr;
	delete[] testPtArr;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//ClosestPt operator for use with ClosestPointMP version
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ForceDrivenRegistration::ClosestPoint(Pt3D point, Pt3D* testPtArr, float* sArr, float* tArr, float* distanceArr, int i)
{
	Pt3D v0, v1, v2, d0, e0, e1;
	float a, b, c, d, e, det;
	size_t nodeID;
	float s, t, distance;
	Pt3D testPt;

	nodeID = (size_t)m_MatchingSurfaceElements[i].nodeIDs[0];
	v0.x = m_updateNodePositions[3 * nodeID];
	v0.y = m_updateNodePositions[3 * nodeID + 1];
	v0.z = m_updateNodePositions[3 * nodeID + 2];

	nodeID = (size_t)m_MatchingSurfaceElements[i].nodeIDs[1];
	v1.x = m_updateNodePositions[3 * nodeID];
	v1.y = m_updateNodePositions[3 * nodeID + 1];
	v1.z = m_updateNodePositions[3 * nodeID + 2];

	nodeID = (size_t)m_MatchingSurfaceElements[i].nodeIDs[2];
	v2.x = m_updateNodePositions[3 * nodeID];
	v2.y = m_updateNodePositions[3 * nodeID + 1];
	v2.z = m_updateNodePositions[3 * nodeID + 2];

	//d0 = v0 - point
	d0.x = v0.x - point.x;  d0.y = v0.y - point.y;  d0.z = v0.z - point.z;

	//e0 = edge between v1 and v0
	e0.x = v1.x - v0.x;  e0.y = v1.y - v0.y;  e0.z = v1.z - v0.z;

	//e1 = edge between v2 and v0
	e1.x = v2.x - v0.x;  e1.y = v2.y - v0.y;  e1.z = v2.z - v0.z;

	//a = dot product between e0 and e0
	a = e0.x * e0.x + e0.y * e0.y + e0.z * e0.z;

	//b = dot product between e0 and e1
	b = e0.x * e1.x + e0.y * e1.y + e0.z * e1.z;

	//c = dot product between e1 and e1
	c = e1.x * e1.x + e1.y * e1.y + e1.z * e1.z;

	//d = dot product between e0 and d0
	d = e0.x * d0.x + e0.y * d0.y + e0.z * d0.z;

	//e = dot product between e1 and d0
	e = e1.x * d0.x + e1.y * d0.y + e1.z * d0.z;

	det = a * c - b * b;
	s = b * e - c * d;
	t = b * d - a * e;

	if (s + t < det)
	{
		if (s < 0.0f)
		{
			if (t < 0.0f)
			{
				if (d < 0.0f)
				{
					s = clamp(-d / a, 0.0f, 1.0f);
					t = 0.0f;
				}
				else
				{
					s = 0.0f;
					t = clamp(-e / c, 0.0f, 1.0f);
				}
			}
			else
			{
				s = 0.0f;
				t = clamp(-e / c, 0.0f, 1.0f);
			}
		}
		else if (t < 0.0f)
		{
			s = clamp(-d / a, 0.0f, 1.0f);
			t = 0.0f;
		}
		else
		{
			s /= det;
			t /= det;
		}
	}
	else
	{
		if (s < 0.0f)
		{
			float tmp0 = b + d;
			float tmp1 = c + e;
			if (tmp1 > tmp0)
			{
				float numer = tmp1 - tmp0;
				float denom = a - 2.0f * b + c;
				s = clamp(numer / denom, 0.0f, 1.0f);
				t = 1 - s;
			}
			else
			{
				t = clamp(-e / c, 0.0f, 1.0f);
				s = 0.0f;
			}
		}
		else if (t < 0.0f)
		{
			if (a + d > b + e)
			{
				float numer = c + e - b - d;
				float denom = a - 2.0f * b + c;
				s = clamp(numer / denom, 0.0f, 1.0f);
				t = 1.0f - s;
			}
			else
			{
				s = clamp(-e / c, 0.0f, 1.0f);
				t = 0.0f;
			}
		}
		else
		{
			float numer = c + e - b - d;
			float denom = a - 2.0f * b + c;
			s = clamp(numer / denom, 0.0f, 1.0f);
			t = 1.0f - s;
		}
	}


	testPt.x = v0.x + e0.x * s + e1.x * t;
	testPt.y = v0.y + e0.y * s + e1.y * t;
	testPt.z = v0.z + e0.z * s + e1.z * t;
	distance = (testPt.x - point.x) * (testPt.x - point.x) + (testPt.y - point.y) * (testPt.y - point.y) +
		(testPt.z - point.z) * (testPt.z - point.z);

	distanceArr[i] = distance;
	sArr[i] = s;
	tArr[i] = t;
	testPtArr[i].x = testPt.x;
	testPtArr[i].y = testPt.y;
	testPtArr[i].z = testPt.z;

}

void ForceDrivenRegistration::SaveUpdatedModel_test(string pathIn, string pathOut)
{
	

	m_u = m_solver.solve(m_f);
	//get undeformed liver surface and generate deformed surface model
	FILE* fi;

	string filename = pathOut + "deformed_liver_volumetric.off";
	//FILE* fi;
	fopen_s(&fi, filename.c_str(), "w");
	fprintf_s(fi, "OFF\n");
	fprintf_s(fi, "%d %d %d\n", m_PreOpModel.nNodes, m_PreOpModel.nElements, 0);

	for (int i = 0; i < m_PreOpModel.nNodes; i++)
	{
		fprintf_s(fi, "%f %f %f\n", m_updateNodePositions[3 * i], m_updateNodePositions[3 * i + 1], m_updateNodePositions[3 * i + 2]);
	}

	for (int i = 0; i < m_PreOpModel.nElements; i++)

		fprintf_s(fi, "%d %d %d %d %d\n", 4, m_PreOpModel.elements[i].nodeIDs[0], m_PreOpModel.elements[i].nodeIDs[1],
			m_PreOpModel.elements[i].nodeIDs[2], m_PreOpModel.elements[i].nodeIDs[3]);

	fclose(fi);

	
	string path = pathIn;

	
	string fileName = pathIn;
	fopen_s(&fi, fileName.c_str(), "r");

	fscanf_s(fi, "OFF\n");
	int nNodes, nElements, junk;
	fscanf_s(fi, "%d %d %d\n", &nNodes, &nElements, &junk);

	Element* SurfaceElements = new Element[nElements];
	int nSurfaceElements = nElements;

	int* SurfaceVertices = new int[nNodes];
	int nSurfaceVertices = nNodes;

	//convert vertices to match preop mesh node index
	float nodeX, nodeY, nodeZ;
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
		bool bFlag = false;
		for (int j = 0; j < m_PreOpModel.nNodes; j++)
		{
			if (nodeX == m_PreOpModel.nodesX[j] && nodeY == m_PreOpModel.nodesY[j] && nodeZ == m_PreOpModel.nodesZ[j])
			{
				SurfaceVertices[i] = j;
				bFlag = true;
				break;
			}

		}
		if (bFlag == false)
			printf("%d\n", i);
	}

	int vert0, vert1, vert2;

	for (int i = 0; i < nElements; i++)
	{
		fscanf_s(fi, "%d %d %d %d\n", &junk, &vert0, &vert1, &vert2);
		SurfaceElements[i].nodeIDs[0] = vert0;
		SurfaceElements[i].nodeIDs[1] = vert1;
		SurfaceElements[i].nodeIDs[2] = vert2;
	}

	fclose(fi);

	for (size_t i = 0; i < nNodes; i++)
	{
		int j = SurfaceVertices[i];
		m_updateNodePositions[3 * i] = m_PreOpModel.nodesX[j] + (float)m_u[3 * j];
		m_updateNodePositions[3 * i + 1] = m_PreOpModel.nodesY[j] + (float)m_u[3 * j + 1];
		m_updateNodePositions[3 * i + 2] = m_PreOpModel.nodesZ[j] + (float)m_u[3 * j + 2];
	}

	//location to save data 
	//path = "D:\\LiverMeshesDatabase\\Simulations\\KellysPhantomLiver\\Sim1\\Results\\26\\";
	path = pathOut;

	//save deformed liver surface
	filename = path + "deformed_liver_surface.off";
	//FILE* fi;
	fopen_s(&fi, filename.c_str(), "w");
	fprintf_s(fi, "OFF\n");
	fprintf_s(fi, "%d %d %d\n", nNodes, nElements, 0);

	for (int i = 0; i < nNodes; i++)
	{
		fprintf_s(fi, "%f %f %f\n", m_updateNodePositions[3 * i], m_updateNodePositions[3 * i + 1], m_updateNodePositions[3 * i + 2]);
	}

	for (int i = 0; i < nElements; i++)

		fprintf_s(fi, "%d %d %d %d\n", 3, SurfaceElements[i].nodeIDs[0], SurfaceElements[i].nodeIDs[1],
			SurfaceElements[i].nodeIDs[2]);

	fclose(fi);
	printf("Finished!");
}


