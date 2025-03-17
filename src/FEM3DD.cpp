#include "..\eigen-eigen-323c052e1731/Eigen/Dense"
#include "..\eigen-eigen-323c052e1731/Eigen/Sparse"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <time.h>
#include "Common.h"
//#include "tools.h"
//#include "ImageIO.h"
#include "FEM3DD.h"

//void compute_BN_matrix(Eigen::MatrixXd& B, Eigen::MatrixXd& F, Eigen::MatrixXd& BN);

////////////////////////////////////////////////////////////////////////////////////////
//Modulus is defined at vertices
////////////////////////////////////////////////////////////////////////////////////////
FEM3DD::FEM3DD()
{
	m_nodesX = NULL;
	m_nodesY = NULL;
	m_modulus = NULL;
	m_elements = NULL;
	m_boundaryElements = NULL;
	m_displacementConstraints = NULL;
	m_tractionConstraints = NULL;
	
	m_poissonRatio = 0.49f;

	m_nDims = 3;
}

////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////
FEM3DD::~FEM3DD()
{

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FEM3DD::setup(float* nodesX, float* nodesY, float* nodesZ, Element* elements, int nNodes, int nElements,
	float modulus, float poissonRatio, bool bLinearElasticMaterial)
{
	m_nodesX = nodesX;
	m_nodesY = nodesY;
	m_nodesZ = nodesZ;

	m_modulus = (double)modulus;
	m_poissonRatio = (double)poissonRatio;
	m_nodesCount = (double)nNodes;

	m_elements = elements;
	m_elementCount = nElements;

	m_bTetrahedron = true;
	m_bPiecewiseModulus = true;
	m_bLinearElasticMaterial = bLinearElasticMaterial;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FEM3DD::setup(Mesh &mesh, float modulus, float poissonRatio, bool bLinearElasticMaterial)
{
	m_nodesX = mesh.nodesX;
	m_nodesY = mesh.nodesY;
	m_nodesZ = mesh.nodesZ;

	m_modulus = (double)modulus;
	m_poissonRatio = (double)poissonRatio;
	m_nodesCount = (double)mesh.nNodes;

	m_elements = mesh.elements;
	m_elementCount = mesh.nElements;

	m_bTetrahedron = true;
	m_bPiecewiseModulus = true;
	m_bLinearElasticMaterial = bLinearElasticMaterial;

	/*m_u.resize(3 * m_nodesCount);
	m_u.setZero();

	m_f.resize(3 * m_nodesCount);
	m_f.setZero();

	m_K.resize(3 * m_nodesCount, 3 * m_nodesCount);
	m_K.setZero();

	calculateStiffnessMatrix();

	SetInitialConditions();

	Eigen::SparseLU<Eigen::SparseMatrix<double> > solver0;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;

	solver.compute(m_K);

	for (int i = 0; i < 3 * m_nodesCount; i++)
		m_u(i) = 1.0;

	m_f = m_K * m_u;*/
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//  
//////////////////////////////////////////////////////////////////////////////////////////////////////
void FEM3DD::calculateStiffnessMatrix(Eigen::VectorXd& u, Eigen::SparseMatrix<double>& stiffnessMatrix)
{
	//std::vector<Eigen::Triplet<float> > triplets;

	m_triplets.clear();

	for (int i = 0; i < m_elementCount; i++)
	{
		CalculateElementMatrix_Tetrahedron(u, i);
		//printf("\n");
		//CalculateElementMatrix_Tetra0(u, i);
	}
	

	stiffnessMatrix.resize(3 * m_nodesCount, 3 * m_nodesCount);
	stiffnessMatrix.setZero();

	stiffnessMatrix.setFromTriplets(m_triplets.begin(), m_triplets.end());

	m_triplets.clear();
}



//////////////////////////////////////////////////////////////////////////////////////////////////////
////K(12,12) : TangentStiffnessMatrix
//////////////////////////////////////////////////////////////////////////////////////////////////////
void FEM3DD::CalculateElementMatrix_Tetrahedron(Eigen::VectorXd& displacement, int whichElement)
{
	int nDims = 3;

	Eigen::Matrix<double, 12, 12> K;

	K.setZero();

	int node[4];
	node[0] = m_elements[whichElement].nodeIDs[0];
	node[1] = m_elements[whichElement].nodeIDs[1];
	node[2] = m_elements[whichElement].nodeIDs[2];
	node[3] = m_elements[whichElement].nodeIDs[3];

	Eigen::MatrixXd x(3, 4);
	x << (double)m_nodesX[node[0]], (double)m_nodesX[node[1]], (double)m_nodesX[node[2]], (double)m_nodesX[node[3]],
		(double)m_nodesY[node[0]], (double)m_nodesY[node[1]], (double)m_nodesY[node[2]], (double)m_nodesY[node[3]],
		(double)m_nodesZ[node[0]], (double)m_nodesZ[node[1]], (double)m_nodesZ[node[2]], (double)m_nodesZ[node[3]];

	double X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4;
	X1 = (double)m_nodesX[node[0]]; X2 = (double)m_nodesX[node[1]];  X3 = (double)m_nodesX[node[2]]; X4 = (double)m_nodesX[node[3]];
	Y1 = (double)m_nodesY[node[0]]; Y2 = (double)m_nodesY[node[1]];  Y3 = (double)m_nodesY[node[2]]; Y4 = (double)m_nodesY[node[3]];
	Z1 = (double)m_nodesZ[node[0]]; Z2 = (double)m_nodesZ[node[1]];  Z3 = (double)m_nodesZ[node[2]]; Z4 = (double)m_nodesZ[node[3]];

	//Eigen::MatrixXd matS(4, 4);
	//matS.setZero();
	//matS << 1., X1, Y1, Z1,
	//	1., X2, Y2, Z2,
	//	1., X2, Y2, Z2,
	//	1., X3, Y3, Z3;

	//float vol0 = matS.determinant() / 6.;

	double detJ = -X1 * Y2 * Z3 + X1 * Y2 * Z4 + X1 * Y3 * Z2 - X1 * Y3 * Z4 - X1 * Y4 * Z2 + X1 * Y4 * Z3
		+ X2 * Y1 * Z3 - X2 * Y1 * Z4 - X2 * Y3 * Z1 + X2 * Y3 * Z4 + X2 * Y4 * Z1 - X2 * Y4 * Z3
		- X3 * Y1 * Z2 + X3 * Y1 * Z4 + X3 * Y2 * Z1 - X3 * Y2 * Z4 - X3 * Y4 * Z1 + X3 * Y4 * Z2
		+ X4 * Y1 * Z2 - X4 * Y1 * Z3 - X4 * Y2 * Z1 + X4 * Y2 * Z3 + X4 * Y3 * Z1 - X4 * Y3 * Z2;

	double vol = detJ / 6.0;

	if (vol < 0.0)
				printf("negative jacobian for element %d\n", whichElement);

	//printf("element vol : %f\n", vol);

	Eigen::MatrixXd dN_dX(4,3);

	dN_dX <<   -Y2 * Z3 + Y2 * Z4 + Y3 * Z2 - Y3 * Z4 - Y4 * Z2 + Y4 * Z3,
			    X2 * Z3 - X2 * Z4 - X3 * Z2 + X3 * Z4 + X4 * Z2 - X4 * Z3,
			   -X2 * Y3 + X2 * Y4 + X3 * Y2 - X3 * Y4 - X4 * Y2 + X4 * Y3,
			    Y1 * Z3 - Y1 * Z4 - Y3 * Z1 + Y3 * Z4 + Y4 * Z1 - Y4 * Z3,
			   -X1 * Z3 + X1 * Z4 + X3 * Z1 - X3 * Z4 - X4 * Z1 + X4 * Z3,
			    X1 * Y3 - X1 * Y4 - X3 * Y1 + X3 * Y4 + X4 * Y1 - X4 * Y3,
			   -Y1 * Z2 + Y1 * Z4 + Y2 * Z1 - Y2 * Z4 - Y4 * Z1 + Y4 * Z2,
			    X1 * Z2 - X1 * Z4 - X2 * Z1 + X2 * Z4 + X4 * Z1 - X4 * Z2,
			   -X1 * Y2 + X1 * Y4 + X2 * Y1 - X2 * Y4 - X4 * Y1 + X4 * Y2,
			    Y1 * Z2 - Y1 * Z3 - Y2 * Z1 + Y2 * Z3 + Y3 * Z1 - Y3 * Z2,
			   -X1 * Z2 + X1 * Z3 + X2 * Z1 - X2 * Z3 - X3 * Z1 + X3 * Z2,
			    X1 * Y2 - X1 * Y3 - X2 * Y1 + X2 * Y3 + X3 * Y1 - X3 * Y2;

	dN_dX /= detJ;

	//deformation gradient ( 3 x 3 )
	Eigen::MatrixXd F(3, 3);
	F = Eigen::Matrix3d::Identity();

	//create nonlinear displacement-strain matrix BN (6x12)
	Eigen::MatrixXd BN(6, 12);
	compute_BN_matrix(dN_dX, F, BN);

	Eigen::MatrixXd C_SE(6, 6);
	LinearElastic(C_SE, m_modulus);

	//float vol = detJ / 6.0f;
	K = BN.transpose() * C_SE * BN * vol;

	//triplets for assembling global stiffness matrix K
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{

			Eigen::Triplet<double> trplt11(3 * node[i] + 0, 3 * node[j] + 0, K(3 * i + 0, 3 * j + 0));
			Eigen::Triplet<double> trplt12(3 * node[i] + 0, 3 * node[j] + 1, K(3 * i + 0, 3 * j + 1));
			Eigen::Triplet<double> trplt13(3 * node[i] + 0, 3 * node[j] + 2, K(3 * i + 0, 3 * j + 2));

			Eigen::Triplet<double> trplt21(3 * node[i] + 1, 3 * node[j] + 0, K(3 * i + 1, 3 * j + 0));
			Eigen::Triplet<double> trplt22(3 * node[i] + 1, 3 * node[j] + 1, K(3 * i + 1, 3 * j + 1));
			Eigen::Triplet<double> trplt23(3 * node[i] + 1, 3 * node[j] + 2, K(3 * i + 1, 3 * j + 2));

			Eigen::Triplet<double> trplt31(3 * node[i] + 2, 3 * node[j] + 0, K(3 * i + 2, 3 * j + 0));
			Eigen::Triplet<double> trplt32(3 * node[i] + 2, 3 * node[j] + 1, K(3 * i + 2, 3 * j + 1));
			Eigen::Triplet<double> trplt33(3 * node[i] + 2, 3 * node[j] + 2, K(3 * i + 2, 3 * j + 2));

			m_triplets.push_back(trplt11);
			m_triplets.push_back(trplt12);
			m_triplets.push_back(trplt13);
			m_triplets.push_back(trplt21);
			m_triplets.push_back(trplt22);
			m_triplets.push_back(trplt23);
			m_triplets.push_back(trplt31);
			m_triplets.push_back(trplt32);
			m_triplets.push_back(trplt33);

		}
	}


}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FEM3DD::LinearElastic(Eigen::MatrixXd& C_SE, double modulus)
{
	double lam = m_poissonRatio * modulus / ((1.0 + m_poissonRatio) * (1.0 - 2.0 * m_poissonRatio));
	double mu = modulus / (2.0f * (1.0f + m_poissonRatio));

	C_SE << lam + 2.0f * mu, lam, lam, 0.0f, 0.0f, 0.0f,
		lam, lam + 2.0f * mu, lam, 0.0f, 0.0f, 0.0f,
		lam, lam, lam + 2.0f * mu, 0.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 0.0f, mu, 0.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 0.0f, mu, 0.0f,
		0.0f, 0.0f, 0.0f, 0.0f, 0.0f, mu;

}

////////////////////////////////////////////////////////////////////////////////////////////////////
//create nonlinear displacement-strain matrix BN 
////////////////////////////////////////////////////////////////////////////////////////////////////
void compute_BN_matrix(Eigen::MatrixXd& B, Eigen::MatrixXd& F, Eigen::MatrixXd& BN)
{
	int nNodes = (int)B.rows();
	int nDims = (int)B.cols();

	double F11 = F(0, 0);
	double F12 = F(0, 1);
	double F21 = F(1, 0);
	double F22 = F(1, 1);

	if (nDims == 2)
	{
		for (int i = 0; i < nNodes; i++)
		{
			//ETAXX
			BN(0, 2 * i) = F11 * B(i, 0);
			BN(0, 2 * i + 1) = F21 * B(i, 0);

			//ETAYY
			BN(1, 2 * i) = F12 * B(i, 1);
			BN(1, 2 * i + 1) = F22 * B(i, 1);

			//ETAXY
			BN(2, 2 * i) = F11 * B(i, 1) + F12 * B(i, 0);
			BN(2, 2 * i + 1) = F21 * B(i, 1) + F22 * B(i, 0);
		}
	}
	else
	{
		double F13 = F(0, 2);
		double F31 = F(2, 0);
		double F23 = F(1, 2);
		double F32 = F(2, 1);
		double F33 = F(2, 2);

		for (int i = 0; i < nNodes; i++)
		{
			//ETAXX
			BN(0, 3 * i) = F11 * B(i, 0);
			BN(0, 3 * i + 1) = F21 * B(i, 0);
			BN(0, 3 * i + 2) = F31 * B(i, 0);

			//ETAYY
			BN(1, 3 * i) = F12 * B(i, 1);
			BN(1, 3 * i + 1) = F22 * B(i, 1);
			BN(1, 3 * i + 2) = F32 * B(i, 1);

			//ETAZZ
			BN(2, 3 * i) = F13 * B(i, 2);
			BN(2, 3 * i + 1) = F23 * B(i, 2);
			BN(2, 3 * i + 2) = F33 * B(i, 2);

			//ETAXY
			BN(3, 3 * i) = F11 * B(i, 1) + F12 * B(i, 0);
			BN(3, 3 * i + 1) = F21 * B(i, 1) + F22 * B(i, 0);
			BN(3, 3 * i + 2) = F31 * B(i, 1) + F32 * B(i, 0);

			//ETAYZ
			BN(4, 3 * i) = F12 * B(i, 2) + F13 * B(i, 1);
			BN(4, 3 * i + 1) = F22 * B(i, 2) + F23 * B(i, 1);
			BN(4, 3 * i + 2) = F32 * B(i, 2) + F33 * B(i, 1);

			//ETAXZ
			BN(5, 3 * i) = F11 * B(i, 2) + F13 * B(i, 0);
			BN(5, 3 * i + 1) = F21 * B(i, 2) + F23 * B(i, 0);
			BN(5, 3 * i + 2) = F31 * B(i, 2) + F33 * B(i, 0);
		}
	}
}

//void TestFEMDD()
//{
//	FILE* fi, * fi1;
//
//	string path = "D:\\LiverMeshesDatabase\\Simulations\\KellysPhantomLiver\\Sim1\\";
//
//	string fileName = path + "undeformed_liver_volumetric.dat";
//	fopen_s(&fi, fileName.c_str(), "r");
//
//	int nNodes, nElements, nSurfaceNodes;
//
//	fscanf_s(fi, "%d %d\n", &nNodes, &nSurfaceNodes);
//
//	float* nodesX = new float[nNodes];
//	float* nodesY = new float[nNodes];
//	float* nodesZ = new float[nNodes];
//
//	float node1, node2, node3;
//	for (int i = 0; i < nNodes; i++)
//	{
//		fscanf_s(fi, "%f %f %f\n", &node1, &node2, &node3);
//		node1 = floor(node1 * 100.0f);
//		node1 /= 100.0f;
//		node2 = floor(node2 * 100.0f);
//		node2 /= 100.0f;
//		node3 = floor(node3 * 100.0f);
//		node3 /= 100.0f;
//		nodesX[i] = node1;
//		nodesY[i] = node2;
//		nodesZ[i] = node3;
//	}
//
//	fscanf_s(fi, "%d\n", &nElements);
//
//	Element* elements = new Element[nElements];
//
//	int v1, v2, v3, v4;
//	for (int i = 0; i < nElements; i++)
//	{
//		fscanf_s(fi, "%d %d %d %d\n", &v1, &v2, &v3, &v4);
//		//change oreder of nodes for the FEM to have positive volume
//		elements[i].nodeIDs[0] = v1 - 1;
//		elements[i].nodeIDs[2] = v2 - 1;
//		elements[i].nodeIDs[1] = v3 - 1;
//		elements[i].nodeIDs[3] = v4 - 1;
//	}
//
//	fclose(fi);
//
//	Mesh mesh;
//	mesh.nElements = nElements;
//	mesh.nNodes = nNodes;
//	mesh.nodesX = nodesX;
//	mesh.nodesY = nodesY;
//	mesh.nodesZ = nodesZ;
//	mesh.elements = elements;
//
//	FEM3DD myFEM;
//
//	float modulus = 1.0f;
//	float poissonRatio = 0.49f;
//	bool bLinearElasticMaterial = true;
//	myFEM.setup(mesh, modulus, poissonRatio, bLinearElasticMaterial);
//
//}