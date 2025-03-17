#pragma once
#include <string>
#include "..\eigen-eigen-323c052e1731/Eigen/Dense"
#include "..\eigen-eigen-323c052e1731/Eigen/Sparse"
#include "Common.h"

using namespace std;
class FEM3DD
{

public:

protected:

	double                       m_poissonRatio;

	int                         m_forceCount;
	int                         m_nDisplacementConstraints;
	int                         m_nTractionConstraints;
	int                         m_elementCount;
	int                      	m_nodesCount;
	int                         m_boundaryElementCount;
	int                         m_nDims;
	bool                        m_bTetrahedron;
	bool                        m_bPiecewiseModulus;
	bool                        m_bLinearElasticMaterial;

	float* m_nodesX;
	float* m_nodesY;
	float* m_nodesZ;

	double m_modulus;

	Element* m_elements;
	Element* m_boundaryElements;
	Constraint* m_displacementConstraints;
	Constraint* m_tractionConstraints; 

	//Eigen::VectorXd m_u;
	//Eigen::VectorXd m_f;
	//Eigen::SparseMatrix<double>  m_K;

	std::vector<Eigen::Triplet<double> > m_triplets;

	void CalculateElementMatrix_Tetrahedron(Eigen::VectorXd &u, int whichElement);
	
	void LinearElastic(Eigen::MatrixXd& C_SE, double modulus);

public:
	FEM3DD();
	~FEM3DD();

	void setup(float* nodesX, float* nodesY, float* nodesZ, Element* elements, int nNodes, int nElements,
		float modulus, float poissonRatio, bool bLinearElasticMaterial);

	void setup(Mesh &mesh, float modulus, float poissonRatio, bool bLinearElasticMaterial);

	void calculateStiffnessMatrix(Eigen::VectorXd& u, Eigen::SparseMatrix<double>& stiffnessMatrix);

	int GetNumberElements() { return m_elementCount; };
	int GetNumberNodes() { return m_nodesCount; };

	float* GetNodesX() { return m_nodesX; };
	float* GetNodesY() { return m_nodesY; };
	float* GetNodesZ() { return m_nodesZ; };
	Element* GetElements() { return m_elements; };


};


void compute_BN_matrix(Eigen::MatrixXd& B, Eigen::MatrixXd& F, Eigen::MatrixXd& BN);