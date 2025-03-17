#pragma once

#include "Common.h"
#include "FEM3DD.h"
#include "ImageIO.h"
#include "Eigen/Sparse"
#include <omp.h>
#include <math.h>
#include "float.h"
#include "Matrix3D.h"
#include "pts_tools.h"
#include <cstdlib>
#include <ctime>
#include "downsamplePointCloud.h"
#include <chrono>

class ForceDrivenRegistration
{
public:

	static struct InputData
	{
		Mesh PreOpModel;
		int* ForceVertices;
		Element* ForceElements;
		Element *MatchingSurfaceElements;
		int* MatchingSurfaceVertices;
		Pt3D* IntraOpPoints;
		Pt3D* PreOpFiducials;
		Pt3D* IntraOpFiducials;
		Pt3D* DeformedVolumeGTVertices;
		bool* ZeroBoundaryVertices;
		int nForceVertices;
		int nForceElements;
		int nMatchingSurfaceElements;
		int nMatchingSurfaceVertices;
		int nIntraOpPoints;
		int nFiducialPoints;
		int nZeroBoundaryVertices;
		int nDeformedVolumeGTVertices;
		double alpha;
		double beta;
		double gamma;
		double tau;
		double modulus; // 1.0;
		double poissonRatio;//0.49;
		bool bUseNesterov;
	};

private:

	Mesh m_PreOpModel;

	int* m_ForceVertices;
	int m_nForceVertices;

	Element* m_ForceElements;
	int m_nForceElements;

	bool* m_ZeroBoundaryVertices;
	int m_nZeroBoundaryVertices;

	Element* m_MatchingSurfaceElements;
	int m_nMatchingSurfaceElements;

	int* m_MatchingSurfaceVertices;
	int m_nMatchingSurfaceVertices;

	Pt3D* m_IntraOpPoints;
	int m_nIntraOpPoints;

	Pt3D* m_DeformedVolumeGTVertices;
	int m_nDeformedVolumeGTVertices;

	double m_alpha, m_beta, m_gamma;

	FEM3DD m_FEM;

	double m_modulus;
	double m_poissonRatio;
	bool m_bLinearElasticMaterial;

	Eigen::SparseMatrix<double>  m_K;

	Eigen::VectorXd m_u;
	Eigen::VectorXd m_f;
	Eigen::VectorXd m_fBar;
	Eigen::VectorXd m_grad;
	Eigen::VectorXd m_dJdf;
	Eigen::VectorXf m_updateNodePositions;

	//arays for updating alpha
	double* m_dist;
	float* m_baryCentric;
	int* m_nodeIDs;

	int m_nFemNodes;

	int** m_Neighbors;
	int* m_nNeighbors;

	double m_lambda;

	double m_sigma2;

	double m_alp;

	double m_tau;

	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > m_solver;

	Pt3D* m_deformedPts;
	Pt3D* m_undeformedPts;

	double* m_DataTermErr;
	double* m_RegistrationErr;
	int m_nAlgorithmIters;

	bool m_bUseDisplacement;
	bool m_bUseNesterov;

	void Regularizer(float alp, int nIters);
	void Adjacency_Matrix();
	void SetZeroBoundaryConditions();
	void SetInitialConditions();
	void ClosestPoint(Pt3D point, Pt3D& closestPt, float* barycentric, float& minDistance, int& faceID);
	void ClosestPointMP(Pt3D point, Pt3D& closestPt, float* barycentric, float& minDistance, int& faceID);
	void ClosestPoint(Pt3D point, Pt3D* testPtArr, float* sArr, float* tArr, float* distanceArr, int i);
	void gradientDataTerm(Eigen::VectorXd &f, Eigen::VectorXd& dJdf, double&DataTerm);
	void updateAlpha(double& alp);

public:
	ForceDrivenRegistration();
	~ForceDrivenRegistration();

	void setup(InputData data);
	void SaveUpdatedModel_dev(string pathIn, string pathOut); // Also save estimated forces.
	void solve(int nIters, int nRegularizerIters, int nDataIters);
	void GroundTruthError(double& avg, double& sd, double& max, bool bRMS);

	void SaveUpdatedModel_test(string pathIn, string pathOut);
};
