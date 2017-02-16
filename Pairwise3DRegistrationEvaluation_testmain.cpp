#include "Pairwise3DRegistrationEvaluation.h"

#include <stdexcept>


using namespace std;



/** \brief Demonstrative wrapper of the algorithm to test. Implement PairwiseRegistration() method.
*/
class TestPairwise3DRegistrationAlgorithm : public Pairwise3DRegistrationEvaluation::IPairwise3DRegistrationAlgorithm
{
	public:

	vtkTransform* PairwiseRegistration(vtkPolyData* cloud_trg, vtkPolyData* cloud_ref, const int idx_trg, const int idx_ref, double &secCPUtime)
	{
		//access cloud_trg so to fill your data
		float* ptrPoints_trg = Pairwise3DRegistrationEvaluation::GetPolyDataPointsPointer(cloud_trg);	//x0,y0,z0,x1,y1,z1,x2,y2,z2. etc...
		vtkIdType nPoints_trg = Pairwise3DRegistrationEvaluation::GetPolyDataNumberOfPoints(cloud_trg);

		float* ptrNormals_trg = Pairwise3DRegistrationEvaluation::GetPolyDataPointNormalsPointer(cloud_trg);	//nx0,ny0,nz0,nx1,ny1,nz1,nx2,ny2,nz2. etc...

		vtkIdType* ptrTrianglesIds_trg = Pairwise3DRegistrationEvaluation::GetPolyDataTrianglesPointer(cloud_trg); //3,t00,t01,t02,3,t10,t11,t12,3,t20,t21,t22...etc each triangle is a quadruple!!!
		vtkIdType nTriangles_trg = cloud_trg->GetNumberOfPolys();

		//access cloud_ref so to fill your data
		float* ptrPoints_ref = Pairwise3DRegistrationEvaluation::GetPolyDataPointsPointer(cloud_ref);	//x0,y0,z0,x1,y1,z1,x2,y2,z2. etc...
		vtkIdType nPoints_ref = Pairwise3DRegistrationEvaluation::GetPolyDataNumberOfPoints(cloud_ref);

		float* ptrNormals_ref = Pairwise3DRegistrationEvaluation::GetPolyDataPointNormalsPointer(cloud_ref);	//nx0,ny0,nz0,nx1,ny1,nz1,nx2,ny2,nz2. etc...

		vtkIdType* ptrTrianglesIds_ref = Pairwise3DRegistrationEvaluation::GetPolyDataTrianglesPointer(cloud_ref); //3,t00,t01,t02,3,t10,t11,t12,3,t20,t21,t22...etc each triangle is a quadruple!!!
		vtkIdType nTriangles_ref = cloud_ref->GetNumberOfPolys();


		Pairwise3DRegistrationEvaluation::CpuTimeProfiler cpuTime;
		//perform registration here!!!
		secCPUtime = cpuTime.GetElapsedSecs();

		//if(registration failed)
		//	return NULL;

		//fill the final rigid motion
		vtkTransform* rigidMotion = vtkTransform::New();
		rigidMotion->GetMatrix()->Element[0][0] = 1.0;
		rigidMotion->GetMatrix()->Element[0][1] = 0.0;
		rigidMotion->GetMatrix()->Element[0][2] = 0.0;
		rigidMotion->GetMatrix()->Element[0][3] = 0.0;
		rigidMotion->GetMatrix()->Element[1][0] = 0.0;
		rigidMotion->GetMatrix()->Element[1][1] = 1.0;
		rigidMotion->GetMatrix()->Element[1][2] = 0.0;
		rigidMotion->GetMatrix()->Element[1][3] = 0.0;
		rigidMotion->GetMatrix()->Element[2][0] = 0.0;
		rigidMotion->GetMatrix()->Element[2][1] = 0.0;
		rigidMotion->GetMatrix()->Element[2][2] = 1.0;
		rigidMotion->GetMatrix()->Element[2][3] = 0.0;
		rigidMotion->GetMatrix()->Element[3][0] = 0.0;
		rigidMotion->GetMatrix()->Element[3][1] = 0.0;
		rigidMotion->GetMatrix()->Element[3][2] = 0.0;
		rigidMotion->GetMatrix()->Element[3][3] = 1.0;

		return rigidMotion;	
	};
};


int main(int argc, char** argv)
{
	string absDatasetPath = argv[1];	//folder containing the partial views to register
	string absResultFilename = argv[2];	//output csv file reporting the result of the registration

	//Instantiate the wrapper of the algorithm to test
	Pairwise3DRegistrationEvaluation::IPairwise3DRegistrationAlgorithm* ptrTestAlgorithm = new TestPairwise3DRegistrationAlgorithm();

	//Instantiate the benchmark
	Pairwise3DRegistrationEvaluation::Pairwise3DRegistrationBenchmark benchmark; 

	//Set the parameters
	benchmark.SetDatasetPath(absDatasetPath);
	benchmark.SetAlgorithm(ptrTestAlgorithm);
	//optional params
	//vector<pair<int, int>> vViewPairs;
	//vViewPairs.push_back(pair<int,int>(1,2));
	//vViewPairs.push_back(pair<int,int>(1,3));
	//vViewPairs.push_back(pair<int,int>(2,3));
	//benchmark.SetViewPairs(vViewPairs);
	//benchmark.SetMeshFileExtension("ply");
	//benchmark.m_lowMemory = true;
	//benchmark.m_RMSEthresh = 5.0;
	//benchmark.m_skipICP = false;
	//benchmark.m_ICPnMaxIters = 10;
	//benchmark.m_ICPeps = 0.1;
	//benchmark.m_ICPmaxRadius = 8.0;
	//benchmark.m_ICPsamplingPerc = 1.0;

	//Load partial views and ground truth
	benchmark.PrepareEvaluation();

	//Run the benchmark on all the view pairs
	benchmark.Evaluate();

	benchmark.PrintResults(absResultFilename, benchmark.GetNRegistrations(), benchmark.GetRMSE(), benchmark.GetTotalCPUtime() );

	//benchmark.ComputeOverlappingAreas(absResultFilename);	//Use this function for computing the overlapping area of all the view pairs after the ground truth alignment

	delete ptrTestAlgorithm;
}

