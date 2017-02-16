#include "Pairwise3DRegistrationEvaluation_FineRegistration.h"

#include "vtkMath.h"
#include <stdexcept>



#ifndef GICP_IS_NOT_INCLUDED
	#include "gicp.h"
	using namespace dgc::gicp;
#endif
	
using namespace std;

Pairwise3DRegistrationEvaluation::GICPPairwiseRegistration::GICPPairwiseRegistration(int numMaxIter, double epsilon, float maxSearchRadius, bool verbose)
:PairwiseRegistration(false, false, verbose),
m_numMaxIter(numMaxIter),
m_epsilon(epsilon),
m_maxSearchRadius(maxSearchRadius)
{
}


vtkTransform* Pairwise3DRegistrationEvaluation::GICPPairwiseRegistration::Register(vtkPolyData* polyData1In, vtkPolyData* polyData2In)
{
#ifdef GICP_IS_NOT_INCLUDED
	throw runtime_error("Generalized ICP not included. Include GICP (http://www.robots.ox.ac.uk/~avsegal/generalized_icp.html) or define GICP_IS_NOT_INCLUDED and use ICPPairwiseRegistration");
#else

	if(m_transform)
	{
		m_transform->Delete();
	}
	m_transform = vtkTransform::New();
	m_transform->PostMultiply();


	vtkPolyData* transformedPoly = Duplicate(polyData2In);


	GICPPointSet p1, p2;
	dgc_transform_t t_base, t0, t1;

	// set up the transformations
	dgc_transform_identity(t_base);
	dgc_transform_identity(t0);
	dgc_transform_identity(t1);


	//p2 get polyData1In and vice versa tohave the rototranslation that move p2 to p1
	p2.SetPoints(GetPolyDataPointsPointer(polyData1In), polyData1In->GetNumberOfPoints());
	p1.SetPoints(GetPolyDataPointsPointer(polyData2In), polyData2In->GetNumberOfPoints());

	// build kdtrees and normal matrices
	if(m_verbose)
	{
		cout << "Building KDTree and computing surface normals/matrices..." << endl;
	}

	//double gicp_epsilon = 1e-3;
	//p1.SetGICPEpsilon(gicp_epsilon);
	//p2.SetGICPEpsilon(gicp_epsilon);
	p1.SetEpsilon(m_epsilon);
	p2.SetEpsilon(m_epsilon);  
	p1.BuildKDTree();
	p1.ComputeMatrices();
	p2.BuildKDTree();
	p2.ComputeMatrices();

	// align the point clouds
	if(m_verbose)
	{
		cout << "Aligning point cloud..." << endl;
	}
	dgc_transform_copy(t1, t0);
	//p2.SetMaxIterationInner(8);
	p2.SetMaxIteration(m_numMaxIter);
	int iterations = p2.AlignScan(&p1, t_base, t1, m_maxSearchRadius);

	// print the result
	if(m_verbose)
	{
		cout << "Converged. Num Itarations: " << iterations << endl;
		dgc_transform_print(t1, "Rototranslation Matrix");
	}

	transformedPoly->Delete();
	
	m_transform->SetMatrix( (double*)t1 );
#endif

	return m_transform;
}



Pairwise3DRegistrationEvaluation::ICPPairwiseRegistration::ICPPairwiseRegistration(const float maxSearchRadius, const double maxRmse, const int numMaxIter, const float samplingPerc, const bool verbose, const bool calcScale, const bool calcRMSE)
:PairwiseRegistration(calcScale, calcRMSE, verbose),
m_polyData1BoundaryPointsIds(NULL),
m_polyData1BoundaryPointsIdsAssigned(false),
m_polyData1KdTree(NULL),
m_polyData1KdTreeAssigned(false),
m_numMaxIter(numMaxIter),
m_maxRMSE(maxRmse),
m_maxSearchRadius(maxSearchRadius),
m_samplingPerc(samplingPerc)
{
}

Pairwise3DRegistrationEvaluation::ICPPairwiseRegistration::~ICPPairwiseRegistration()
{
	if(!m_polyData1BoundaryPointsIdsAssigned)
	{
		if(m_polyData1BoundaryPointsIds)
		{
			delete[] m_polyData1BoundaryPointsIds;
			m_polyData1BoundaryPointsIds = NULL;
		}
	}
	if(!m_polyData1KdTreeAssigned)
	{
		if(m_polyData1KdTree)
		{
			delete m_polyData1KdTree;
			m_polyData1KdTree = NULL;
		}
	}
}

void Pairwise3DRegistrationEvaluation::ICPPairwiseRegistration::SetPolyData1BoundaryPointsIds(bool* polyData1BoundaryPointsIds)
{
	if(!m_polyData1BoundaryPointsIdsAssigned)
	{
		if(m_polyData1BoundaryPointsIds)
		{
			delete[] m_polyData1BoundaryPointsIds;
			m_polyData1BoundaryPointsIds = NULL;
		}
	}
	m_polyData1BoundaryPointsIds = polyData1BoundaryPointsIds?polyData1BoundaryPointsIds:NULL;
	m_polyData1BoundaryPointsIdsAssigned = polyData1BoundaryPointsIds?true:false;
}

void Pairwise3DRegistrationEvaluation::ICPPairwiseRegistration::SetPolyData1KdTree(KdTree* polyData1KdTree)
{
	if(!m_polyData1KdTreeAssigned)
	{
		if(m_polyData1KdTree)
		{
			delete[] m_polyData1KdTree;
			m_polyData1KdTree = NULL;
		}
	}
	m_polyData1KdTree = polyData1KdTree?polyData1KdTree:NULL;
	m_polyData1KdTreeAssigned = polyData1KdTree?true:false;
}

vtkTransform* Pairwise3DRegistrationEvaluation::ICPPairwiseRegistration::Register(vtkPolyData* polyData1In, vtkPolyData* polyData2In)
{
	if(m_transform)
	{
		m_transform->Delete();
	}
	m_transform = vtkTransform::New();
	m_transform->PostMultiply();

	vtkPolyData* poly1_sampl;
	vtkPolyData* poly2_sampl;
	if(m_samplingPerc < 1.0f)
	{
		poly1_sampl = CreatePolyData_RandomSampling(polyData1In, m_samplingPerc);
		poly2_sampl = CreatePolyData_RandomSampling(polyData2In, m_samplingPerc);
	}
	else
	{
		poly1_sampl = polyData1In;
		poly2_sampl = Duplicate(polyData2In);
	}

	if(!m_polyData1BoundaryPointsIdsAssigned)
	{
		if(m_polyData1BoundaryPointsIds)
		{
			delete[] m_polyData1BoundaryPointsIds;
			m_polyData1BoundaryPointsIds = NULL;
		}
		GetBoundaryPoints(poly1_sampl, m_polyData1BoundaryPointsIds);
	}

	if(!m_polyData1KdTreeAssigned)
	{
		if(m_polyData1KdTree)
		{
			delete m_polyData1KdTree;
			m_polyData1KdTree = NULL;
		}
		m_polyData1KdTree = new KdTree(poly1_sampl);
	}


	int nPoints2_sampl = poly2_sampl->GetNumberOfPoints();

	int nearestPointIdx;
	int numIter = 0;
	double rmse = std::numeric_limits<double>::max();
	int numICPpoints;
	float* ICPpoints1 = new float[nPoints2_sampl*3];
	float* ICPpoints2 = new float[nPoints2_sampl*3];
	float* ptrICPpoints1;
	float* ptrICPpoints2;

	float* ptrPoly1_sampl = GetPolyDataPointsPointer(poly1_sampl);
	float* ptrPoly2_sampl = GetPolyDataPointsPointer(poly2_sampl);


	HornOrientation horn(false, true);
	double* matrTransform;

	double dist;
	while((numIter < m_numMaxIter)&&(rmse > m_maxRMSE))
	{
		numICPpoints = 0;
		ptrICPpoints1 = ICPpoints1;
		ptrICPpoints2 = ICPpoints2;
		float* ptrptrPoly2_sampl = ptrPoly2_sampl;
		for(int po=0;po<nPoints2_sampl; po++)
		{
			nearestPointIdx = m_polyData1KdTree->FindNearestPointWithinRadius( ptrptrPoly2_sampl, m_maxSearchRadius, dist );
			if((nearestPointIdx != -1) && (!m_polyData1BoundaryPointsIds[nearestPointIdx]))
			{
				memcpy(ptrICPpoints1, ptrPoly1_sampl + 3*nearestPointIdx, sizeof(float)*3);
				memcpy(ptrICPpoints2, ptrPoly2_sampl + 3*po, sizeof(float)*3);
				ptrICPpoints1+=3;
				ptrICPpoints2+=3;
				numICPpoints++;
			}
			ptrptrPoly2_sampl+=3;
		}
		if(numICPpoints>=3)
		{
			int ret = horn.Orientate(ICPpoints1, ICPpoints2, numICPpoints, matrTransform);
			if(ret < 0)
			{
				delete[] ICPpoints1;
				delete[] ICPpoints2;

				m_rmse = std::numeric_limits<double>::max();
			
				return NULL;
			}

			rmse =  horn.GetRMSE();
		}
		else
		{
			delete[] ICPpoints1;
			delete[] ICPpoints2;

			m_rmse = std::numeric_limits<double>::max();
			
			return NULL;
		}
		ApplyTransform(ptrPoly2_sampl, nPoints2_sampl, matrTransform);
		m_transform->Concatenate(matrTransform);

		numIter++;
	}

	m_rmse = rmse;


	if(m_samplingPerc < 1.0f)
	{
		poly1_sampl->Delete();
		poly2_sampl->Delete();
	}
	else
	{
		poly2_sampl->Delete();
	}

	delete[] ICPpoints1;
	delete[] ICPpoints2;

	return m_transform;
}


Pairwise3DRegistrationEvaluation::PairwiseRegistration::PairwiseRegistration(bool calcScale, bool calcRMSE, bool verbose)
	:m_calcRMSE(calcRMSE)
	,m_calcScale(calcScale)
	,m_verbose(verbose)
{
	m_transform = vtkTransform::New();
	m_rmse = -1.0;
}

double Pairwise3DRegistrationEvaluation::PairwiseRegistration::CalcRMSE(float* points1, float* points2, int nPoints)
{
	double rmse = calcDistance2BetweenPoints(points1, points2, nPoints);
	rmse /=(double)nPoints;
	rmse = sqrt(rmse);

	return rmse;
}

double Pairwise3DRegistrationEvaluation::PairwiseRegistration::CalcRMSE(vtkPolyData* polyData1, vtkPolyData* polyData2)
{
	if(polyData1->GetNumberOfPoints() != polyData2->GetNumberOfPoints())
	{
		return -1.0;
	}

	return CalcRMSE( (float*)polyData1->GetPoints()->GetVoidPointer(0), (float*)polyData2->GetPoints()->GetVoidPointer(0), polyData1->GetNumberOfPoints() );
}

Pairwise3DRegistrationEvaluation::PairwiseRegistration::~PairwiseRegistration()
{
	if(m_transform)
	{
		m_transform->Delete();
		m_transform = NULL;
	}
}



Pairwise3DRegistrationEvaluation::HornOrientation::HornOrientation(bool calcScale, bool calcRMSE)
	:AbsoluteOrientation(calcScale, calcRMSE)
{
}

int Pairwise3DRegistrationEvaluation::HornOrientation::Orientate(float *points1, float* points2, int nPoints, double* &matrTransform)
{
	if(nPoints < 3)
	{
		return -1;
	}

	float* points1Out = new float[nPoints*3];
	float* points2Out = new float[nPoints*3];
	memcpy(points1Out, points1, sizeof(float)*nPoints*3);
	memcpy(points2Out, points2, sizeof(float)*nPoints*3);

	float barycenter1[3];
	float barycenter2[3];
	GetCentroid(points1Out, nPoints, barycenter1);
	GetCentroid(points2Out, nPoints, barycenter2);

	//translate points in barycenters
	float invbarycenter1[3];
	float invbarycenter2[3];
	invbarycenter1[0] = -barycenter1[0];
	invbarycenter1[1] = -barycenter1[1];
	invbarycenter1[2] = -barycenter1[2];
	invbarycenter2[0] = -barycenter2[0];
	invbarycenter2[1] = -barycenter2[1];
	invbarycenter2[2] = -barycenter2[2];
	Translate(points1Out, nPoints, invbarycenter1);
	Translate(points2Out, nPoints, invbarycenter2);

	float* ptr1;
	float* ptr2;

	//calc scale
	float scale=1.0;
	if(m_calcScale)
	{
		double sum1 = 0.0;
		double sum2 = 0.0;
		ptr1 = points1Out;
		ptr2 = points2Out;
		for(int po=0; po<nPoints; po++)
		{
			sum1 += sqrt( (ptr1[0] * ptr1[0]) + (ptr1[1] * ptr1[1]) + (ptr1[2] * ptr1[2]) );
			sum2 += sqrt( (ptr2[0] * ptr2[0]) + (ptr2[1] * ptr2[1]) + (ptr2[2] * ptr2[2]) );
			ptr1+=3;
			ptr2+=3;
		}
		if(!IsZero( sum2 ))
		{
			scale = sum1 / sum2;
		}
	}

	//calc rotation
	float Sxx=0, Sxy=0, Sxz=0;
	float Syx=0, Syy=0, Syz=0;
	float Szx=0, Szy=0, Szz=0;
	ptr1 = points1Out;
	ptr2 = points2Out;
	for (int po=0; po<nPoints ;po++)
	{
		Sxx += ptr1[0] * ptr2[0];
		Sxy += ptr1[0] * ptr2[1];
		Sxz += ptr1[0] * ptr2[2];
		Syx += ptr1[1] * ptr2[0];
		Syy += ptr1[1] * ptr2[1];
		Syz += ptr1[1] * ptr2[2];
		Szx += ptr1[2] * ptr2[0];
		Szy += ptr1[2] * ptr2[1];
		Szz += ptr1[2] * ptr2[2];
		ptr1 +=3;
		ptr2 +=3;
	}

	////calc N
	//double n[] = { Sxx + Syy + Szz,  Syz - Szy      , Szx - Sxz         , Sxy - Syx,
	//	Syz - Szy      , Sxx - Syy - Szz , Sxy + Syx         , Szx + Sxz,
	//	Szx - Sxz      , Sxy + Syx       , - Sxx + Syy - Szz , Syz + Szy,
	//	Sxy - Syx      , Szx + Sxz       , Syz + Szy         , - Sxx - Syy + Szz };

	//double evect[16];
	//double evalues[] = {0, 0, 0, 0};

	//CvMat N, EigenVectors, Eigenvalues;

	//cvInitMatHeader( &N           , 4, 4, CV_64FC1, n, 2147483647);
	//cvInitMatHeader( &EigenVectors, 4, 4, CV_64FC1, evect, 2147483647);
	//cvInitMatHeader( &Eigenvalues , 4, 1, CV_64FC1, evalues, 2147483647);

	//cvEigenVV(&N, &EigenVectors, &Eigenvalues, 0 );




	//calc N
	double *n[4];
	n[0] = new double[4];
	n[1] = new double[4];
	n[2] = new double[4];
	n[3] = new double[4];
	n[0][0] = Sxx + Syy + Szz;
	n[0][1] = Syz - Szy;
	n[0][2] = Szx - Sxz;
	n[0][3] = Sxy - Syx;
	n[1][0] = Syz - Szy;
	n[1][1] = Sxx - Syy - Szz;
	n[1][2] = Sxy + Syx;
	n[1][3] = Szx + Sxz;
	n[2][0] = Szx - Sxz;
	n[2][1] = Sxy + Syx;
	n[2][2] = - Sxx + Syy - Szz;
	n[2][3] = Syz + Szy;
	n[3][0] = Sxy - Syx;
	n[3][1] = Szx + Sxz;
	n[3][2] = Syz + Szy;
	n[3][3] = - Sxx - Syy + Szz;

	double *evect[4];
	evect[0] = new double[4];
	evect[1] = new double[4];
	evect[2] = new double[4];
	evect[3] = new double[4];
	memset(evect[0], 0.0, sizeof(double)*4);
	memset(evect[1], 0.0, sizeof(double)*4);
	memset(evect[2], 0.0, sizeof(double)*4);
	memset(evect[3], 0.0, sizeof(double)*4);

	double evalues[] = {0, 0, 0, 0};

	int ret = vtkMath::JacobiN(n, 4, evalues, evect);
	if(ret == 0)
	{
		delete[] points1Out;
		delete[] points2Out;

		for (int i = 0; i < 4; i++)
		{
			delete [] evect[i];
			delete [] n[i];
		}

		return -1;
	}

	//CvMat N, EigenVectors, Eigenvalues;

	//cvInitMatHeader( &N           , 4, 4, CV_64FC1, n, 2147483647);
	//cvInitMatHeader( &EigenVectors, 4, 4, CV_64FC1, evect, 2147483647);
	//cvInitMatHeader( &Eigenvalues , 4, 1, CV_64FC1, evalues, 2147483647);

	//cvEigenVV(&N, &EigenVectors, &Eigenvalues, 0 );



	//calc final transform
	//initialize last row
	for(int co=0; co<3; co++)
	{
		m_matrTransform[3*4 + co] = 0;
	}
	m_matrTransform[3*4 + 3] = 1;

	//translate to barycenter1
	for(int ro=0; ro<3; ro++)
	{
		m_matrTransform[ro*4+3] = barycenter1[ro];
	}

	//first eigenvector is the rotation quaternion
	//scale and rotate
	double rotMatr[3][3];
	//evect[0] = - evect[0];
	//vtkMath::QuaternionToMatrix3x3(evect, rotMatr);

	//as  JacobiN returns evect in column wise:
	double maxEvect[] = { - evect[0][0], evect[1][0], evect[2][0], evect[3][0]};
	vtkMath::QuaternionToMatrix3x3(maxEvect, rotMatr);

	for(int ro=0; ro<3; ro++)
	{
		for(int co=0; co<3; co++)
		{
			m_matrTransform[ro*4+co] = rotMatr[ro][co] * scale;
		}
	}

	//translate to invbarycenter2 (It's like m_transform->Translate(invbarycenter2) )
	double matrixTrasl[16];
	vtkMatrix4x4::Identity(matrixTrasl);
	matrixTrasl[4*0+3] = invbarycenter2[0];
	matrixTrasl[4*1+3] = invbarycenter2[1];
	matrixTrasl[4*2+3] = invbarycenter2[2];
	vtkMatrix4x4::Multiply4x4(m_matrTransform, matrixTrasl, m_matrTransform);


	if(m_calcRMSE)
	{
		Translate(points2Out, nPoints, barycenter2);
		Translate(points1Out, nPoints, barycenter1);
		ApplyTransform(points2Out, nPoints, m_matrTransform);
		m_rmse = CalcRMSE(points1Out, points2Out, nPoints);
	}

	delete[] points1Out;
	delete[] points2Out;

	for (int i = 0; i < 4; i++)
	{
		delete [] evect[i];
		delete [] n[i];
	}

	matrTransform = m_matrTransform;
	return 0;

}



Pairwise3DRegistrationEvaluation::AbsoluteOrientation::AbsoluteOrientation(bool calcScale, bool calcRMSE)
	:m_calcRMSE(calcRMSE)
	,m_calcScale(calcScale)
{
	m_transform = vtkTransform::New();
	m_rmse = -1.0;
}

vtkTransform* Pairwise3DRegistrationEvaluation::AbsoluteOrientation::Orientate(vtkPolyData* polyData1In, vtkPolyData* polyData2In )
{
	if(polyData1In->GetNumberOfPoints() != polyData2In->GetNumberOfPoints())
	{
		if(m_transform)
		{
			m_transform->Delete();
			m_transform = NULL;
		}
		return NULL;
	}

	double* ptrMatr = NULL;
	int ret = Orientate( (float*)(polyData1In->GetPoints()->GetVoidPointer(0)), (float*)(polyData2In->GetPoints()->GetVoidPointer(0)), polyData1In->GetNumberOfPoints(), ptrMatr);

	if(m_transform)
	{
		m_transform->Delete();
		m_transform = NULL;
	}
	if(ret>=0)
	{
		m_transform = vtkTransform::New();
		m_transform->SetMatrix(ptrMatr);
	}

	return m_transform;
}

double Pairwise3DRegistrationEvaluation::AbsoluteOrientation::CalcRMSE(float* points1, float* points2, int nPoints)
{
	double rmse = 0;
	float* ptrPoints1 = points1;
	float* ptrPoints2 = points2;
	for(int po=0; po<nPoints; po++)
	{
		rmse += vtkMath::Distance2BetweenPoints(ptrPoints1, ptrPoints2);
		ptrPoints1+=3;
		ptrPoints2+=3;
	}
	rmse /=(double)nPoints;
	rmse = sqrt(rmse);

	return rmse;
}

double Pairwise3DRegistrationEvaluation::AbsoluteOrientation::CalcRMSE(vtkPolyData* polyData1, vtkPolyData* polyData2)
{
	if(polyData1->GetNumberOfPoints() != polyData2->GetNumberOfPoints())
	{
		return -1.0;
	}

	return CalcRMSE( (float*)polyData1->GetPoints()->GetVoidPointer(0), (float*)polyData2->GetPoints()->GetVoidPointer(0), polyData1->GetNumberOfPoints() );
}


Pairwise3DRegistrationEvaluation::AbsoluteOrientation::~AbsoluteOrientation()
{
	if(m_transform)
	{
		m_transform->Delete();
		m_transform = NULL;
	}
}



