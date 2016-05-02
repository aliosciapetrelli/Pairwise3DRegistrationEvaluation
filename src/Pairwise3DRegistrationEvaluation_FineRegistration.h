#ifndef PAIRWISE_3D_REGISTRATION_EVALUATION_FINE_REGISTRATION_H
#define PAIRWISE_3D_REGISTRATION_EVALUATION_FINE_REGISTRATION_H


#include "Pairwise3DRegistrationEvaluation_utils.h"

#include "vtkTransform.h"
#include "vtkPolyData.h"

namespace Pairwise3DRegistrationEvaluation
{
	class PairwiseRegistration
	{
	protected:
		vtkTransform*		m_transform;
		double				m_rmse;
		bool				m_calcRMSE;
		bool				m_calcScale;
		bool				m_verbose;

	public:

		PairwiseRegistration(bool calcScale = true, bool calcRMSE = false, bool verbose = false);
		virtual ~PairwiseRegistration();

		vtkTransform*		GetTransform(){return m_transform;}

		vtkTransform*	Register(vtkPolyData* polyData1In, vtkPolyData* polyData2In){return m_transform;};

		static double CalcRMSE(float* points1, float* points2, int nPoints);
		static double CalcRMSE(vtkPolyData* polyData1, vtkPolyData* polyData2);

		double GetRMSE(){return m_rmse;}
		void CalcScale(bool calcScale){m_calcScale = calcScale;}
		void CalcRMSE(bool calcRMSE){m_calcRMSE = calcRMSE;}
		bool GetCalcScale(){return m_calcScale;}
		bool GetCalcRMSE(){return m_calcRMSE;}
	};


	class AbsoluteOrientation
	{
	protected:
		vtkTransform*		m_transform;
		double				m_matrTransform[16];
		double				m_rmse;
		bool				m_calcRMSE;
		bool				m_calcScale;

	public:
		AbsoluteOrientation(bool calcScale = true, bool calcRMSE = false);
		virtual				~AbsoluteOrientation();

		vtkTransform*		GetTransform(){return m_transform;}
		double*				GetMatrTransform(){return m_matrTransform;}

		//every polydata must be allocated, vtkTransform mustn't be deallocated
		virtual vtkTransform*		Orientate(vtkPolyData* polyData1In, vtkPolyData* polyData2In);

		virtual int Orientate(float *points1, float* points2, int nPoints, double* &matrTransform) = 0;

		static double CalcRMSE(float* points1, float* points2, int nPoints);
		static double CalcRMSE(vtkPolyData* polyData1, vtkPolyData* polyData2);

		double GetRMSE(){return m_rmse;}
		void CalcScale(bool calcScale){m_calcScale = calcScale;}
		void CalcRMSE(bool calcRMSE){m_calcRMSE = calcRMSE;}
		bool GetCalcScale(){return m_calcScale;}
		bool GetCalcRMSE(){return m_calcRMSE;}
	};



	class HornOrientation : public AbsoluteOrientation
	{
	public:
		HornOrientation(bool calcScale = true, bool calcRMSE = false);

		virtual int Orientate(float *points1, float* points2, int nPoints, double* &matrTransform);

	};



	class ICPPairwiseRegistration : public PairwiseRegistration
	{
	protected:
		bool*							m_polyData1BoundaryPointsIds;
		bool							m_polyData1BoundaryPointsIdsAssigned;

		int								m_numMaxIter;
		double							m_maxRMSE;
		float							m_maxSearchRadius;
		float							m_samplingPerc;

		KdTree*							m_polyData1KdTree;
		bool							m_polyData1KdTreeAssigned;

	public:
		ICPPairwiseRegistration(const float maxSearchRadius = 0.2f, const double maxRmse = 0.0001, const int numMaxIter = 10, const float samplingPerc = 1.0, const bool verbose = false, const bool calcScale = false, const bool calcRMSE = false);
		virtual ~ICPPairwiseRegistration();
		virtual vtkTransform*	Register(vtkPolyData* polyData1In, vtkPolyData* polyData2In);
		virtual void SetPolyData1BoundaryPointsIds(bool* polyData1BoundaryPointsIds);
		virtual void SetPolyData1KdTree(KdTree* polyData1KdTree);
		virtual void SetNumMaxIter(int numMaxIter){m_numMaxIter = numMaxIter;}
		virtual void SetMaxRMSE(double maxRmse){m_maxRMSE = maxRmse;}
		virtual void SetSamplingPerc(float samplingPerc){m_samplingPerc = samplingPerc;}
		virtual void SetMaxSearchRadius(double maxSearchRadius){m_maxSearchRadius = maxSearchRadius;}

	};


	class GICPPairwiseRegistration : public PairwiseRegistration
	{
	protected:
		int								m_numMaxIter;
		double							m_epsilon;
		float							m_maxSearchRadius;

	public:
		GICPPairwiseRegistration(int numMaxIter = 5, double epsilon = 0.0001, float maxSearchRadius = 0.2, bool verbose = false);
		virtual vtkTransform*	Register(vtkPolyData* polyData1In, vtkPolyData* polyData2In);
		virtual void SetNumMaxIter(int numMaxIter){m_numMaxIter = numMaxIter;}
		virtual void SetEpsilon(double epsilon){m_epsilon = epsilon;}
		virtual void SetMaxSearchRadius(double maxSearchRadius){m_maxSearchRadius = maxSearchRadius;}
	};



}

#endif



