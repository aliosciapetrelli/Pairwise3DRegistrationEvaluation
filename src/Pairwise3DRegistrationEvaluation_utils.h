#ifndef PAIRWISE_3D_REGISTRATION_EVALUATION_UTILS_H
#define PAIRWISE_3D_REGISTRATION_EVALUATION_UTILS_H


#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <stdint.h>

#include "vtkPolyData.h"
#include "vtkAbstractTransform.h"
#include "vtkPolyDataCollection.h"
#include "vtkTransform.h"
#include "vtkAbstractPointLocator.h"


#define BOOST_IS_NOT_INCLUDED

#ifdef _MSC_VER
	#ifndef _CRT_SECURE_NO_WARNINGS
		#define _CRT_SECURE_NO_WARNINGS
	#endif
	#define NOMINMAX
	#include "windows.h"
#endif

#ifdef __GNUC__
    #include <time.h>
    #include <sys/times.h>
#endif // __GNUC__


namespace Pairwise3DRegistrationEvaluation
{

	class CpuTimeProfiler
	{
	protected:
#if defined (_MSC_VER)
		LARGE_INTEGER m_frequency;

		ULARGE_INTEGER m_actualTime;
		FILETIME m_actualUserTime;
		FILETIME m_actualKernelTime;
		FILETIME m_actualCreateTime;
		FILETIME m_actualExitTime;
#elif defined __GNUC__
		clockid_t m_clk_id;
		timespec m_actualTime;
#endif

	public:
#ifdef __GNUC__
		/// Typical choices: CLOCK_MONOTONIC, CLOCK_PROCESS_CPUTIME_ID, CLOCK_THREAD_CPUTIME_ID
		CpuTimeProfiler(clockid_t clk_id = CLOCK_PROCESS_CPUTIME_ID);
#else
		CpuTimeProfiler();
#endif
		double GetElapsedHours();
		double GetElapsedMins();
		double GetElapsedSecs();
		double GetElapsedMilli();
		void GetActualTime();
	};



	/*!
	* \brief Compute minimum.
	*
	* \param a,b instances to compare.
	*/
	template <typename Type> inline Type Min(Type a, Type b) { return (a <= b)? a : b; };

	/*!
	* \brief Compute maximum.
	*
	* \param a,b instances to compare.
	*/
	template <typename Type> inline Type Max(Type a, Type b) { return (a >= b)? a : b; };



	std::string ChangeExtension(const std::string &filename, const std::string &newExt);

	bool ExistsDir(const std::string &directory);
	bool ExistFiles(const std::vector<std::string> &absFiles);
	bool ExistsFile(const std::string &filename_abs);

	std::string RemoveExt(const std::string &filename);

	std::string GetRelativeName(const std::string &filename);
	std::string GetRelativeName_WithoutExt(const std::string &filename);

	vtkPolyData* Duplicate(vtkPolyData* source);

	void ApplyTransform(vtkPolyData *polyData, vtkAbstractTransform* transform);
	void ApplyTransform(float* points, int nPoints, const double* matrTransf);

	double calcDistance2BetweenPoints(float* points1, float* points2, int nPoints);

	double CalcRMSE(float* points1, float* points2, int nPoints, double meshRes);
	double CalcRMSE(vtkPolyData* polyData1, vtkPolyData* polyData2, double meshRes);

	/*!
	* \brief Compute squared euclidean distance between two 3D points.
	*
	* \param point1 3D point
	* \param point2 3D point
	* \return squared euclidean distance
	*/
	inline float EuclideanDistance2_3D(float* point1, float* point2){return (point1[0]-point2[0])*(point1[0]-point2[0]) + (point1[1]-point2[1])*(point1[1]-point2[1]) + (point1[2]-point2[2])*(point1[2]-point2[2]);  };

	/*!
	* \brief Compute squared euclidean distance between two 3D points.
	*
	* \param point1 3D point
	* \param point2 3D point
	* \return squared euclidean distance
	*/
	inline double EuclideanDistance2_3D(double* point1, double* point2){return (point1[0]-point2[0])*(point1[0]-point2[0]) + (point1[1]-point2[1])*(point1[1]-point2[1]) + (point1[2]-point2[2])*(point1[2]-point2[2]);  };

	/*!
	* \brief Compute euclidean distance between two 3D points.
	*
	* \param point1 3D point
	* \param point2 3D point
	* \return euclidean distance
	*/
	inline float EuclideanDistance_3D(float* point1, float* point2){ return sqrt(EuclideanDistance2_3D(point1, point2)); };

	/*!
	* \brief Compute euclidean distance between two 3D points.
	*
	* \param point1 3D point
	* \param point2 3D point
	* \return euclidean distance
	*/
	inline double EuclideanDistance_3D(double* point1, double* point2){ return sqrt(EuclideanDistance2_3D(point1, point2)); };

	vtkPolyData* CreatePolyData(const float* points, const vtkIdType nPoints, const unsigned char* colors = NULL, const vtkIdType* triangles = NULL, const vtkIdType nTriangles = 0);

	vtkPolyData* CreatePolyData(vtkPolyData* polyData, vtkIdType* pointIds, vtkIdType nPointIds, const bool insertNormals = true, const bool insertColors = true, const bool insertTriangles = true);
	vtkPolyData* CreatePolyData_RandomSampling(vtkPolyData* cloud, const float samplingPerc);

	void CalcMeshNormals(vtkPolyData *polyData, const double dFeatureAngle = -1.0);
	void CalcMeshNormals(vtkPolyDataCollection *polyDataCollection, const double dFeatureAngle = -1.0);

	void RemoveNotPolysCells(vtkPolyData *polyData);

	void CleanPolyData(vtkPolyData *polyData, const bool mergePoints = true, const double tolerance = 0.0, bool removeNotPolysCells = true);
	void CleanPolyData(vtkPolyDataCollection *polyDataCollection, const bool mergePoints = true, const double tolerance = 0.0, bool removeNotPolysCells = true);


	float* GetPolyDataPointsPointer(vtkPolyData* polyData, const unsigned int iPoint=0);
	float* GetPolyDataPointNormalsPointer(vtkPolyData* polyData, const unsigned int iNormal=0);
	vtkIdType* GetPolyDataTrianglesPointer(vtkPolyData* polyData, const unsigned int iTriangle=0);
	unsigned char* GetPolyDataColorsPointer(vtkPolyData* polyData, const unsigned int iColor=0);
	vtkIdType GetPolyDataNumberOfPoints(vtkPolyData* polyData);

	float* AllocateNormals(vtkPolyData *polyData);
	float* AllocatePoints(vtkPolyData *polyData, const vtkIdType nPoints);
	unsigned char* AllocateColors(vtkPolyData *polyData, unsigned char* color = NULL, bool* mask = NULL);
	vtkIdType* AllocateTriangles(vtkPolyData *polyData, const vtkIdType nTriangles);

	vtkPolyData* ReadMesh(const std::string &absMeshFileName, bool convertToVtk = false);
	vtkPolyDataCollection* ReadMeshes(const std::vector<std::string> &vAbsMeshFileNames, bool convertToVtk = false );
	vtkPolyDataCollection* ReadMeshes(const std::string &meshesPath, std::vector<std::string> &vAbsMeshFileNames, const std::string meshExt = "ply", bool convertToVtk = false );

	vtkPolyData* LoadPolyData(const std::string &absFilename);
	vtkPolyDataCollection* LoadPolyData(const std::vector<std::string> &filenames);

	vtkPolyData* LoadPly(const char* filename);
	vtkPolyData* LoadOff(const char* filename);
	vtkPolyData* LoadOBJ(const char* filename);
	vtkPolyData* LoadSTL(const char* filename);
	vtkPolyData* LoadVtk(const char* filename);
	vtkPolyData* LoadVertTri(const char* filenameWithoutExtension);
	vtkPolyData* LoadPcd(const char* filename);
	vtkPolyData* LoadSimpleTxt(const char* filename);
	vtkPolyData* LoadXYZ(const char* fileName);

	bool WriteVtk(const std::string &absFileName, vtkPolyData* polyData, const bool fileTypeASCII = false );
	bool WriteVtk(const std::vector<std::string> &absFiles, vtkPolyDataCollection* polyDataCollection, const bool fileTypeASCII = false );

	void Preprocessing(vtkPolyData* &polyData);
	void Preprocessing(vtkPolyDataCollection* &polyDataCollection);

	template<typename T> T LexicalCast(const std::string& s)
	{
		std::stringstream ss(s);

		T result;
		if ((ss >> result).fail() || !(ss >> std::ws).eof())
		{
			//throw std::bad_cast();
			cout << "ERROR:Impossible to cast " << s;
			getchar();
			exit(-1);
		}

		return result;
	}

	//logical sorting (e.g. Windows explorer)
	class StringCompare_Smart_Incr
	{
	public:
		inline bool operator() (const std::string& a, const std::string& b) const
		{
			unsigned posStr = 0;
			while( (posStr < a.size() ) && (posStr < b.size()) )
			{
				unsigned tkn_idx_a = (unsigned int) a.find_first_of("0123456789", posStr);
				unsigned tkn_idx_b = (unsigned int) b.find_first_of("0123456789", posStr);
				std::string suba = a.substr(posStr, tkn_idx_a - posStr );
				std::string subb = b.substr(posStr, tkn_idx_b - posStr );
				if(suba == subb)
				{
					//same substring

					if(tkn_idx_a == a.size())
					{
						//end of a and at least of b
						return true;
					}

					if(tkn_idx_a == a.size())
					{
						//end of b but not of a
						return false;
					}

					unsigned numberEnd_a = (unsigned int) a.find_first_not_of("0123456789", tkn_idx_a+1);
					unsigned numberEnd_b = (unsigned int) b.find_first_not_of("0123456789", tkn_idx_b+1);
					//check number
					long long number_a = LexicalCast<long long>(a.substr(tkn_idx_a, numberEnd_a - tkn_idx_a));
					long long number_b = LexicalCast<long long>(b.substr(tkn_idx_b, numberEnd_b - tkn_idx_b));
					//long number_a = std::atol(a.substr(tkn_idx_a).c_str());
					//long number_b = std::atol(b.substr(tkn_idx_b).c_str());
					if(number_a != number_b)
					{
						return (number_a < number_b);
					}
				}
				else
				{
					//different substring
					return (suba < subb);
				}
				posStr = (unsigned int) a.find_first_not_of("0123456789", tkn_idx_a + 1);
			}

			return ( a.size() < b.size() );
		}
	};

	bool FindFilesEndWith(const std::string & path, const std::string & endingString, std::vector<std::string> & foundFiles, bool getOnlyFileName = false, const int nMaxItems = std::numeric_limits<int>::max());

	double ComputeMeshResolution(vtkPolyData* cloud);
	double ComputeMeshResolution(vtkPolyDataCollection* polyDataCollection, const std::string &absMeshResFile = "");

	template<typename T> void SaveSingleValue(const std::string &absFileName, T value, const bool binary = true)
	{
		if(binary)
		{
			fstream binary_file;

			binary_file.open(absFileName, ios::out|ios::binary);

			if (!binary_file.is_open() )
			{
				cout << "ERROR: Impossible to open file " << absFileName << "\n";
				getchar();
				return;
			}

			binary_file.write((char*)(&value), sizeof(value));
			binary_file.close();

		}
		else
		{

			ofstream file;
			file.open(absFileName, ios_base::out);

			if (!file.is_open() )
			{
				cout << "ERROR: Impossible to open file " << absFileName << "\n";
				getchar();
				return;
			}

			file << value;

			file.close();
		}
	}

	template<typename T> void LoadSingleValue(const std::string &absFileName, T &value, const bool binary = true)
	{
		if(binary)
		{
			fstream binary_file(absFileName, ios::binary|ios::in);

			if (!binary_file.is_open() )
			{
				cout << "ERROR: Impossible to open file " << absFileName << "\n";
				getchar();
				return;
			}

			binary_file.read((char*)(&value), sizeof(value) );

			binary_file.close();
		}
		else
		{

			ifstream file; 
			file.open(absFileName, ios_base::in);

			if (!file.is_open() )
			{
				cout << "ERROR: Impossible to open file " << absFileName << "\n";
				getchar();
				return;
			}	

			file >> value;

			file.close();
		}
	}

	void LoadGroundTruth(const std::string &strGTFile, std::vector<std::string> &vViewFileNames, std::vector<vtkTransform*> &vTransforms);
	void LoadGroundTruth(const std::string &strGTFile, std::vector<std::string> &vMeshFileNames, std::vector<vtkTransform*> &vTransforms, vtkPolyDataCollection* meshes );

	vtkTransform* CreateInverseTransform(vtkTransform* transform);
	void DeleteTransforms(std::vector<vtkTransform*> &transforms);

	vtkTransform* ComputeGroundTruth(vtkTransform* &transf_trg, vtkTransform* &transf_ref);

	void ComputeGroundTruth_All(std::vector<vtkTransform*> &vTransforms_GT, const std::vector<std::pair<int,int>> &vMeshPairs, std::vector<vtkTransform*> &vPairwiseTransforms_GT);

	void CalcAllViewPairs(std::vector<std::pair<int, int> > &vMatchPairs, int nElements);


	class KdTree
	{
	protected:
		vtkAbstractPointLocator*	m_kdTree;
		vtkIdList*				m_pointsList;
		vtkPolyData*			m_polyData;

	public:
		KdTree(const bool useTrueKdTree = false);
		KdTree(vtkPolyData* polyData, const bool useTrueKdTree = false);
		~KdTree();

		void SetPolyData(vtkPolyData* polyData);
		vtkPolyData* GetPolyData(){return m_polyData;};

		vtkAbstractPointLocator* GetKdTree(){return m_kdTree;}
		vtkIdList* GetFoundPoints(){return m_pointsList;}

		//FindPointsWithinRadius
		vtkIdList* FindPointsWithinRadius(const float* const point, float radius);
		vtkIdList* FindPointsWithinRadius(const double* point, double radius);
		/*!
		* \param pointIndex index of m_polydata point.
		*/
		vtkIdList* FindPointsWithinRadius(const int pointIndex, float radius);

		//FindNearestPoint
		int FindNearestPoint(double* point);
		int FindNearestPoint(float* point);
		/*!
		* \param pointIndex index of m_polydata point.
		*/
		int FindNearestPoint(const int pointIndex);

		//FindNearestPoint, return distance
		double FindNearestPoint(double* point, int &nearestPointId);
		double FindNearestPoint(float* point, int &nearestPointId);
		/*!
		* \param pointIndex index of m_polydata point.
		*/
		double FindNearestPoint(const int pointIndex, int &nearestPointId);

		//FindNearestPointWithinRadius
		int FindNearestPointWithinRadius(double* point, double radius, double & dist);
		/*!
		* \param point
		* \param radius
		* \param dist output parameter
		* \return the index of the nearest point in the mesh of the polydata, -1 if not found
		*/
		int FindNearestPointWithinRadius(float* point, float radius, double & dist);
		/*!
		* \param pointIndex index of m_polydata point.
		*/
		int FindNearestPointWithinRadius(const int pointIndex, float radius, double & dist);

		//FindNearestNPoints
		vtkIdList* FindNearestNPoints(int n, const double* point);
		vtkIdList* FindNearestNPoints(int n, const float* point);
		vtkIdList* FindNearestNPoints(int n, const int pointIndex);
	};


	int GetBoundaryPoints(vtkPolyData *polydata, bool* &boundaryPointsIds);

	void GetCentroid(float* points, int nPoints, float centroid[3]);
	void Translate(float* points, int nPoints, const float point[3]);

	const float zeroFloatEps8 = 1E-8f;

	inline bool IsZero(float val, float zeroFloatEps = zeroFloatEps8){return (abs(val)<zeroFloatEps);};


	double CalcOverlappingAreaMax(vtkPolyData* polyData1, vtkPolyData* polyData2, double maxDistance, double absNormalDotThresh = 0.9, bool getOverlappedPoints = false, bool* overlappedPoints1 = NULL, bool* overlappedPoints2 = NULL);

	inline int NextRandom(int i){ return std::rand()%i;}

	size_t Sample_random(std::vector<bool> &vData, const float samplingFactor, const unsigned int seed = std::numeric_limits<unsigned int>::max());
}

#endif



