#ifndef PAIRWISE_3D_REGISTRATION_EVALUATION
#define PAIRWISE_3D_REGISTRATION_EVALUATION

#include "Pairwise3DRegistrationEvaluation_utils.h"

/** \brief Framework for the evaluation of pairwise 3D registration algorithms as presented in :
*  Petrelli A., Di Stefano L., "Pairwise registration by local orientation cues", Computer Graphics Forum, 2015
*
* \author Alioscia Petrelli
*/
namespace Pairwise3DRegistrationEvaluation
{
    /** \brief Interface that wraps the algorithm to test:
    *
    * Implement the PairwiseRegistration() method
    * \author Alioscia Petrelli
    */
	class IPairwise3DRegistrationAlgorithm
	{
	public:
		/** \brief Perform the rigid motion that aligns cloud_ref against cloud_trg
		*
		* The user must implement this method. See Pairwise3DRegistrationEvaluation_testmain.cpp for an example.
		* Indices of meshes in m_cMeshes are also provided. The user can use them to refer additional data regarding the meshes (as, for example, already computed feature points and local descriptors)
		* \param cloud_trg target mesh
		* \param cloud_ref reference mesh. The returned rigid motion should estimate the rigid motion that moves cloud_ref to cloud_trg
		* \param idx_trg index of cloud_trg in m_cMeshes
		* \param idx_ref index of cloud_ref in m_cMeshes
		* \return the rigid motion that aligns cloud_ref against cloud_trg. If the registration fails, return NULL
		*/
		virtual vtkTransform* PairwiseRegistration(vtkPolyData* cloud_trg, vtkPolyData* cloud_ref, const int idx_trg, const int idx_ref, double &secCPUtime) = 0;
	};


    /** \brief Perform the evaluation of the algorithm implementing IPairwise3DRegistrationAlgorithm.
    *
    * \author Alioscia Petrelli
    */
	class Pairwise3DRegistrationBenchmark
	{
	public:
		Pairwise3DRegistrationBenchmark();
		virtual ~Pairwise3DRegistrationBenchmark();

		/** \brief Set the dataset to test as the folder containing the views and the ground truth file
		*/
		void SetDatasetPath(const std::string &absDatasetPath ){m_absDatasetPath = absDatasetPath; };

		/** \brief Return the meshes loaded from absDatasetPath. Run PrepareEvaluation first.*/
		vtkPolyDataCollection* GetMeshes(){return m_cMeshes;};

		/** \brief Return the average mesh resolution of meshes loaded from absDatasetPath. Run PrepareEvaluation first.*/
		double GetMeanMeshRes(){return m_meanMeshRes;};
		
		/** \brief Set the algorithm to test
		*/
		void SetAlgorithm(IPairwise3DRegistrationAlgorithm* ptrAlgorithm){m_ptrAlgorithm = ptrAlgorithm;};
		
		/** \brief Set the pair of views on which the evaluation is performed.
		*
		* If not set, the evaluation is performed on each view pair.
		*/
		void SetViewPairs(std::vector< std::pair<int,int> > &vViewPairs){m_vViewPairs = vViewPairs;};

		/** \brief Set the extension of the files storing the partial views. Default: ply
		*/
		void SetMeshFileExtension(std::string &meshFileExt ){m_meshFileExt =  meshFileExt;};

		/** \brief Load partial views and ground truth.	*/
		void PrepareEvaluation();

		/** \brief Perform the evaluation for each pair of views.
		*
		* Apply the algorithm for all the view pairs and computes the figure of merits.
		*/
		void Evaluate();

		/** \brief Return the number of correctly aligned view pairs. */
		int GetNRegistrations(){return m_nRegistrations;}

		/** \brief Return the number of evaluated view pairs. */
		int GetNumViewPairs_Evaluated(){return (int)m_vViewPairs.size();}

		/** \brief Return the accuracy of the correctly aligned view pairs. */
		double GetRMSE(){return m_meanRMSE_coarse;}

		/** \brief Return the CPU time spent to align all the view pairs. Divide by GetNumViewPairs() to get the average CPU time.*/
		double GetTotalCPUtime(){return m_cpuTime_Tot;}

		/** \brief Save in a csv file the result of the evaluation
		 *
		 *  \param absCsvFileName absolute path of the csv file
		 *  \param NRegistrations number of correctly aligned view pairs
		 *  \param RMSE accuracy of the correctly aligned view pairs
		 *  \param CPUTime_Tot CPU time spent to align all the view pairs
		 */
		void PrintResults(const std::string &absCsvFileName, const int NRegistrations, const double RMSE, const double CPUTime_Tot);	

		/** \brief Compute the percentages of overlapping area between each view pair (after the ground truth rigid motion is applied) and store them in absOutFilename.
		*/
		void ComputeOverlappingAreas(const std::string &absOutFilename, const double overlappingAreaMaxDistance = 1.0, const double overlappingAreaAbsNormalDotThresh = 0.9);

	protected:
		/** \brief Path containing the views and the ground truth file.*/
		std::string								m_absDatasetPath;

		/** \brief Algorithm to test.*/
		IPairwise3DRegistrationAlgorithm*		m_ptrAlgorithm;

		/** \brief List of view pairs considered in the evaluation. If not set, all the view pairs are used.*/
		std::vector< std::pair<int,int> >		m_vViewPairs;

		/** \brief Extension of the files storing the partial views. Default: ply*/
		std::string								m_meshFileExt;


		/** \brief List of partial views comprising the dataset*/
		vtkPolyDataCollection*					m_cMeshes;

		/** \brief List of absolute path of the files storing the partial views*/
		std::vector<std::string>				m_vMeshAbsFileNames;

		/** \brief Average mesh resolution (average length of adjacent vertices) fo the partial views*/
		double									m_meanMeshRes;

		/** \brief M ground truth rigid motions that align the M partial views in a global reference frame*/
		std::vector<vtkTransform*>				m_vTransforms_GT;

		/** \brief N = M*(M-1)/2 ground truth rigid motions that aligns each pair of the M partial views*/
		std::vector<vtkTransform*>				m_vPairwiseTransforms_GT;

		/** \brief N = M*(M-1)/2 rigid motions estimated by the tested algorithm*/
		std::vector<vtkTransform*>				m_vPairwiseTransforms;

		/** \brief N = M*(M-1)/2 cpu times (secs) spent for the registration of the N = M*(M-1)/2 view pairs*/
		std::vector<double>						m_vPairwiseCpuTimes;

		/** \brief Figure of merit denoting the number of correctly aligned view pairs*/
		int										m_nRegistrations;

		/** \brief Figure of merit denoting the accuracy of the correctly aligned view pairs*/
		double									m_meanRMSE_coarse;

		/** \brief Figure of merit denoting the average cpu time spent to perform the registration of a view pair*/
		double									m_cpuTime_Tot;
		
		/** \brief Compute the percentages of overlapping area between each view pair (after the ground truth rigid motion is applied).*/
		static void computeOverlappingAreas(vtkPolyDataCollection* dataset, std::vector<vtkTransform*> &vPairwiseTransforms, std::vector<std::pair<int,int>> &vViewPairs, std::vector<double> &vOverAreas, const double maxDistance, const double absNormalDotThresh);

		/** \brief Save as filename file the percentages of overlapping area between each view pair (after the ground truth rigid motion is applied).*/
		static void printOverlappingAreas(const std::string &filename, const std::vector<std::string> &vMeshFileNames, const std::vector<std::pair<int,int>> &vMeshPairs, const std::vector<double> &vOverAreas);

	public:
		/** \brief If false, the framework loads all the partial views at once and keeps them in memory as long as the evaluation is performed. Otherwise it loads only two views at time on demand of the framework*/
		bool									m_lowMemory;

		/** \brief Threshold used to check if an estimated rigid motion is correct (after ICP is performed or not, depending on m_skipICP)*/
		double									m_RMSEthresh;

		/** \brief If true, the framework skips the ICP registration and directly checks the correctness of the rigid motion estimated by the algorithm.*/
		bool									m_skipICP;

		/** \brief ICP param: max number of iterations.*/
		int										m_ICPnMaxIters;

		/** \brief ICP param: epsilon.*/
		float									m_ICPeps;

		/** \brief ICP param: maxi radius used in neighbour search.*/
		float									m_ICPmaxRadius;

		/** \brief ICP param: uniform sampling ratio.*/
		float									m_ICPsamplingPerc;

	};

}

#endif