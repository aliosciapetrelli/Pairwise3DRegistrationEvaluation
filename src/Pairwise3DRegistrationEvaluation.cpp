#include "Pairwise3DRegistrationEvaluation.h"

#include "Pairwise3DRegistrationEvaluation_FineRegistration.h"

#include <stdexcept>

using namespace::std;

Pairwise3DRegistrationEvaluation::Pairwise3DRegistrationBenchmark::Pairwise3DRegistrationBenchmark():
m_ptrAlgorithm(NULL),
m_meshFileExt("ply"),
m_cMeshes(NULL),
m_meanMeshRes(0),
m_nRegistrations(0),
m_meanRMSE_coarse(0.0),
m_cpuTime_Tot(0.0),
m_lowMemory(false),
m_RMSEthresh(5.0f),
m_skipICP(false),
m_ICPnMaxIters(10),
m_ICPeps(0.1f),
m_ICPmaxRadius(8.0f),
m_ICPsamplingPerc(1.0f)
{

}


Pairwise3DRegistrationEvaluation::Pairwise3DRegistrationBenchmark::~Pairwise3DRegistrationBenchmark()
{
	DeleteTransforms(m_vTransforms_GT);
	DeleteTransforms(m_vPairwiseTransforms_GT);
	DeleteTransforms(m_vPairwiseTransforms);

	if(m_cMeshes)
	{
		m_cMeshes->Delete();
		m_cMeshes = NULL;
	}
}

void Pairwise3DRegistrationEvaluation::Pairwise3DRegistrationBenchmark::PrepareEvaluation()
{
	//check if the benchmark is correctly set

	if(m_absDatasetPath == "")
	{
		throw runtime_error("ERROR: absDatasetPath not set. Use SetDatasetPath() for setting the path containing partial views");
	}

	if(!ExistsDir(m_absDatasetPath))
	{
		throw runtime_error("ERROR: absDatasetPath does not exist. Use SetDatasetPath() for setting a correct path");
	}

	if(m_ptrAlgorithm == NULL)
	{
		throw runtime_error("ERROR: algorithm not set. Use SetAlgorithm() for setting the algorithm to test");
	}



	//read meshes

	if(m_cMeshes)
	{
		m_cMeshes->Delete();
		m_cMeshes = NULL;
	}
	if(m_lowMemory)
	{
		FindFilesEndWith(m_absDatasetPath, m_meshFileExt, m_vMeshAbsFileNames); 
	}
	else
	{
		m_cMeshes = ReadMeshes(m_absDatasetPath, m_vMeshAbsFileNames, m_meshFileExt, true  );
	}

	//compute mesh resolution
	m_meanMeshRes = ComputeMeshResolution(m_cMeshes, m_absDatasetPath + "/MeshRes.txt");


	//read ground truth

	string groundTruthFile = m_absDatasetPath + "/groundtruth.txt";
	if(GetRelativeName(m_absDatasetPath).find("TUM_room") == string::npos)
	{
		LoadGroundTruth(groundTruthFile, m_vMeshAbsFileNames, m_vTransforms_GT, m_cMeshes );
	}

	if(m_vViewPairs.size() == 0)
	{
		CalcAllViewPairs(m_vViewPairs, m_vMeshAbsFileNames.size());
	}

	if(GetRelativeName(m_absDatasetPath).find("TUM_room") != string::npos )
	{
		string groundTruthFile_pairwise = m_absDatasetPath + "/groundTruth_pairwise.txt";
		vector<string> vFakeNames(m_vViewPairs.size());
		for(size_t pa=0; pa<m_vViewPairs.size(); pa++)
		{
			stringstream strste;
			strste << pa << "_" << m_vViewPairs[pa].first << "_" << m_vViewPairs[pa].second << ".ply";
			vFakeNames[pa] = strste.str();
		}
		LoadGroundTruth(groundTruthFile_pairwise, vFakeNames, m_vPairwiseTransforms_GT );
	}
	else
	{
		ComputeGroundTruth_All(m_vTransforms_GT, m_vViewPairs, m_vPairwiseTransforms_GT);
	}


	m_vPairwiseTransforms.resize(m_vViewPairs.size(), NULL);
	m_vPairwiseCpuTimes.resize(m_vViewPairs.size(), 0.0);

}

void Pairwise3DRegistrationEvaluation::Pairwise3DRegistrationBenchmark::Evaluate()
{

	//perform pairwise registrations
	m_cpuTime_Tot = 0.0;
	for(size_t pa=0; pa<m_vViewPairs.size(); pa++)
	{
		int idx_trg = m_vViewPairs[pa].first;
		int idx_ref = m_vViewPairs[pa].second;

		cout << "Align pair: " << pa << " / " <<  m_vViewPairs.size() << "\ttrg: " << idx_trg << "\tref: " << idx_ref << " ... ";

		vtkPolyData* mesh_trg = NULL;
		vtkPolyData* mesh_ref = NULL;
		if(m_lowMemory)
		{
			mesh_trg = ReadMesh(m_vMeshAbsFileNames[idx_trg], true);
			mesh_ref = ReadMesh(m_vMeshAbsFileNames[idx_ref], true);
		}
		else
		{
			mesh_trg = (vtkPolyData*)(m_cMeshes->GetItemAsObject(idx_trg));
			mesh_ref = (vtkPolyData*)(m_cMeshes->GetItemAsObject(idx_ref));
		}

		m_vPairwiseTransforms[pa] = m_ptrAlgorithm->PairwiseRegistration(mesh_trg,  mesh_ref, idx_trg, idx_ref, m_vPairwiseCpuTimes[pa]);

		cout << "Pairwise CpuTime :" << m_vPairwiseCpuTimes[pa] << endl;

		m_cpuTime_Tot += m_vPairwiseCpuTimes[pa];


		if(m_lowMemory)
		{
			mesh_ref->Delete();
			mesh_trg->Delete();
		}
	}


	//Compute the 3 figure of merits


	if( !m_cMeshes ) 
	{
		if(!m_lowMemory)
		{
			m_cMeshes = ReadMeshes(m_absDatasetPath, m_vMeshAbsFileNames, m_meshFileExt, true);
		}
	}

	#ifdef GICP_IS_NOT_INCLUDED
		cout << "WARNING: GICP is not included: Use standard ICP instead. Accordingly, results turn out different from that reported in Petrelli A., Di Stefano L., \"Pairwise registration by local orientation cues\", Computer Graphics Forum, 2015" << endl;
		ICPPairwiseRegistration icp(m_ICPmaxRadius*m_meanMeshRes, m_ICPeps*m_meanMeshRes, m_ICPnMaxIters);
	#else
		GICPPairwiseRegistration gicp(m_ICPnMaxIters, m_ICPeps*m_meanMeshRes, m_ICPmaxRadius*m_meanMeshRes, false);
	#endif


	double RMSE_coarse;
	double RMSE_icp;
	m_nRegistrations = 0;
	m_meanRMSE_coarse = 0.0;
	double RMSE_ToCheck;
	vector<float> vRMSE_coarse(m_vPairwiseTransforms.size(), -1.0f);
	vector<float> vRMSE_fine(m_vPairwiseTransforms.size(), -1.0f);
	for(size_t pa=0; pa<m_vViewPairs.size(); pa++)
	{
		int idx_trg = m_vViewPairs[pa].first;
		int idx_ref = m_vViewPairs[pa].second;

		cout << "Evaluation pair: " << pa << " / " <<  m_vViewPairs.size() << "\ttrg: " << idx_trg << "\tref: " << idx_ref << endl;

		vtkPolyData* mesh_trg;
		vtkPolyData* mesh_ref;
		if(m_lowMemory)
		{
			mesh_trg = ReadMesh(m_vMeshAbsFileNames[idx_trg], true);
			mesh_ref = ReadMesh(m_vMeshAbsFileNames[idx_ref], true);
		}
		else
		{
			mesh_trg = (vtkPolyData*)(m_cMeshes->GetItemAsObject(idx_trg));
			mesh_ref = (vtkPolyData*)(m_cMeshes->GetItemAsObject(idx_ref));
		}

		if(!m_vPairwiseTransforms[pa])
		{
			if(m_lowMemory)
			{
				mesh_ref->Delete();
				mesh_trg->Delete();
			}

			continue;
		}

		vtkPolyData* mesh_ref_roto = Duplicate( mesh_ref );
		vtkPolyData* mesh_ref_GT = Duplicate( mesh_ref );
		ApplyTransform(mesh_ref_GT, m_vPairwiseTransforms_GT[pa]);
		ApplyTransform(mesh_ref_roto, m_vPairwiseTransforms[pa]);


		RMSE_coarse = CalcRMSE(mesh_ref_GT, mesh_ref_roto, m_meanMeshRes);

		vRMSE_coarse[pa] = RMSE_coarse;

		RMSE_ToCheck = RMSE_coarse;
		if(!m_skipICP)
		{

			vtkTransform* ICPtransform = NULL;

			if(m_ICPsamplingPerc < 1.0)
			{
				vtkPolyData* mesh_trg_sampl = CreatePolyData_RandomSampling(mesh_trg, m_ICPsamplingPerc);
				vtkPolyData* mesh_ref_roto_sampl = CreatePolyData_RandomSampling(mesh_ref_roto, m_ICPsamplingPerc);
				#ifdef GICP_IS_NOT_INCLUDED
					ICPtransform = icp.Register(mesh_trg_sampl, mesh_ref_roto_sampl);
				#else
					ICPtransform = gicp.Register(mesh_trg_sampl, mesh_ref_roto_sampl);
				#endif
				mesh_trg_sampl->Delete();
				mesh_ref_roto_sampl->Delete();
			}
			else
			{
				#ifdef GICP_IS_NOT_INCLUDED
					ICPtransform = icp.Register(mesh_trg, mesh_ref_roto);
				#else
					ICPtransform = gicp.Register(mesh_trg, mesh_ref_roto);
				#endif
			}

			if(!ICPtransform)
			{
				mesh_ref_roto->Delete();
				mesh_ref_GT->Delete();

				continue;		

			}

			ApplyTransform(mesh_ref_roto, ICPtransform); 



			//Check RMSE
			RMSE_icp = CalcRMSE(mesh_ref_GT, mesh_ref_roto, m_meanMeshRes);
			vRMSE_fine[pa] = RMSE_icp;

			RMSE_ToCheck = RMSE_icp;
			//cout << "RMSE_coarse: " << RMSE_coarse << "   RMSE_icp: " << RMSE_icp;
		}

		if(RMSE_ToCheck <= m_RMSEthresh)
		{
			m_nRegistrations++;
			m_meanRMSE_coarse += RMSE_coarse;
		}

		mesh_ref_roto->Delete();
		mesh_ref_GT->Delete();

		if(m_lowMemory)
		{
			mesh_ref->Delete();
			mesh_trg->Delete();
		}


	}

	if(m_meanRMSE_coarse > 0.0)
	{
		m_meanRMSE_coarse /= double(m_nRegistrations);
	}

}

void Pairwise3DRegistrationEvaluation::Pairwise3DRegistrationBenchmark::PrintResults(const std::string &absCsvFileName, const int NRegistrations, const double RMSE, const double CPUTime_Tot)
{
	std::ofstream	outFile( absCsvFileName, std::ios_base::out );

	if (!outFile.is_open() )
	{
		throw runtime_error("ERROR: Impossible to open file " + absCsvFileName);
	}	

	//print header
	outFile << "Dataset;";
	outFile << "N View Pairs (Tot);N View Pairs (Evaluated);";
	outFile << "N° Registrations;N° Registrations(%);";
	outFile << "RMSE (mr);";
	outFile << "CPU time (sec);";
	outFile << endl;

	int nViewPairs_Tot = (int)((m_vMeshAbsFileNames.size() * (m_vMeshAbsFileNames.size()-1))/2.0);
	int nViewPairs_Evaluated = GetNumViewPairs_Evaluated();

	outFile << m_absDatasetPath << ";";

	outFile << nViewPairs_Tot << ";";
	outFile << nViewPairs_Evaluated << ";";

	outFile << NRegistrations << ";";
	outFile << NRegistrations/(double)nViewPairs_Evaluated << ";";

	outFile << RMSE << ";";

	outFile << (CPUTime_Tot/double(nViewPairs_Evaluated)) << ";";
}


void Pairwise3DRegistrationEvaluation::Pairwise3DRegistrationBenchmark::ComputeOverlappingAreas(const std::string &absOutFilename, const double overlappingAreaMaxDistance, const double overlappingAreaAbsNormalDotThresh)
{
	//check if the benchmark is correctly set

	if(m_absDatasetPath == "")
	{
		throw runtime_error("ERROR: absDatasetPath not set. Use SetDatasetPath() for setting the path containing partial views");
	}

	if(!ExistsDir(m_absDatasetPath))
	{
		throw runtime_error("ERROR: absDatasetPath does not exist. Use SetDatasetPath() for setting a correct path");
	}

	//read meshes
	if(m_cMeshes)
	{
		m_cMeshes->Delete();
		m_cMeshes = NULL;
	}
	m_cMeshes = ReadMeshes(m_absDatasetPath, m_vMeshAbsFileNames, m_meshFileExt, true  );

	//compute mesh resolution
	m_meanMeshRes = ComputeMeshResolution(m_cMeshes, m_absDatasetPath + "/MeshRes.txt");

	//read ground truth
	string groundTruthFile = m_absDatasetPath + "/groundtruth.txt";
	if(GetRelativeName(m_absDatasetPath).find("TUM_room") != 0 )
	{
		LoadGroundTruth(groundTruthFile, m_vMeshAbsFileNames, m_vTransforms_GT, m_cMeshes );
	}

	if(m_vViewPairs.size() == 0)
	{
		CalcAllViewPairs(m_vViewPairs, m_cMeshes->GetNumberOfItems());
	}


	if(GetRelativeName(m_absDatasetPath).find("TUM_room") == 0 )
	{
		string groundTruthFile_pairwise = m_absDatasetPath + "/groundTruth_pairwise.txt";
		vector<string> vFakeNames(m_vViewPairs.size());
		for(size_t pa=0; pa<m_vViewPairs.size(); pa++)
		{
			stringstream strste;
			strste << pa << "_" << m_vViewPairs[pa].first << "_" << m_vViewPairs[pa].second << ".ply";
			vFakeNames[pa] = strste.str();
		}
		LoadGroundTruth(groundTruthFile_pairwise, vFakeNames, m_vPairwiseTransforms_GT );
	}
	else
	{
		ComputeGroundTruth_All(m_vTransforms_GT, m_vViewPairs, m_vPairwiseTransforms_GT);
	}


	vector<double> vOverArea;
	computeOverlappingAreas(m_cMeshes, m_vPairwiseTransforms_GT, m_vViewPairs, vOverArea, overlappingAreaMaxDistance * m_meanMeshRes, overlappingAreaAbsNormalDotThresh);
	printOverlappingAreas(absOutFilename, m_vMeshAbsFileNames, m_vViewPairs, vOverArea);
}


void Pairwise3DRegistrationEvaluation::Pairwise3DRegistrationBenchmark::computeOverlappingAreas(vtkPolyDataCollection* dataset, std::vector<vtkTransform*> &vPairwiseTransforms, std::vector<std::pair<int,int>> &vViewPairs, std::vector<double> &vOverAreas, const double maxDistance, const double absNormalDotThresh)
{
	vOverAreas.clear();

	int nViews = dataset->GetNumberOfItems();
	double overlappingArea;
	for(size_t pa=0; pa<vViewPairs.size(); pa++)
	{
		int iView1 = vViewPairs[pa].first;
		int iView2 = vViewPairs[pa].second;

		cout << iView1 << "   " << iView2 << endl;

		vtkPolyData* mesh1 = (vtkPolyData*)(dataset->GetItemAsObject(iView1));
		vtkPolyData* mesh2 = (vtkPolyData*)(dataset->GetItemAsObject(iView2));
		vtkPolyData* mesh2GT = Duplicate(mesh2);
		ApplyTransform(mesh2GT, vPairwiseTransforms[pa]);

		overlappingArea = CalcOverlappingAreaMax(mesh1, mesh2GT, maxDistance, absNormalDotThresh);

		vOverAreas.push_back(overlappingArea);

		mesh2GT->Delete();
	}

}

void Pairwise3DRegistrationEvaluation::Pairwise3DRegistrationBenchmark::printOverlappingAreas(const std::string &filename, const std::vector<std::string> &vMeshFileNames, const std::vector<std::pair<int,int>> &vMeshPairs, const std::vector<double> &vOverAreas)
{
	std::ofstream outFile(filename);

	if(!outFile.is_open())
	{
		cout << "ERROR (PrintOverlappingAreas): impossible to open file for writing " << filename;
		getchar();
		exit(-1);
	}

	outFile << "Num Meshes:" << endl;
	outFile << vMeshFileNames.size() << endl;

	outFile << "Num View pairs:" << endl;
	outFile << vMeshPairs.size() << endl;

	outFile << "view pair ID; first view ID; second view ID; first view name; second view name; overlapping area percentage" << endl;

	for(size_t pa=0; pa<vMeshPairs.size(); pa++)
	{
		outFile << pa << " " << vMeshPairs[pa].first << " " << vMeshPairs[pa].second << " "; 
		outFile << GetRelativeName_WithoutExt( vMeshFileNames[ vMeshPairs[pa].first ] ) << " ";
		outFile << GetRelativeName_WithoutExt( vMeshFileNames[ vMeshPairs[pa].second ] ) << " ";
		outFile << vOverAreas[pa] << endl;
	}

	outFile.close();
}