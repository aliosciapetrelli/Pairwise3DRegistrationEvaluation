#include "Pairwise3DRegistrationEvaluation_utils.h"


#include <algorithm>


using namespace std;


#ifdef _MSC_VER
    #ifndef _CRT_SECURE_NO_WARNINGS
    #define _CRT_SECURE_NO_WARNINGS
    #endif
    #define NOMINMAX
    #include "windows.h"
#endif

#include <io.h>   // For access().
#include <sys/types.h>  // For stat().
#include <sys/stat.h>   // For stat().

#include <cassert>


#include "vtkTransformFilter.h"
#include "vtkIdTypeArray.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkPolyDataWriter.h"
#include "vtkErrorCode.h"
#include "vtkPLYReader.h"
#include "vtkOBJReader.h"
#include "vtkSTLReader.h"
#include "vtkPolyDataReader.h"
#include "vtkSimplePointsReader.h"
#include "vtkPolyDataNormals.h"
#include "vtkCleanPolyData.h"
#include "vtkKdTreePointLocator.h"
#include "vtkPointLocator.h"
#include "vtkMath.h"
#include "vtkPointData.h"
#include "vtkTriangle.h"


bool Pairwise3DRegistrationEvaluation::ExistsDir(const string &directory)
{
	if ( _access( directory.c_str(), 0 ) == 0 )
    {
        struct stat status;
        stat( directory.c_str(), &status );

        if ( status.st_mode & S_IFDIR )
        {
            return true;
        }

		//The path you entered is a file.
		return false;
    }

	return false;
}



string Pairwise3DRegistrationEvaluation::GetRelativeName(const string &filename)
{
	std::string strFinal = filename;
	size_t posSlash = strFinal.rfind('/');
	size_t posBackSlash = strFinal.rfind('\\');
	if( (posSlash == std::string::npos) && (posBackSlash == std::string::npos) )
	{
		return strFinal;
	}

	if (posSlash == std::string::npos)
	{
		return strFinal.substr(posBackSlash+1);
	}

	if (posBackSlash == std::string::npos)
	{
		return strFinal.substr(posSlash+1);
	}

	return strFinal.substr(Max(posSlash+1, posBackSlash+1));
}

vtkPolyData* Pairwise3DRegistrationEvaluation::Duplicate(vtkPolyData* source)
{
	vtkPolyData* poly = vtkPolyData::New();
	poly->DeepCopy(source);
	return poly;
}


void Pairwise3DRegistrationEvaluation::ApplyTransform(vtkPolyData *polyData, vtkAbstractTransform* transform)
{
	vtkTransformFilter* transformFilter = vtkTransformFilter::New();
	//transformFilter->SetInputConnection(polyData->GetProducerPort());
	transformFilter->SetInput(polyData);
	transformFilter->SetTransform(transform);
	transformFilter->Update();
	polyData->ShallowCopy(transformFilter->GetPolyDataOutput());
	transformFilter->Delete();
}

double Pairwise3DRegistrationEvaluation::CalcRMSE(float* points1, float* points2, int nPoints, double meshRes)
{
	double rmse = Pairwise3DRegistrationEvaluation::calcDistance2BetweenPoints(points1, points2, nPoints);
	rmse /=(double)nPoints;
	rmse /= (meshRes*meshRes);
	rmse = sqrt(rmse);

	return rmse;
}

double Pairwise3DRegistrationEvaluation::CalcRMSE(vtkPolyData* polyData1, vtkPolyData* polyData2, double meshRes)
{
	if(polyData1->GetNumberOfPoints() != polyData2->GetNumberOfPoints())
	{
		return -1.0;
	}

	return CalcRMSE( (float*)polyData1->GetPoints()->GetVoidPointer(0), (float*)polyData2->GetPoints()->GetVoidPointer(0), polyData1->GetNumberOfPoints(), meshRes );
}


double Pairwise3DRegistrationEvaluation::calcDistance2BetweenPoints(float* points1, float* points2, int nPoints)
{
	double rmse = 0;
	float* ptrPoints1 = points1;
	float* ptrPoints2 = points2;
	for(int po=0; po<nPoints; po++)
	{
		rmse += EuclideanDistance2_3D(ptrPoints1, ptrPoints2);
		ptrPoints1 += 3;
		ptrPoints2 += 3;
	}

	return rmse;
}

vtkPolyData* Pairwise3DRegistrationEvaluation::CreatePolyData_RandomSampling(vtkPolyData* cloud, const float samplingPerc)
{

	vector<bool> vSamples(cloud->GetNumberOfPoints(), false);
	int nSampledPoints = (int)Sample_random(vSamples, samplingPerc, 1);

	vtkIdType* pointIds = new vtkIdType[nSampledPoints];
	vtkIdType* ptrpointIds = pointIds;
	for(vtkIdType po = 0; po < cloud->GetNumberOfPoints(); po++)
	{
		if(vSamples[po])
		{
			*ptrpointIds = po;
			++ptrpointIds;
		}
	}

	vtkPolyData* sampledCloud = CreatePolyData(cloud, pointIds, nSampledPoints, false, false, false);
	
	delete[] pointIds;

	return sampledCloud;
}


vtkPolyData* Pairwise3DRegistrationEvaluation::CreatePolyData(vtkPolyData* polyData, vtkIdType* pointIds, vtkIdType nPointIds, const bool insertNormals, const bool insertColors, const bool insertTriangles)
{
	vtkPolyData* poly = vtkPolyData::New();

	vtkIdType* newPointsPos = new vtkIdType[polyData->GetNumberOfPoints()];
	for(int po=0; po<polyData->GetNumberOfPoints(); po++)
	{
		newPointsPos[po] = -1;
	}

	//insert points
	AllocatePoints(poly, nPointIds);
	float* ptrNewPoints = (float*)(poly->GetPoints()->GetVoidPointer(0));
	float* ptrOldPoints = (float*)(polyData->GetPoints()->GetVoidPointer(0));
	for(int po=0; po< nPointIds; po++)
	{
		memcpy(ptrNewPoints, ptrOldPoints+ 3 * pointIds[po], sizeof(float)*3);
		ptrNewPoints+=3;
		newPointsPos[(pointIds[po])] = po;
	}

	//insert point normals
	if(insertNormals)
	{
		float* polyDataNormals = GetPolyDataPointNormalsPointer(polyData);
		if(polyDataNormals != NULL)
		{
			float* polyNormals = AllocateNormals(poly);
			for(int i=0; i<nPointIds; i++)
			{
				polyNormals[3*i] = polyDataNormals[3*pointIds[i]];
				polyNormals[3*i + 1] = polyDataNormals[3*pointIds[i] + 1];
				polyNormals[3*i + 2] = polyDataNormals[3*pointIds[i] + 2];
			}
		}
	}

	//insert point colors
	if(insertColors)
	{
		unsigned char* polyDataColors = GetPolyDataColorsPointer(polyData);
		if(polyDataColors != NULL)
		{
			unsigned char* polyColors = AllocateColors(poly);
			for(int i=0; i<nPointIds; i++)
			{
				polyColors[3*i] = polyDataColors[3*pointIds[i]];
				polyColors[3*i + 1] = polyDataColors[3*pointIds[i] + 1];
				polyColors[3*i + 2] = polyDataColors[3*pointIds[i] + 2];
			}
		}
	}


	//insert triangles
	if(insertTriangles)
	{
		if(polyData->GetNumberOfPolys() > 0)
		{
			vtkIdTypeArray* newCells = vtkIdTypeArray::New();
			newCells->SetNumberOfComponents(4);
			newCells->SetNumberOfTuples(polyData->GetNumberOfPolys());
			vtkIdType* ptrOldCells = (vtkIdType *)(polyData->GetPolys()->GetPointer());
			vtkIdType* ptrNewCells = (vtkIdType *)(newCells->GetPointer(0));
			int nNewCells = 0;
			for(int ce=0; ce<polyData->GetNumberOfPolys(); ce++)
			{
				if( (newPointsPos[ptrOldCells[1]]!=-1) && (newPointsPos[ptrOldCells[2]]!=-1) && (newPointsPos[ptrOldCells[3]]!=-1) )
				{
					*ptrNewCells = 3;
					ptrNewCells++;
					*ptrNewCells = newPointsPos[ptrOldCells[1]];
					ptrNewCells++;
					*ptrNewCells = newPointsPos[ptrOldCells[2]];
					ptrNewCells++;
					*ptrNewCells = newPointsPos[ptrOldCells[3]];
					ptrNewCells++;
					nNewCells++;
				}
				ptrOldCells += 4;
			}
			newCells->Resize(nNewCells);

			vtkCellArray* cellArray = vtkCellArray::New();
			cellArray->SetCells(nNewCells, newCells);
			AllocateTriangles(poly, nNewCells);
			poly->SetPolys(cellArray);
			cellArray->Delete();
			newCells->Delete();
		}
	}


	delete[] newPointsPos;
	return poly;
}



float* Pairwise3DRegistrationEvaluation::AllocatePoints(vtkPolyData *polyData, vtkIdType nPoints)
{
	vtkFloatArray* floatArray = vtkFloatArray::New();
	floatArray->SetNumberOfComponents(3);
	floatArray->SetNumberOfTuples(nPoints);

	vtkPoints* points3D = vtkPoints::New();
	points3D->SetData(floatArray);

	polyData->SetPoints(points3D);

	floatArray->Delete();
	points3D->Delete();

	return (float*)(polyData->GetPoints()->GetVoidPointer(0));

}

float* Pairwise3DRegistrationEvaluation::GetPolyDataPointNormalsPointer(vtkPolyData* polyData, const unsigned int iNormal)
{
	vtkDataArray* Array = polyData->GetPointData()->GetNormals();

	if(Array == NULL)
		return NULL;
	else
		return (float*)(Array->GetVoidPointer(iNormal*3));
}


float* Pairwise3DRegistrationEvaluation::AllocateNormals(vtkPolyData *polyData)
{

	vtkFloatArray* normals = vtkFloatArray::New();
	normals->SetNumberOfComponents( 3 );
	normals->SetNumberOfTuples( polyData->GetNumberOfPoints() );
	normals->SetName("Normals");
	polyData->GetPointData()->SetNormals(normals);
	normals->Delete();
	return (float *)(polyData->GetPointData()->GetNormals()->GetVoidPointer(0));
}

unsigned char* Pairwise3DRegistrationEvaluation::GetPolyDataColorsPointer(vtkPolyData* polyData, const unsigned int iColor)
{
	vtkDataArray* Array = polyData->GetPointData()->GetArray("RGB");

	if(Array == NULL)
		return NULL;
	else
		return (unsigned char*)(Array->GetVoidPointer(iColor*3));
}

unsigned char* Pairwise3DRegistrationEvaluation::AllocateColors(vtkPolyData *polyData, unsigned char* color, bool* mask)
{
	vtkUnsignedCharArray* colorArray = vtkUnsignedCharArray::New();
	colorArray->SetNumberOfComponents( 3 );
	colorArray->SetName("RGB");
	colorArray->SetNumberOfTuples( polyData->GetNumberOfPoints() );
	if(color)
	{
		unsigned char* ptrColors = (unsigned char *)(colorArray->GetVoidPointer(0));
		unsigned char noColor[3] = {200, 200, 200};
		for(int co=0; co<polyData->GetNumberOfPoints(); co++)
		{
			if(mask && !mask[co])
			{
				memcpy(ptrColors, noColor, sizeof(unsigned char)*3);
			}
			else
			{
				memcpy(ptrColors, color, sizeof(unsigned char)*3);
			}
			ptrColors+=3;
		}
	}
	polyData->GetPointData()->AddArray(colorArray);
	polyData->GetPointData()->SetActiveScalars("RGB");
	colorArray->Delete();
	return (unsigned char *)(polyData->GetPointData()->GetScalars("RGB")->GetVoidPointer(0));
}


vtkIdType* Pairwise3DRegistrationEvaluation::AllocateTriangles(vtkPolyData *polyData, const vtkIdType nTriangles)
{
	vtkIdTypeArray *idTypeTrianglesArray = vtkIdTypeArray::New();
	idTypeTrianglesArray->SetNumberOfComponents(4);
	idTypeTrianglesArray->SetNumberOfTuples(nTriangles);

	vtkCellArray *trianglesCellArray = vtkCellArray::New();
	trianglesCellArray->SetCells( nTriangles, idTypeTrianglesArray);

	polyData->SetPolys(trianglesCellArray);

	idTypeTrianglesArray->Delete();
	trianglesCellArray->Delete();

	return polyData->GetPolys()->GetPointer();

}


vtkPolyDataCollection* Pairwise3DRegistrationEvaluation::ReadMeshes(const string &meshesPath, vector<string> &vAbsMeshFileNames, const string meshExt, bool convertToVtk )
{
	FindFilesEndWith(meshesPath, meshExt, vAbsMeshFileNames, false);
	return ReadMeshes(vAbsMeshFileNames, convertToVtk);
}

vtkPolyDataCollection* Pairwise3DRegistrationEvaluation::ReadMeshes(const vector<string> &vAbsMeshFileNames, bool convertToVtk )
{
	vtkPolyDataCollection* meshes;

	vector<string> vAbsMeshFileNames_vtk;
	vAbsMeshFileNames_vtk.resize(vAbsMeshFileNames.size());
	for(size_t fi=0; fi<vAbsMeshFileNames.size(); fi++)
	{
		vAbsMeshFileNames_vtk[fi] = ChangeExtension(vAbsMeshFileNames[fi].c_str(), "vtk");
	}

	bool filesExist_ply = ExistFiles(vAbsMeshFileNames);
	bool filesExist_vtk = ExistFiles(vAbsMeshFileNames_vtk);

	if(!filesExist_ply)
	{
		cout << "No ply files to load"; 
		getchar();
		exit(-1);
	}
	if( filesExist_vtk && convertToVtk )
	{
		//load meshes in vtk format
		meshes = LoadPolyData(vAbsMeshFileNames_vtk);
		if(!meshes)
		{
			throw runtime_error("ERROR: Impossible to load vtk files"); 
		}
	}
	else
	{
		//load meshes in ply format
		cout << "Loading clouds in ply format...";
		meshes = LoadPolyData(vAbsMeshFileNames);
		if(!meshes)
		{
			throw runtime_error("ERROR: Impossible to load ply files"); 
		}

		Preprocessing(meshes);

		if( convertToVtk )
		{
			//save meshes in vtk format
			WriteVtk(vAbsMeshFileNames_vtk, meshes); 
		}
	}

	return meshes;
}


bool Pairwise3DRegistrationEvaluation::FindFilesEndWith(const std::string & path, const std::string & endingString, std::vector<std::string> & foundFiles, bool getOnlyFileName, const int nMaxItems )
{
	foundFiles.clear();

#ifdef _MSC_VER
	// Required structs for searching for files and directories
	WIN32_FIND_DATA FindFileData;
	HANDLE hFind = INVALID_HANDLE_VALUE;

	// Build the file search string...
	char searchDir[2048] = {0};
	char fullpath[2048] = {0};

	// ...if it already is a path that ends with \ or /, add '*'...
	if(path.at(path.length() - 1) == '\\' || path.at(path.length() - 1) == '/')
	{
		_snprintf(searchDir, 2047, "%s*", path.c_str());
		_snprintf(fullpath, 2047, "%s", path.c_str()); // just copy path
	}
	// ...otherwise, add '\*' to the end of the path.
	else
	{
		_snprintf(searchDir, 2047, "%s/*", path.c_str());
		_snprintf(fullpath, 2047, "%s/", path.c_str()); // copy path and add slash (required when building filenames)
	}

	// Find the first file in the directory.
	hFind = FindFirstFile(searchDir, &FindFileData);

	// If there is no file, return
	if (hFind == INVALID_HANDLE_VALUE)
	{
		return false;
	}

	// loop
	do
	{
		// Skip ".", ".." and all directories
		if( ((FindFileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) != 0) || strcmp(FindFileData.cFileName, ".") == 0 || strcmp(FindFileData.cFileName, "..") == 0)
		{
			continue;
		}

		if(endingString.size() > std::string(FindFileData.cFileName).size())
		{
			continue;
		}

		// Store filename into the vector if it starts with the given string
		if(endingString.size() > 0)
		{
			if(std::string(FindFileData.cFileName).rfind(endingString) == (std::string(FindFileData.cFileName).size() - endingString.size()) )
			{
				if(getOnlyFileName)
				{
					foundFiles.push_back(FindFileData.cFileName);
				}
				else
				{
					// File found: create a path to the file
					char filePath[2048] = {0};
					_snprintf(filePath, 2047, "%s%s", fullpath, FindFileData.cFileName);
					// Add it to vector of files found
					foundFiles.push_back(filePath);
				}
			}
		}
		else // Always store filename if a starting string has not been provided
		{
			if(getOnlyFileName)
			{
				foundFiles.push_back(FindFileData.cFileName);
			}
			else
			{
				// File found: create a path to the file
				char filePath[2048] = {0};
				_snprintf(filePath, 2047, "%s%s", fullpath, FindFileData.cFileName);
				// Add it to vector of files found
				foundFiles.push_back(filePath);
			}
		}
	}
	// Loop while we find more files
	while(FindNextFile(hFind, &FindFileData) != 0);

	// Release
	FindClose(hFind);

	sort(foundFiles.begin(), foundFiles.end(), StringCompare_Smart_Incr());

	foundFiles.resize( Min(nMaxItems, (int)foundFiles.size() ) );

	return true;
#else // not _MSC_VER
    DIR* directory = opendir(path.c_str());
    if(directory)
    {
        string parent(path);
        if(parent[parent.length()-1] != '/')
            parent.append("/");

        struct dirent dirEntry;
        struct dirent* res = &dirEntry;
        while((readdir_r(directory, &dirEntry, &res) == 0) && (res)) // thread-safe
            if((dirEntry.d_type == DT_REG) &&
                    (strncmp(dirEntry.d_name+(d_namlen-endingString.size()), endingString.c_str(), endingString.length()) == 0) &&
                    (strcmp(dirEntry.d_name, ".") != 0) &&
                    (strcmp(dirEntry.d_name, "..") != 0))
			if(getOnlyFileName)
			{
				foundFiles.push_back(dirEntry.d_name);
			}
			else
			{
                foundFiles.push_back(parent + dirEntry.d_name);
			}
        closedir(directory);

		sort(foundFiles.begin(), foundFiles.end(), StringCompare_Smart_Incr());

		foundFiles.resize( Min(nMaxItems, (int)foundFiles.size() ) );

        return true;
    }

    return false;
#endif // _MSC_VER
}

std::string Pairwise3DRegistrationEvaluation::ChangeExtension(const std::string &filename, const std::string &newExt)
{
	string strTemp = filename;
	size_t dotPos = strTemp.find_last_of(".");
	if(dotPos == std::string::npos)
	{
		strTemp = strTemp + "." + newExt;
	}
	else
	{
		strTemp.replace(dotPos+1, strTemp.size(), newExt);
	}
	return strTemp;
}


bool Pairwise3DRegistrationEvaluation::ExistFiles(const vector<string> &absFiles)
{
	for(size_t fi=0; fi<absFiles.size(); fi++)
	{
		if(!ExistsFile(absFiles[fi]))
		{
			return false;
		}
	}

	return true;
}


bool Pairwise3DRegistrationEvaluation::ExistsFile(const string &filename_abs)
{
	#ifdef _MSC_VER
		return !(GetFileAttributes(filename_abs.c_str()) == INVALID_FILE_ATTRIBUTES); 			
	#else	
		#ifndef BOOST_IS_NOT_INCLUDED
			return boost::filesystem::exists( filename_abs.c_str() );
		#else		
			#ifdef __GNUC__
				//assume that on Linux the compiler is gcc 
				struct stat st;
				return (stat(filename_abs.c_str(), &st) == 0);
			#else			
				cout << "ERROR (FileExists): unknown environment";
				getchar();
				exit(-1);

				return true;

			#endif

		#endif

	#endif
}

bool Pairwise3DRegistrationEvaluation::WriteVtk(const vector<string> &absFiles, vtkPolyDataCollection* polyDataCollection, const bool fileTypeASCII)
{
	if(absFiles.size() != polyDataCollection->GetNumberOfItems())
	{
		return false;
	}
	vtkPolyData* poly;
	polyDataCollection->InitTraversal();
	for(size_t fi=0; fi<absFiles.size(); fi++)
	{
		poly = polyDataCollection->GetNextItem();
		WriteVtk(absFiles[fi].c_str(), poly);
	}
	return true;
}


bool Pairwise3DRegistrationEvaluation::WriteVtk(const std::string &absFileName, vtkPolyData* polyData, const bool fileTypeASCII )
{
	if(polyData == NULL)
	{
		return false;
	}

	vtkPolyDataWriter* polydataWriter = vtkPolyDataWriter::New();
	polydataWriter->SetInputConnection(polyData->GetProducerPort());
	if(fileTypeASCII)
	{
		polydataWriter->SetFileTypeToASCII();
	}
	else
	{
		polydataWriter->SetFileTypeToBinary();
	}
	polydataWriter->SetFileName(absFileName.c_str());
	polydataWriter->Update();
	int errCod = polydataWriter->GetErrorCode();

	polydataWriter->Delete();

	if(errCod != vtkErrorCode::NoError)
	{
		return false;
	}
	else
	{
		return true;
	}
}


vtkPolyDataCollection* Pairwise3DRegistrationEvaluation::LoadPolyData(const std::vector<std::string> &filenames)
{
	vtkPolyDataCollection* polyDataCollection = vtkPolyDataCollection::New();
	vtkPolyData* poly;
	for(std::vector<std::string>::const_iterator iter = filenames.begin(); iter!= filenames.end(); ++iter)
	{
		poly = LoadPolyData(iter->data());
		if(poly == NULL)
		{
			polyDataCollection->Delete();
			return NULL;
		}
		polyDataCollection->AddItem(poly);
		poly->Delete();
	}

	return polyDataCollection;
}



vtkPolyData* Pairwise3DRegistrationEvaluation::LoadPolyData(const std::string &absFilename)
{
	std::string strFileName = absFilename;
	#ifndef BOOST_IS_NOT_INCLUDED
		boost::filesystem::path filepath(filename);
    #if defined (BOOST_FILESYSTEM_VERSION) && BOOST_FILESYSTEM_VERSION >= 3
      std::string ext = filepath.extension().string();
    #else
      std::string ext = filepath.extension();
    #endif
	#else
		std::string ext = strFileName.substr(strFileName.find_last_of("."));
	#endif
	std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

	if (ext == ".ply")
		return LoadPly(strFileName.c_str());
	else if (ext == ".off")
		return LoadOff(strFileName.c_str());
	else if (ext == ".obj")
		return LoadOBJ(strFileName.c_str());
	else if (ext == ".stl")
		return LoadSTL(strFileName.c_str());
	else if (ext == ".vtk")
		return LoadVtk(strFileName.c_str());
	else if (ext == ".vert")
		return LoadVertTri(strFileName.substr(0, strFileName.length() - ext.length()).c_str());
	else if (ext == ".tri")
		return LoadVertTri(strFileName.substr(0, strFileName.length() - ext.length()).c_str());
	else if (ext == ".pcd")
		return LoadPcd(strFileName.c_str());
	else if (ext == ".txt")
		return LoadSimpleTxt(strFileName.c_str());
	else if (ext == ".xyz")
		return LoadXYZ(strFileName.c_str());
	else
		throw runtime_error("ERROR Type not supported");

	return NULL;
}



vtkPolyData* Pairwise3DRegistrationEvaluation::LoadPly(const char* filename)
{
	vtkPLYReader  *plyReader = vtkPLYReader::New();
	plyReader->SetFileName(filename);
	plyReader->Update();
	int errCod = plyReader->GetErrorCode();
	if(errCod != vtkErrorCode::NoError)
	{
		const char* errString = vtkErrorCode::GetStringFromErrorCode(errCod);
		printf("\nERROR LoadPly: %s\n", errString);
		plyReader->Delete();
		return NULL;
	}

	vtkPolyData* polyData = vtkPolyData::New();
	polyData->ShallowCopy(plyReader->GetOutput());
	plyReader->Delete();
	return polyData;
}

vtkPolyData* Pairwise3DRegistrationEvaluation::LoadPcd(const char* filename)
{

#ifndef BOOST_IS_NOT_INCLUDED

	std::ifstream in(filename);
	if(!in.is_open())
	{
		perror("Problems opening Pcd file");
		return NULL;
	}



	bool rgb = false;
	int width;
	int height;
	bool createTriangles = true;
	int nPoints;


	std::string line;

	int specified_channel_count = 0;


	// field_sizes represents the size of one element in a field (e.g., float = 4, char = 1)
	// field_counts represents the number of elements in a field (e.g., x = 1, normal_x = 1, fpfh = 33)
	std::vector<int> field_sizes, field_counts;
	// field_types represents the type of data in a field (e.g., F = float, U = unsigned)
	std::vector<char> field_types;
	std::vector<std::string> st;

	// Read the header and fill it in with wonderful values
	try
	{
		while (!in.eof ())
		{
			getline (in, line);
			// Ignore empty lines
			if (line == "")
				continue;

			// Tokenize the line
			boost::trim (line);
			boost::split (st, line, boost::is_any_of ("\t\r "), boost::token_compress_on);

			std::string line_type = st.at (0);

			// Ignore comments
			if (line_type.substr (0, 1) == "#")
				continue;

			// Version numbers are not needed for now, but we are checking to see if they're there
			if (line_type.substr (0, 7) == "VERSION")
				continue;

			// Get the field indices (check for COLUMNS too for backwards compatibility)
			if ( (line_type.substr (0, 6) == "FIELDS") || (line_type.substr (0, 7) == "COLUMNS") )
			{
				if (std::find(st.begin(), st.end(), "rgb") != st.end())
					rgb = true;
				//specified_channel_count = st.size () - 1;

				//// Allocate enough memory to accommodate all fields
				//cloud.fields.resize (specified_channel_count);
				//for (int i = 0; i < specified_channel_count; ++i)
				//{
				//  std::string col_type = st.at (i + 1);
				//  cloud.fields[i].name = col_type;
				//}

				//// Default the sizes and the types of each field to float32 to avoid crashes while using older PCD files
				//int offset = 0;
				//for (int i = 0; i < specified_channel_count; ++i, offset += 4)
				//{
				//  cloud.fields[i].offset   = offset;
				//  cloud.fields[i].datatype = sensor_msgs::PointField::FLOAT32;
				//  cloud.fields[i].count    = 1;
				//}
				//cloud.point_step = offset;
				continue;
			}

			// Get the field sizes
			if (line_type.substr (0, 4) == "SIZE")
			{
				//specified_channel_count = st.size () - 1;

				//// Allocate enough memory to accommodate all fields
				//if (specified_channel_count != (int)cloud.fields.size ())
				//  throw "The number of elements in <SIZE> differs than the number of elements in <FIELDS>!";

				//// Resize to accommodate the number of values
				//field_sizes.resize (specified_channel_count);

				//int offset = 0;
				//for (int i = 0; i < specified_channel_count; ++i)
				//{
				//  int col_type = atoi (st.at (i + 1).c_str ());
				//  cloud.fields[i].offset = offset;                // estimate and save the data offsets
				//  offset += col_type;
				//  field_sizes[i] = col_type;                      // save a temporary copy
				//}
				//cloud.point_step = offset;
				////if (cloud.width != 0)
				//  //cloud.row_step   = cloud.point_step * cloud.width;
				continue;
			}

			// Get the field types
			if (line_type.substr (0, 4) == "TYPE")
			{
				/*if (field_sizes.empty ())
				throw "TYPE of FIELDS specified before SIZE in header!";*/

				//specified_channel_count = st.size () - 1;

				//// Allocate enough memory to accommodate all fields
				//if (specified_channel_count != (int)cloud.fields.size ())
				//  throw "The number of elements in <TYPE> differs than the number of elements in <FIELDS>!";

				//// Resize to accommodate the number of values
				//field_types.resize (specified_channel_count);

				//for (int i = 0; i < specified_channel_count; ++i)
				//{
				//  field_types[i] = st.at (i + 1).c_str ()[0];
				//  cloud.fields[i].datatype = getFieldType (field_sizes[i], field_types[i]);
				//}
				continue;
			}

			// Get the field counts
			if (line_type.substr (0, 5) == "COUNT")
			{
				/*if (field_sizes.empty () || field_types.empty ())
				throw "COUNT of FIELDS specified before SIZE or TYPE in header!";

				specified_channel_count = st.size () - 1;*/

				//// Allocate enough memory to accommodate all fields
				//if (specified_channel_count != (int)cloud.fields.size ())
				//  throw "The number of elements in <COUNT> differs than the number of elements in <FIELDS>!";

				//field_counts.resize (specified_channel_count);

				//int offset = 0;
				//for (int i = 0; i < specified_channel_count; ++i)
				//{
				//  cloud.fields[i].offset = offset;
				//  int col_count = atoi (st.at (i + 1).c_str ());
				//  cloud.fields[i].count = col_count;
				//  offset += col_count * field_sizes[i];
				//}
				//// Adjust the offset for count (number of elements)
				//cloud.point_step = offset;
				continue;
			}

			// Get the width of the data (organized point cloud dataset)
			if (line_type.substr (0, 5) == "WIDTH")
			{
				width = atoi (st.at (1).c_str ());
				//cloud.width = atoi (st.at (1).c_str ());
				//if (cloud.point_step != 0)
				//  cloud.row_step = cloud.point_step * cloud.width;      // row_step only makes sense for organized datasets
				continue;
			}

			// Get the height of the data (organized point cloud dataset)
			if (line_type.substr (0, 6) == "HEIGHT")
			{
				height = atoi (st.at (1).c_str ());
				continue;
			}

			// Get the acquisition viewpoint
			if (line_type.substr (0, 9) == "VIEWPOINT")
			{
				/* pcd_version = PCD_V7;
				if (st.size () < 8)
				throw "Not enough number of elements in <VIEWPOINT>! Need 7 values (tx ty tz qw qx qy qz).";

				origin      = Eigen::Vector4f (atof (st.at (1).c_str ()), atof (st.at (2).c_str ()), atof (st.at (3).c_str ()), 0);
				orientation = Eigen::Quaternionf (atof (st.at (4).c_str ()), atof (st.at (5).c_str ()), atof (st.at (6).c_str ()), atof (st.at (7).c_str ()));*/
				continue;
			}

			// Get the number of points
			if (line_type.substr (0, 6) == "POINTS")
			{
				nPoints = atoi (st.at (1).c_str ());
				// Need to allocate: N * point_step
				//cloud.data.resize (nr_points * cloud.point_step);
				continue;
			}

			// Read the header + comments line by line until we get to <DATA>
			if (line_type.substr (0, 4) == "DATA")
			{
				//data_idx = fs.tellg ();
				assert (st.at (1).substr (0, 6) == "ascii");
				break;
			}
		}
	}
	catch (const char *)
	{
		return NULL;
	}


	if( (height == 1) || (width == 1) )
	{
		createTriangles = false;
	}



	vtkPolyData* poly = vtkPolyData::New();
	//read points
	AllocatePoints(poly, nPoints);
	float* ptrPoints = GetPolyDataPointsPointer(poly);


	unsigned char* ptrRgb = NULL;
	float color;
	if(rgb)
	{
		ptrRgb = AllocateColors(poly);
	}

	//find valid point indexes
	int n = 0;
	int* pIndex = new int[width*height];	//map with indexes of valid points
	for(int j=0; j<height; j++)
	{
		for(int i=0; i<width; i++)
		{
			in >> line;
			if(line == "nan")
			{
				std::getline(in, line);
				pIndex[j*width+i] = -1;
				continue;
			}

			*ptrPoints++ = atof(line.c_str());
			in >> *ptrPoints++;
			in >> *ptrPoints++;

			if(rgb)
			{
				//conversion from float to int32 and bgr
				in >> color;
				//if(color == 0)
				//{
				//	pIndex[j*width+i] = -1;
				//	ptrPoints -= 3;
				//	continue;
				//}
				int rgb = *reinterpret_cast<int*>(&color);
				*ptrRgb++ = ((rgb >> 16) & 0xff);
				*ptrRgb++ = ((rgb >> 8) & 0xff);
				*ptrRgb++ = (rgb & 0xff);


				//memcpy(ptrRgb++, (unsigned char*)(&color) + 2 , sizeof(unsigned char));
				//memcpy(ptrRgb++, (unsigned char*)(&color) + 1 , sizeof(unsigned char));
				//memcpy(ptrRgb++, (unsigned char*)(&color) + 0 , sizeof(unsigned char));
			}

			pIndex[j*width+i] = n;
			n++;
		}
	}


	if(createTriangles)
	{
		vtkIdType* triangles;
		int nTriangles = CVPipes::Core::Impl::CreateTrianglesFromIndexMap(pIndex, width, height, triangles);

		vtkIdType* ptrPolyTriangles = AllocateTriangles(poly, nTriangles);
		memcpy(ptrPolyTriangles, triangles, sizeof(vtkIdType)*4*nTriangles);

		delete[] triangles;
	}

	delete[] pIndex;

	return poly;
#else

	throw runtime_error("ERROR cannot use LoadPcd without boost");

	return NULL;
#endif

}


vtkPolyData* Pairwise3DRegistrationEvaluation::LoadSTL(const char* filename)
{
	vtkSTLReader *STLReader = vtkSTLReader::New();
	STLReader->SetFileName(filename);
	STLReader->Update();
	int errCod = STLReader->GetErrorCode();
	if(errCod != vtkErrorCode::NoError)
	{
		const char* errString = vtkErrorCode::GetStringFromErrorCode(errCod);
		printf("\nERROR LoadSTL: %s\n", errString);
		STLReader->Delete();
		return NULL;
	}

	vtkPolyData* polyData = vtkPolyData::New();
	polyData->ShallowCopy(STLReader->GetOutput());
	STLReader->Delete();
	return polyData;
}


vtkPolyData* Pairwise3DRegistrationEvaluation::LoadOBJ(const char* filename)
{
	vtkOBJReader *OBJReader = vtkOBJReader::New();
	OBJReader->SetFileName(filename);
	OBJReader->Update();
	int errCod = OBJReader->GetErrorCode();
	if(errCod != vtkErrorCode::NoError)
	{
		printf("ERROR loading OBJ %s\nerror Code: %d, Description: %s\n", OBJReader->GetFileName(), errCod, vtkErrorCode::GetStringFromErrorCode(errCod));
		OBJReader->Delete();
		return NULL;
	}

	vtkPolyData* polyData = vtkPolyData::New();
	polyData->ShallowCopy(OBJReader->GetOutput());
	OBJReader->Delete();
	return polyData;
}

vtkPolyData* Pairwise3DRegistrationEvaluation::LoadOff(const char* filename)
{
	std::ifstream in(filename);
	if(!in.is_open())
	{
		perror("Problems opening OFF file");
		return NULL;
	}

	int n_pts, n_polys, nAux;

	std::string line;
	std::getline(in, line);
	assert(line == "OFF" || line == "COFF");


	/*char lineStart;
	fpos_t position;

	fgetpos (in, &position);
	fscanf(in, "%c", &lineStart);

  */
	std::getline(in, line);
	while (line.find_first_of('#') == 0 || line == "" )
	{
		/*while (lineStart != '\n')
			fscanf(in, "%c", &lineStart);

		fgetpos (in, &position);
		fscanf(in, "%c", &lineStart);*/
		std::getline(in, line);
	}


	//char n_pts_str[10];
	//char n_pts_start[2] = {lineStart, '\0'};

	std::stringstream lineStream(line);

	lineStream >> n_pts;
	lineStream >> n_polys;
	lineStream >> nAux;
	//fscanf(in, "%s %d %d\n", n_pts_str, &n_polys, &nAux);

	/*std::string n_pts_std_str = n_pts_start;
	n_pts_std_str += n_pts_str;
	n_pts = atoi(n_pts_std_str.c_str());*/

	vtkPolyData *poly = vtkPolyData::New();
	vtkPoints *points = vtkPoints::New();
	//points->SetDataTypeToDouble();
	vtkCellArray *polys = vtkCellArray::New();

	// Read the points from the file
	double p[3];
	int i,j;
	for(i = 0; i < n_pts; i++)
	{
		//std::getline(in, line);
		//while(line.find_first_of('#') == 0)
		//	std::getline(in, line);

		char c = '#';
		do
		{
			c = in.get();
			in.putback(c);
			if (c == '#' || c == ' ')
			{
				std::getline(in, line);
			}
		} while (c == '#' || c == ' ');

		/*lineStream.str(line);
		lineStream.clear();*/

		in >> p[0];
		in >> p[1];
		in >> p[2];
		//fscanf(in, "%f %f %f\n", &p[0], &p[1], &p[2]);
		points->InsertPoint(i,p);

		// skip comments at the end of the line
		std::getline(in, line);
	}

	assert(points->GetNumberOfPoints() == n_pts);

	// Read the triangles from the file
	for(i = 0; i < n_polys; i++)
	{
		std::getline(in, line);
		while(line.find_first_of('#') == 0 || line == "")
			std::getline(in, line);

		lineStream.str(line);
		lineStream.clear();

		vtkIdType n_vertices;
		lineStream >> n_vertices;
		//fscanf(in, "%d ", &n_vertices);
		polys->InsertNextCell(n_vertices);
		vtkIdType vertex_id;
		for(j = 0; j < n_vertices; j++)
		{

			lineStream >> vertex_id;
		   //fscanf(in, "%d ", &vertex_id);
			polys->InsertCellPoint(vertex_id);
		}

		//fscanf(in, "%d ", &n_vertices);
	}

	assert(polys->GetNumberOfCells() == n_polys);

	poly->SetPoints(points);
	points->Delete();
	poly->SetPolys(polys);
	polys->Delete();
	//fclose(in);
	return poly;
}


vtkPolyData* Pairwise3DRegistrationEvaluation::LoadVtk(const char* filename)
{
	vtkPolyDataReader *polydataReader = vtkPolyDataReader::New();
	polydataReader->SetFileName(filename);
	polydataReader->Update();
	int errCod = polydataReader->GetErrorCode();
	if(errCod != vtkErrorCode::NoError)
	{
		const char* errString = vtkErrorCode::GetStringFromErrorCode(errCod);
		printf("\nERROR LoadVtk: %s\n", errString);
		polydataReader->Delete();
		return NULL;
	}

	vtkPolyData* polyData = vtkPolyData::New();
	polyData->DeepCopy(polydataReader->GetOutput());
	polydataReader->Delete();
	return polyData;
}


vtkPolyData* Pairwise3DRegistrationEvaluation::LoadVertTri(const char* filenameWithoutExtension)
{
	std::string filenameVert(filenameWithoutExtension);
	std::string filenameTri(filenameWithoutExtension);
	filenameVert += ".vert";
	filenameTri += ".tri";

	std::ifstream inVert, inTri;
	inVert.open(filenameVert.c_str(), fstream::in);
	inTri.open(filenameTri.c_str(), fstream::in);
	if (!inVert.is_open())
	{
		cerr << "file " << filenameVert << " could not be opened" << endl;
		return NULL;
	}
	if (!inTri.is_open())
	{
		cerr << "file " << filenameTri << " could not be opened" << endl;
		return NULL;
	}

	vtkPolyData *poly = vtkPolyData::New();
	vtkPoints *points = vtkPoints::New();
	vtkCellArray *polys = vtkCellArray::New();

	double p[3], t[3];
	int i=0;

	while (!inVert.eof())
	{
		inVert >> p[0] >> p[1] >> p[2];
		points->InsertPoint(i,p);
		i++;
	}

	while (!inTri.eof())
	{
		inTri >> t[0] >> t[1] >> t[2];
		polys->InsertNextCell(3);
		polys->InsertCellPoint(t[0] - 1);
		polys->InsertCellPoint(t[1] - 1);
		polys->InsertCellPoint(t[2] - 1);
	}

	poly->SetPoints(points);
	points->Delete();
	poly->SetPolys(polys);
	polys->Delete();

	inVert.close();
	inTri.close();
	return poly;
}


vtkPolyData* Pairwise3DRegistrationEvaluation::LoadSimpleTxt(const char* filename)
{
	vtkSimplePointsReader *simplePointsReader = vtkSimplePointsReader::New();
	simplePointsReader->SetFileName(filename);
	simplePointsReader->Update();
	int errCod = simplePointsReader->GetErrorCode();
	if(errCod != vtkErrorCode::NoError)
	{
		simplePointsReader->Delete();
		return NULL;
	}

	vtkPolyData* polyData = vtkPolyData::New();
	polyData->ShallowCopy(simplePointsReader->GetOutput());
	simplePointsReader->Delete();
	return polyData;
}



vtkPolyData* Pairwise3DRegistrationEvaluation::LoadXYZ(const char* fileName)
{
	std::ifstream file(fileName);

	if (!file.is_open() )
	{
		printf("\nERROR: Wrong XYZ fileName\n");
		return NULL;
	}

	std::vector<float> vPoints;
	std::vector<float> vNormals;

	float floatVal;
	int nPoints=0;
	while(!file.eof())
	{
		file >> floatVal;
		vPoints.push_back(floatVal);
		if(file.eof())
		{
			break;
		}
		file >> floatVal;
		vPoints.push_back(floatVal);
		if(file.eof())
		{
			break;
		}
		file >> floatVal;
		vPoints.push_back(floatVal);
		if(file.eof())
		{
			break;
		}

		file >> floatVal;
		vNormals.push_back(floatVal);
		if(file.eof())
		{
			break;
		}
		file >> floatVal;
		vNormals.push_back(floatVal);
		if(file.eof())
		{
			break;
		}
		file >> floatVal;
		vNormals.push_back(floatVal);
		nPoints++;
	}
	file.close();

	vtkPolyData* poly = NULL;

	if(nPoints)
	{
		poly = CreatePolyData(&vPoints[0], nPoints);
		float* ptrNormals = AllocateNormals(poly);
		memcpy(ptrNormals, &vNormals[0], sizeof(float)*3*nPoints);
	}

	return poly;
}


vtkPolyData* Pairwise3DRegistrationEvaluation::CreatePolyData(const float* points, const vtkIdType nPoints, const unsigned char* colors, const vtkIdType* triangles, const vtkIdType nTriangles)
{
	vtkPolyData* polyData = vtkPolyData::New();

	//insert points
	AllocatePoints(polyData, nPoints);
	memcpy(polyData->GetPoints()->GetVoidPointer(0), points, nPoints*3 * sizeof(float) );

	if(colors)
	{
		//insert colors
		AllocateColors(polyData);
		memcpy(polyData->GetPointData()->GetArray("RGB")->GetVoidPointer(0), colors, nPoints* 3 * sizeof(unsigned char) );
	}

	if(triangles)
	{
		//insert triangles
		AllocateTriangles(polyData, nTriangles);

		vtkIdType* ptrTrianglesArray = (vtkIdType *)polyData->GetPolys()->GetPointer();
		const vtkIdType* ptrTriangles = triangles;
		for(vtkIdType tr=0; tr<nTriangles; tr++)
		{
			*(ptrTrianglesArray++)= 3;
			memcpy(ptrTrianglesArray, ptrTriangles, sizeof(vtkIdType)*3);
			ptrTriangles+=3;
			ptrTrianglesArray+=3;
		}
	}
	return polyData;
}



void Pairwise3DRegistrationEvaluation::Preprocessing(vtkPolyDataCollection* &polyDataCollection)
{
	cout << "Computing normals...";
	CalcMeshNormals(polyDataCollection);
	CleanPolyData(polyDataCollection, true, 0.0);
	cout << "OK" << endl;
}


void Pairwise3DRegistrationEvaluation::CalcMeshNormals(vtkPolyDataCollection *polyDataCollection, const double dFeatureAngle)
{
	vtkPolyData* iter=NULL;
	polyDataCollection->InitTraversal();
	while((iter = polyDataCollection->GetNextItem()))
	{
		CalcMeshNormals(iter, dFeatureAngle);
	}
}

void Pairwise3DRegistrationEvaluation::CalcMeshNormals(vtkPolyData *polyData, const double dFeatureAngle)
{
	if(polyData == NULL)
	{
		return;
	}
	vtkPolyDataNormals* polydataNormals = vtkPolyDataNormals::New();
	if(dFeatureAngle < 0)
	{
		polydataNormals->SplittingOff();
	}
	else
	{
		polydataNormals->SplittingOn();
		polydataNormals->SetFeatureAngle(dFeatureAngle);
	}
	//vtkAlgorithmOutput* prova = polyData->GetProducerPort();
	//polydataNormals->SetInputConnection(prova);
	polydataNormals->SetInput(polyData);
	polydataNormals->ComputePointNormalsOn();
	polydataNormals->AutoOrientNormalsOff();
	polydataNormals->ConsistencyOff();
	polydataNormals->ComputeCellNormalsOff();

	polydataNormals->Update();

	polyData->ShallowCopy( polydataNormals->GetOutput() );
	polydataNormals->Delete();
}


void Pairwise3DRegistrationEvaluation::CleanPolyData(vtkPolyDataCollection *polyDataCollection, const bool mergePoints, const double tolerance, bool removeNotPolysCells)
{
	vtkPolyData* iter=NULL;
	polyDataCollection->InitTraversal();
	while((iter = polyDataCollection->GetNextItem()))
	{
		CleanPolyData(iter, mergePoints, tolerance, removeNotPolysCells);
	}
}

void Pairwise3DRegistrationEvaluation::CleanPolyData(vtkPolyData *polyData, const bool mergePoints, const double tolerance, bool removeNotPolysCells)
{
	if(polyData == NULL)
	{
		return;
	}

	if(removeNotPolysCells)
	{
		RemoveNotPolysCells(polyData);
	}

	vtkCleanPolyData* cleanPolyData = vtkCleanPolyData::New();
	cleanPolyData->SetInputConnection(polyData->GetProducerPort());
	cleanPolyData->SetAbsoluteTolerance(tolerance);
	cleanPolyData->SetPointMerging(mergePoints);
	cleanPolyData->ToleranceIsAbsoluteOn();
	//cleanPolyData->ConvertLinesToPointsOff();
	//cleanPolyData->ConvertPolysToLinesOff();
	//cleanPolyData->ConvertStripsToPolysOff();
	//cleanPolyData->PieceInvariantOff();

	cleanPolyData->Update();
	polyData->ShallowCopy( cleanPolyData->GetOutput() );
	cleanPolyData->Delete();

	if(removeNotPolysCells)
	{
		RemoveNotPolysCells(polyData);
	}
}

void Pairwise3DRegistrationEvaluation::RemoveNotPolysCells(vtkPolyData *polyData)
{
	polyData->GetLines()->Reset();
	polyData->GetStrips()->Reset();
	polyData->GetVerts()->Reset();
}


double Pairwise3DRegistrationEvaluation::ComputeMeshResolution(vtkPolyDataCollection* polyDataCollection, const std::string &absMeshResFile)
{
	double meshRes = 0;

	bool fileExists = ExistsFile(absMeshResFile);

	if( fileExists  && (absMeshResFile != "") )
	{
		LoadSingleValue<double>(absMeshResFile, meshRes);
		return meshRes;
	}

	vtkPolyData* poly;
	polyDataCollection->InitTraversal();
	int n=0;
	while( (poly = polyDataCollection->GetNextItem()) )
	{
		meshRes += ComputeMeshResolution(poly);
		n++;
	}
	meshRes /= (double)polyDataCollection->GetNumberOfItems();
	
	if(absMeshResFile != "")
	{
		SaveSingleValue<double>(absMeshResFile, meshRes);
	}

	return meshRes;
}



double Pairwise3DRegistrationEvaluation::ComputeMeshResolution(vtkPolyData* cloud)
{
	double meshResolution = 0;
	int numdistances = 0;

	if (cloud == NULL)
		return 0.0;

	vtkCellArray* polys = cloud->GetPolys();

	if (polys == NULL)
		return 0.0;

	vtkIdType* ptrCellIds = polys->GetPointer();

	if (ptrCellIds == NULL)
		return 0.0;

	int nEdges;
	int firstVer;
	int secondVer;
	double firstPoint[3], secondPoint[3];
	for(int ce=0; ce<cloud->GetNumberOfPolys(); ce++)
	{
		nEdges = *ptrCellIds;
		ptrCellIds++;
		for(int ed=0; ed<nEdges; ed++)
		{
			firstVer = *(ptrCellIds+ed);
			secondVer = *(ptrCellIds+((ed+1)%nEdges));
			cloud->GetPoint(firstVer, firstPoint);
			cloud->GetPoint(secondVer, secondPoint);
			meshResolution += EuclideanDistance_3D(firstPoint, secondPoint );
			numdistances++;
		}
		ptrCellIds+=nEdges;
	}
	if(numdistances!=0)
	{
		meshResolution/=numdistances;
	}

	return meshResolution;
}



void Pairwise3DRegistrationEvaluation::LoadGroundTruth(const string &strGTFile, vector<string> &vMeshFileNames, vector<vtkTransform*> &vTransforms, vtkPolyDataCollection* meshes)
{
	DeleteTransforms(vTransforms);

	cout << "Loading ground truth...";
	ifstream file(strGTFile.c_str(), std::ios_base::in);

	if (!file.is_open() )
	{
		cout << "Not found: " << strGTFile << "\n";
		return;
	}

	string strText;
	file >> strText;
	file >> strText;
	file >> strText;

	int nClouds;
	file >> nClouds;
	
	vector<string> vMeshNames;
	vMeshNames.resize(vMeshFileNames.size());

	for(int na=0; na<(int)vMeshFileNames.size(); na++)
	{
		vMeshNames[na] =  GetRelativeName(vMeshFileNames[na].c_str());
	}

	vtkTransform* tempTransform = vtkTransform::New();
	vTransforms.resize(vMeshNames.size());
	int compRes;
	std::vector<bool> foundFiles;
	foundFiles.resize(vMeshNames.size());
	string meshName;
	for(int fo=0; fo<(int)foundFiles.size(); fo++)
	{
		foundFiles[fo] = false;
		vTransforms[fo] = NULL;
	}
	for(int tr=0; tr<nClouds; tr++)
	{
		
		file >> meshName;
		int idString = -1;
		for(int st=0; st<(int)vMeshNames.size(); st++)
		{
			compRes = vMeshNames[st].compare(meshName);
			if(compRes==0)
			{
				idString = st;
				break;
			}
		}
		if(idString != -1)
		{
			vTransforms[idString] = vtkTransform::New();
			foundFiles[idString] = true;
			file >> vTransforms[idString]->GetMatrix()->Element[0][0];
			file >> vTransforms[idString]->GetMatrix()->Element[0][1];
			file >> vTransforms[idString]->GetMatrix()->Element[0][2];
			file >> vTransforms[idString]->GetMatrix()->Element[0][3];
			file >> vTransforms[idString]->GetMatrix()->Element[1][0];
			file >> vTransforms[idString]->GetMatrix()->Element[1][1];
			file >> vTransforms[idString]->GetMatrix()->Element[1][2];
			file >> vTransforms[idString]->GetMatrix()->Element[1][3];
			file >> vTransforms[idString]->GetMatrix()->Element[2][0];
			file >> vTransforms[idString]->GetMatrix()->Element[2][1];
			file >> vTransforms[idString]->GetMatrix()->Element[2][2];
			file >> vTransforms[idString]->GetMatrix()->Element[2][3];
			file >> vTransforms[idString]->GetMatrix()->Element[3][0];
			file >> vTransforms[idString]->GetMatrix()->Element[3][1];
			file >> vTransforms[idString]->GetMatrix()->Element[3][2];
			file >> vTransforms[idString]->GetMatrix()->Element[3][3];
		}
		else
		{
			file >> tempTransform->GetMatrix()->Element[0][0];
			file >> tempTransform->GetMatrix()->Element[0][1];
			file >> tempTransform->GetMatrix()->Element[0][2];
			file >> tempTransform->GetMatrix()->Element[0][3];
			file >> tempTransform->GetMatrix()->Element[1][0];
			file >> tempTransform->GetMatrix()->Element[1][1];
			file >> tempTransform->GetMatrix()->Element[1][2];
			file >> tempTransform->GetMatrix()->Element[1][3];
			file >> tempTransform->GetMatrix()->Element[2][0];
			file >> tempTransform->GetMatrix()->Element[2][1];
			file >> tempTransform->GetMatrix()->Element[2][2];
			file >> tempTransform->GetMatrix()->Element[2][3];
			file >> tempTransform->GetMatrix()->Element[3][0];
			file >> tempTransform->GetMatrix()->Element[3][1];
			file >> tempTransform->GetMatrix()->Element[3][2];
			file >> tempTransform->GetMatrix()->Element[3][3];
		}
	}
	tempTransform->Delete();

	for(int tr=(int)vMeshFileNames.size()-1; tr>=0; tr--)
	{
		if(!foundFiles[tr])
		{
			cout << "View " << vMeshFileNames[tr].c_str() << "removed\n";
			vMeshFileNames.erase( vMeshFileNames.begin()+tr );
			if(vTransforms[tr])
			{
				vTransforms[tr]->Delete();
				vTransforms[tr] = NULL;
			}
			vTransforms.erase( vTransforms.begin()+tr );
			meshes->RemoveItem((int)tr);
		}
	}

	for(int tr = (int)vTransforms.size()-1; tr>=0; tr--)
	{
		if(!vTransforms[tr])
		{
			cout << "transform " << tr << "removed\n";
			vTransforms.erase( vTransforms.begin()+tr );
		}
	}


	file.close();

	cout << "OK\n";
}


void Pairwise3DRegistrationEvaluation::LoadGroundTruth(const string &strGTFile, vector<string> &vViewFileNames, vector<vtkTransform*> &vTransforms)
{
	ifstream file(strGTFile.c_str(), std::ios_base::in);

	if (!file.is_open() )
	{
		cout << "\nERROR: Impossible to open " << strGTFile << "\n";
		getchar();
		exit(-1);
	}

	cout << "Loading ground truth...";
	string strText;
	file >> strText;
	file >> strText;
	file >> strText;

	int nClouds;
	file >> nClouds;
	
	if((size_t)nClouds < vViewFileNames.size()) 
	{
		cout << "ERROR (Registration::LoadGroundTruth): nClouds in groundtruth.txt < vViewFileNames.size()   " << nClouds << "  <  " << vViewFileNames.size() << endl;
		getchar();
		exit(-1);
	}

	vector<string> vMeshNames;
	vMeshNames.resize(vViewFileNames.size());

	for(int na=0; na<(int)vViewFileNames.size(); na++)
	{
		vMeshNames[na] =  GetRelativeName(vViewFileNames[na].c_str());
	}

	vtkTransform* tempTransform = vtkTransform::New();
	vTransforms.resize(vMeshNames.size());
	int compRes;
	std::vector<bool> foundFiles;
	foundFiles.resize(vMeshNames.size());
	string meshName;
	for(int fo=0; fo<(int)foundFiles.size(); fo++)
	{
		foundFiles[fo] = false;
		vTransforms[fo] = NULL;
	}
	for(int tr=0; tr<nClouds; tr++)
	{
		
		file >> meshName;
		int idString = -1;
		for(int st=0; st<(int)vMeshNames.size(); st++)
		{
			compRes = vMeshNames[st].compare(meshName);
			if(compRes==0)
			{
				idString = st;
				break;
			}
		}
		if(idString != -1)
		{
			vTransforms[idString] = vtkTransform::New();
			foundFiles[idString] = true;
			file >> vTransforms[idString]->GetMatrix()->Element[0][0];
			file >> vTransforms[idString]->GetMatrix()->Element[0][1];
			file >> vTransforms[idString]->GetMatrix()->Element[0][2];
			file >> vTransforms[idString]->GetMatrix()->Element[0][3];
			file >> vTransforms[idString]->GetMatrix()->Element[1][0];
			file >> vTransforms[idString]->GetMatrix()->Element[1][1];
			file >> vTransforms[idString]->GetMatrix()->Element[1][2];
			file >> vTransforms[idString]->GetMatrix()->Element[1][3];
			file >> vTransforms[idString]->GetMatrix()->Element[2][0];
			file >> vTransforms[idString]->GetMatrix()->Element[2][1];
			file >> vTransforms[idString]->GetMatrix()->Element[2][2];
			file >> vTransforms[idString]->GetMatrix()->Element[2][3];
			file >> vTransforms[idString]->GetMatrix()->Element[3][0];
			file >> vTransforms[idString]->GetMatrix()->Element[3][1];
			file >> vTransforms[idString]->GetMatrix()->Element[3][2];
			file >> vTransforms[idString]->GetMatrix()->Element[3][3];
		}
		else
		{
			file >> tempTransform->GetMatrix()->Element[0][0];
			file >> tempTransform->GetMatrix()->Element[0][1];
			file >> tempTransform->GetMatrix()->Element[0][2];
			file >> tempTransform->GetMatrix()->Element[0][3];
			file >> tempTransform->GetMatrix()->Element[1][0];
			file >> tempTransform->GetMatrix()->Element[1][1];
			file >> tempTransform->GetMatrix()->Element[1][2];
			file >> tempTransform->GetMatrix()->Element[1][3];
			file >> tempTransform->GetMatrix()->Element[2][0];
			file >> tempTransform->GetMatrix()->Element[2][1];
			file >> tempTransform->GetMatrix()->Element[2][2];
			file >> tempTransform->GetMatrix()->Element[2][3];
			file >> tempTransform->GetMatrix()->Element[3][0];
			file >> tempTransform->GetMatrix()->Element[3][1];
			file >> tempTransform->GetMatrix()->Element[3][2];
			file >> tempTransform->GetMatrix()->Element[3][3];
		}
	}
	tempTransform->Delete();

	for(size_t tr=0; tr<vViewFileNames.size(); tr++)
	{
		if(!foundFiles[tr])
		{
			cout << "ERROR (Registration::LoadGroundTruth): " << vViewFileNames[tr] << " not found in groundtruth.txt" << endl;
			getchar();
			exit(-1);
		}
	}

	file.close();

	cout << "OK\n";
}

void Pairwise3DRegistrationEvaluation::DeleteTransforms(std::vector<vtkTransform*> &transforms)
{
	for(size_t tr=0; tr<transforms.size(); tr++)
	{
		if(transforms[tr])
		{
			transforms[tr]->Delete();
			transforms[tr] = NULL;
		}
	}
}


void Pairwise3DRegistrationEvaluation::CalcAllViewPairs(std::vector<std::pair<int, int> > &vMatchPairs, int nElements)
{
	int nPairs = (int)( (nElements) * ((nElements-1)/2.0) );

	vMatchPairs.resize(nPairs);

	std::vector<std::pair<int, int> >::iterator itPa = vMatchPairs.begin();
	for(int i=0; i<nElements-1; i++)
	{
		for(int j=i+1; j<nElements; j++)
		{
			itPa->first = i;
			itPa->second = j;
			itPa++;
		}
	}
}


vtkPolyData* Pairwise3DRegistrationEvaluation::ReadMesh(const std::string &absMeshFileName, bool convertToVtk )
{
	vtkPolyData *mesh;

	string absMeshFileName_vtk;
	absMeshFileName_vtk = ChangeExtension(absMeshFileName.c_str(), "vtk");
	
	bool filesExist_ply = ExistsFile(absMeshFileName.c_str());
	bool filesExist_vtk = ExistsFile(absMeshFileName_vtk.c_str());

	if(!filesExist_ply)
	{
		cout << "No ply files to load (" << filesExist_ply << ")"; 
		getchar();
		exit(-1);
	}
	if(filesExist_vtk && convertToVtk )
	{
		//load mesh in vtk format
		cout << "Loading cloud in vtk format...";
		mesh = LoadPolyData(absMeshFileName_vtk.c_str());
		cout << "OK\n";
	}
	else
	{
		//load mesh in ply format
		cout << "Loading cloud in ply format...";
		mesh = LoadPolyData(absMeshFileName.c_str());
		cout << "OK\n";

		Preprocessing(mesh);

		//save mesh in vtk format
		if(convertToVtk)
		{
			WriteVtk(absMeshFileName_vtk.c_str(), mesh); 
		}
	}

	return mesh;

}


void Pairwise3DRegistrationEvaluation::ComputeGroundTruth_All(vector<vtkTransform*> &vTransforms_GT, const std::vector<std::pair<int,int>> &vMeshPairs, vector<vtkTransform*> &vPairwiseTransforms_GT)
{
	DeleteTransforms(vPairwiseTransforms_GT);
	vPairwiseTransforms_GT.resize(vMeshPairs.size(), NULL);

	for(size_t pa=0; pa<vMeshPairs.size(); pa++)
	{
		if( (vMeshPairs[pa].first < (int)vTransforms_GT.size()) && (vMeshPairs[pa].second < (int)vTransforms_GT.size()) )
		{
			vPairwiseTransforms_GT[pa] = ComputeGroundTruth(vTransforms_GT[vMeshPairs[pa].first], vTransforms_GT[vMeshPairs[pa].second]);
		}
	}
}


void Pairwise3DRegistrationEvaluation::Preprocessing(vtkPolyData* &polyData)
{
	cout << "Computing normals...";
	CalcMeshNormals(polyData);
	CleanPolyData(polyData, true, 0.0);
	cout << "OK" << endl;
}

vtkTransform* Pairwise3DRegistrationEvaluation::ComputeGroundTruth(vtkTransform* &transf_trg, vtkTransform* &transf_ref)
{
	vtkTransform* pairTransform2 = vtkTransform::New();
	pairTransform2->PostMultiply();
	vtkTransform* inversePairTransform1;

	pairTransform2->Concatenate(transf_ref);
	inversePairTransform1 = CreateInverseTransform(transf_trg);
	pairTransform2->Concatenate( inversePairTransform1);
	inversePairTransform1->Delete();
	return pairTransform2;
}


vtkTransform* Pairwise3DRegistrationEvaluation::CreateInverseTransform(vtkTransform* transform)
{
	//inverse matrix
	vtkMatrix4x4* invTrans2 = vtkMatrix4x4::New();
	vtkMatrix4x4::Invert(transform->GetMatrix(), invTrans2);

	vtkTransform* inverseTransform = vtkTransform::New();
	inverseTransform->PostMultiply();
	inverseTransform->Concatenate(invTrans2);
	invTrans2->Delete();
	return inverseTransform;
}


int Pairwise3DRegistrationEvaluation::GetBoundaryPoints(vtkPolyData *polydata, bool* &boundaryPointsIds)
{
	boundaryPointsIds = new bool[polydata->GetNumberOfPoints()];
	for(int po=0; po<polydata->GetNumberOfPoints(); po++)
	{
		boundaryPointsIds[po] = false;
	}
	vtkIdType* ptrCells = polydata->GetPolys()->GetPointer();
	int nVertex, ve1, ve2;
	vtkIdList* idList = vtkIdList::New();
	polydata->BuildLinks();
	//for every cell in polydata
	for(int ce=0; ce<polydata->GetNumberOfPolys(); ce++)
	{
		nVertex = *ptrCells;
		ptrCells++;
		//for every edge in polydata
		for(int ve=0; ve<nVertex; ve++)
		{
			ve1 = *(ptrCells + ve);
			ve2 = *(ptrCells + ((ve+1)%nVertex) );
			polydata->GetCellEdgeNeighbors(ce, ve1, ve2,idList);
			if(idList->GetNumberOfIds()<1)
			{
				assert(ve1 >= 0 && ve1 < polydata->GetNumberOfPoints());
				assert(ve2 >= 0 && ve2 < polydata->GetNumberOfPoints());
				boundaryPointsIds[ve1] = true;
				boundaryPointsIds[ve2] = true;
			}
		}
		ptrCells += nVertex;
	}
	idList->Delete();
	return polydata->GetNumberOfPoints();
}








Pairwise3DRegistrationEvaluation::KdTree::KdTree(const bool useTrueKdTree)
{
	if(useTrueKdTree)
	{
		m_kdTree = vtkKdTreePointLocator::New();
	}
	else
	{
		m_kdTree = vtkPointLocator::New();
	}
	m_pointsList = vtkIdList::New();
}


Pairwise3DRegistrationEvaluation::KdTree::KdTree(vtkPolyData* polyData, const bool useTrueKdTree)
{
	if(useTrueKdTree)
	{
		m_kdTree = vtkKdTreePointLocator::New();
	}
	else
	{
		m_kdTree = vtkPointLocator::New();
	}
	m_pointsList = vtkIdList::New();
	SetPolyData(polyData);
}

void Pairwise3DRegistrationEvaluation::KdTree::SetPolyData(vtkPolyData* polyData)
{
	m_polyData = polyData;
	m_kdTree->SetDataSet(polyData);
	m_kdTree->Update();
}

Pairwise3DRegistrationEvaluation::KdTree::~KdTree()
{
	m_kdTree->Delete();
	m_pointsList->Delete();

}


vtkIdList* Pairwise3DRegistrationEvaluation::KdTree::FindPointsWithinRadius(const float* const point, float radius )
{
	const double p[3] = {point[0], point[1], point[2]};
	m_kdTree->FindPointsWithinRadius(radius, p, m_pointsList);
	return m_pointsList;
}

vtkIdList* Pairwise3DRegistrationEvaluation::KdTree::FindPointsWithinRadius(const int pointIndex, float radius )
{
	float* point = (float*)(m_polyData->GetPoints()->GetVoidPointer(pointIndex*3));
	const double p[3] = {point[0], point[1], point[2]};
	m_kdTree->FindPointsWithinRadius(radius, p, m_pointsList);
	return m_pointsList;
}



vtkIdList* Pairwise3DRegistrationEvaluation::KdTree::FindPointsWithinRadius(const double* point, double radius )
{
	const double p[3] = {point[0], point[1], point[2]};

	m_kdTree->FindPointsWithinRadius(radius, p, m_pointsList);
	return m_pointsList;
}

int Pairwise3DRegistrationEvaluation::KdTree::FindNearestPoint(double* point)
{
	return m_kdTree->FindClosestPoint(point);
}

int Pairwise3DRegistrationEvaluation::KdTree::FindNearestPoint(float* point)
{
	const double p[3] = {point[0], point[1], point[2]};
	return m_kdTree->FindClosestPoint(p);
}

int Pairwise3DRegistrationEvaluation::KdTree::FindNearestPoint(const int pointIndex)
{
	float* point = (float*)(m_polyData->GetPoints()->GetVoidPointer(pointIndex*3));
	const double p[3] = {point[0], point[1], point[2]};
	return m_kdTree->FindClosestPoint(p);
}

int Pairwise3DRegistrationEvaluation::KdTree::FindNearestPointWithinRadius(double* point, double radius, double & dist)
{
	return m_kdTree->FindClosestPointWithinRadius(radius, point, dist);
}

int Pairwise3DRegistrationEvaluation::KdTree::FindNearestPointWithinRadius(float* point, float radius, double & dist)
{
	double doublePoint[] = {point[0],point[1],point[2]};
	return m_kdTree->FindClosestPointWithinRadius(radius, doublePoint, dist);
}

int Pairwise3DRegistrationEvaluation::KdTree::FindNearestPointWithinRadius(const int pointIndex, float radius, double & dist)
{
	float* point = (float*)(m_polyData->GetPoints()->GetVoidPointer(pointIndex*3));
	double doublePoint[] = {point[0],point[1],point[2]};
	return m_kdTree->FindClosestPointWithinRadius(radius, doublePoint, dist);
}


double Pairwise3DRegistrationEvaluation::KdTree::FindNearestPoint(double* point, int &nearestPointId)
{
	nearestPointId = m_kdTree->FindClosestPoint(point);
	double nearestPoint[3];
	m_polyData->GetPoint(nearestPointId, nearestPoint);
	return sqrt(vtkMath::Distance2BetweenPoints(nearestPoint, point));
}

double Pairwise3DRegistrationEvaluation::KdTree::FindNearestPoint(float* point, int &nearestPointId)
{
	const double p[3] = {point[0], point[1], point[2]};
	nearestPointId = m_kdTree->FindClosestPoint(p);
	float* nearestPoint = GetPolyDataPointsPointer(m_polyData, nearestPointId);
	return sqrt(vtkMath::Distance2BetweenPoints(nearestPoint, point));
}

double Pairwise3DRegistrationEvaluation::KdTree::FindNearestPoint(const int pointIndex, int &nearestPointId)
{
	float* point = (float*)(m_polyData->GetPoints()->GetVoidPointer(pointIndex*3));
	const double p[3] = {point[0], point[1], point[2]};
	nearestPointId = m_kdTree->FindClosestPoint(p);
	float* nearestPoint = GetPolyDataPointsPointer(m_polyData, nearestPointId);
	return sqrt(vtkMath::Distance2BetweenPoints(nearestPoint, point));
}

vtkIdList* Pairwise3DRegistrationEvaluation::KdTree::FindNearestNPoints(int n, const double* point)
{
	const double p[3] = {point[0], point[1], point[2]};
	m_kdTree->FindClosestNPoints(n, p, m_pointsList);
	return m_pointsList;
}

vtkIdList* Pairwise3DRegistrationEvaluation::KdTree::FindNearestNPoints(int n, const float* point)
{
	const double p[3] = {point[0], point[1], point[2]};
	m_kdTree->FindClosestNPoints(n, p, m_pointsList);
	return m_pointsList;
}

vtkIdList* Pairwise3DRegistrationEvaluation::KdTree::FindNearestNPoints(int n, const int pointIndex)
{
	float* point = (float*)(m_polyData->GetPoints()->GetVoidPointer(pointIndex*3));
	const double p[3] = {point[0], point[1], point[2]};

	m_kdTree->FindClosestNPoints(n, p, m_pointsList);
	return m_pointsList;
}

float* Pairwise3DRegistrationEvaluation::GetPolyDataPointsPointer(vtkPolyData* polyData, const unsigned int iPoint)
{
	vtkPoints* Array = polyData->GetPoints();

	if(Array->GetNumberOfPoints() == 0 )
		return NULL;
	else
		return (float*)(Array->GetVoidPointer(iPoint*3));
}


void Pairwise3DRegistrationEvaluation::ApplyTransform(float* points, int nPoints, const double* matrTransf)
{
	float* ptrPoints = points;
	float pointIn[4];
	float pointOut[4];
	pointIn[3] = 1;
	for(int po=0; po<nPoints; po++)
	{
		memcpy(pointIn, ptrPoints, sizeof(float)*3);
		vtkMatrix4x4::MultiplyPoint(matrTransf, pointIn, pointOut);
		memcpy(ptrPoints, pointOut, sizeof(float)*3);
		ptrPoints+=3;
	}
}


void Pairwise3DRegistrationEvaluation::GetCentroid(float* points, int nPoints, float centroid[3])
{
	double centroid_double[] = {0.0, 0.0, 0.0};
	float* ptrPoints = points;
	for(int po=0; po<nPoints; po++)
	{
		centroid_double[0] += *ptrPoints;
		ptrPoints++;
		centroid_double[1] += *ptrPoints;
		ptrPoints++;
		centroid_double[2] += *ptrPoints;
		ptrPoints++;
	}
	centroid[0] = (float)(centroid_double[0]/nPoints);
	centroid[1] = (float)(centroid_double[1]/nPoints);
	centroid[2] = (float)(centroid_double[2]/nPoints);
}

void Pairwise3DRegistrationEvaluation::Translate(float* points, int nPoints, const float point[3])
{
	float* ptrPoints = points;
	for(int po=0; po<nPoints; po++)
	{
		*ptrPoints += point[0];
		ptrPoints++;
		*ptrPoints += point[1];
		ptrPoints++;
		*ptrPoints += point[2];
		ptrPoints++;
	}
}



#ifdef __GNUC__
Pairwise3DRegistrationEvaluation::CpuTimeProfiler::CpuTimeProfiler(clockid_t clk_id)
  : m_clk_id(clk_id)
#else
Pairwise3DRegistrationEvaluation::CpuTimeProfiler::CpuTimeProfiler()
#endif
{
	#if defined (_MSC_VER)
		QueryPerformanceFrequency( &m_frequency );
	#endif
	GetActualTime();	
}

void Pairwise3DRegistrationEvaluation::CpuTimeProfiler::GetActualTime()
{
	#if defined (_MSC_VER)
		GetProcessTimes(GetCurrentProcess(), &m_actualCreateTime, &m_actualExitTime, &m_actualKernelTime, &m_actualUserTime);

		ULARGE_INTEGER actualkernel, actualUser;
		actualkernel.LowPart = m_actualKernelTime.dwLowDateTime;
		actualkernel.HighPart = m_actualKernelTime.dwHighDateTime;

		actualUser.LowPart = m_actualUserTime.dwLowDateTime;
		actualUser.HighPart = m_actualUserTime.dwHighDateTime;

		m_actualTime.QuadPart = actualkernel.QuadPart + actualUser.QuadPart;
  #elif defined __GNUC__
    clock_gettime(m_clk_id, &m_actualTime);
	#endif
}

double Pairwise3DRegistrationEvaluation::CpuTimeProfiler::GetElapsedHours()
{
	return GetElapsedMins() / 60.0;
}

double Pairwise3DRegistrationEvaluation::CpuTimeProfiler::GetElapsedMins()
{
	return GetElapsedSecs() / 60.0;
}

double Pairwise3DRegistrationEvaluation::CpuTimeProfiler::GetElapsedSecs()
{
	#if defined (_MSC_VER)
		FILETIME stopUserTime, stopKernelTime, stopCreateTime, stopExitTime;
		GetProcessTimes(GetCurrentProcess(), &stopCreateTime, &stopExitTime, &stopKernelTime, &stopUserTime);

		ULARGE_INTEGER stopKernel, stopUser, elapsed;
		stopUser.LowPart = stopUserTime.dwLowDateTime;
		stopUser.HighPart = stopUserTime.dwHighDateTime;

		stopKernel.LowPart = stopKernelTime.dwLowDateTime;
		stopKernel.HighPart = stopKernelTime.dwHighDateTime;

		elapsed.QuadPart = stopUser.QuadPart + stopKernel.QuadPart - m_actualTime.QuadPart;

		// the output of GetProcessTimes is expressed in 100-nanosecond time units so we have to divide by 10'000'000 to obtain sec
		return (double)elapsed.QuadPart / 10000000.0 ;
  #elif defined __GNUC__
    return GetElapsedMilli() / 1000;
	#else
		return 0;
	#endif
}

double Pairwise3DRegistrationEvaluation::CpuTimeProfiler::GetElapsedMilli()
{
	#if defined (_MSC_VER)
		FILETIME stopUserTime, stopKernelTime, stopCreateTime, stopExitTime;
		GetProcessTimes(GetCurrentProcess(), &stopCreateTime, &stopExitTime, &stopKernelTime, &stopUserTime);

		ULARGE_INTEGER stopKernel, stopUser, elapsed;
		stopUser.LowPart = stopUserTime.dwLowDateTime;
		stopUser.HighPart = stopUserTime.dwHighDateTime;

		stopKernel.LowPart = stopKernelTime.dwLowDateTime;
		stopKernel.HighPart = stopKernelTime.dwHighDateTime;

		elapsed.QuadPart = stopUser.QuadPart + stopKernel.QuadPart - m_actualTime.QuadPart;

		// the output of GetProcessTimes is expressed in 100-nanosecond time units so we have to divide by 10'000 to obtain ms
		return (double)elapsed.QuadPart / 10000.0 ;
  #elif defined __GNUC__
    timespec stop;
    clock_gettime(m_clk_id, &stop);
    return (difftime(stop.tv_sec, m_actualTime.tv_sec) * 1000.0 + (double)(stop.tv_nsec - m_actualTime.tv_nsec) * 1e-6);
	#else
		return 0;
	#endif
}

string Pairwise3DRegistrationEvaluation::RemoveExt(const string &filename)
{
	std::string strTemp = filename;
	unsigned int dotPos = (int)strTemp.find_last_of(".");
	if(dotPos != std::string::npos)
	{
		strTemp.erase(dotPos);
	}
	return strTemp;
}


string Pairwise3DRegistrationEvaluation::GetRelativeName_WithoutExt(const string &filename)
{
	string relFilename = GetRelativeName(filename);
	return RemoveExt(relFilename);
}

double Pairwise3DRegistrationEvaluation::CalcOverlappingAreaMax(vtkPolyData* polyData1, vtkPolyData* polyData2, double maxDistance, double absNormalDotThresh, bool getOverlappedPoints, bool* overlappedPoints1, bool* overlappedPoints2 )
{
	vtkIdType* ptrFirstTriangles = polyData1->GetPolys()->GetPointer();
	vtkIdType* ptrSecondTriangles = polyData2->GetPolys()->GetPointer();

	//create two polydata with the points that are the centers of the triangles of first and transformedPoly polydata
	//and calc area for every triangle
	vtkPolyData* triangleCenters1 = vtkPolyData::New();
	vtkPolyData* triangleCenters2 = vtkPolyData::New();
	float* ptrTriangleCenters1Points = AllocatePoints(triangleCenters1, polyData1->GetNumberOfPolys());
	float* ptrTriangleCenters2Points = AllocatePoints(triangleCenters2, polyData2->GetNumberOfPolys());
	vtkIdType *ptr1, *ptr2, *ptr3;
	double v1[3], v2[3], v3[3], center[3];
	double area=0, totFirstArea=0, overFirstArea=0, totSecondArea=0, overSecondArea=0;
	for(int tr=0; tr<polyData1->GetNumberOfPolys(); tr++)
	{
		ptrFirstTriangles++;
		ptr1 = ptrFirstTriangles++;
		ptr2 = ptrFirstTriangles++;
		ptr3 = ptrFirstTriangles++;
		polyData1->GetPoint(*ptr1, v1);
		polyData1->GetPoint(*ptr2, v2);
		polyData1->GetPoint(*ptr3, v3);

		//find triangle center
		vtkTriangle::TriangleCenter(v1, v2, v3, center);
		*ptrTriangleCenters1Points++ = center[0];
		*ptrTriangleCenters1Points++ = center[1];
		*ptrTriangleCenters1Points++ = center[2];

		//find triangle area
		area = vtkTriangle::TriangleArea(v1, v2, v3);
		totFirstArea += area;
	}

	for(int tr=0; tr<polyData2->GetNumberOfPolys(); tr++)
	{
		ptrSecondTriangles++;
		ptr1 = ptrSecondTriangles++;
		ptr2 = ptrSecondTriangles++;
		ptr3 = ptrSecondTriangles++;
		polyData2->GetPoint(*ptr1, v1);
		polyData2->GetPoint(*ptr2, v2);
		polyData2->GetPoint(*ptr3, v3);

		//find triangle center
		vtkTriangle::TriangleCenter(v1, v2, v3, center);
		*ptrTriangleCenters2Points++ = center[0];
		*ptrTriangleCenters2Points++ = center[1];
		*ptrTriangleCenters2Points++ = center[2];

		//find triangle area
		area = vtkTriangle::TriangleArea(v1, v2, v3);
		totSecondArea += area;
	}

	KdTree firstKdTree(triangleCenters1);
	KdTree secondKdTree(triangleCenters2);

	//nearTriangles1 and nearTriangles2 indicate triangles nears to other polydata
	bool* nearTriangles1 = new bool[polyData1->GetNumberOfPolys()];
	bool* nearTriangles2 = new bool[polyData2->GetNumberOfPolys()];
	for(int tr=0; tr<polyData1->GetNumberOfPolys(); tr++)
	{
		nearTriangles1[tr] = false;
	}
	for(int tr=0; tr<polyData2->GetNumberOfPolys(); tr++)
	{
		nearTriangles2[tr] = false;
	}

	double normal[3];
	double otherNormal[3];
	double otherv1[3], otherv2[3], otherv3[3];

	vtkIdList* idList;
	ptrTriangleCenters1Points = (float*)triangleCenters1->GetPoints()->GetVoidPointer(0);
	ptrTriangleCenters2Points = (float*)triangleCenters2->GetPoints()->GetVoidPointer(0);
	ptrFirstTriangles = polyData1->GetPolys()->GetPointer();
	ptrSecondTriangles = polyData2->GetPolys()->GetPointer();
	//find every triangle in first that is near to transformedPoly
	for(int tr=0; tr<polyData1->GetNumberOfPolys(); tr++)
	{
		idList = secondKdTree.FindPointsWithinRadius(ptrTriangleCenters1Points, maxDistance);
		if(idList->GetNumberOfIds()>0)
		{
			//check if at least one triangle in the other polydata has similar orientation
			ptr1 = ptrFirstTriangles + 4*tr +1;
			ptr2 = ptrFirstTriangles + 4*tr +2;
			ptr3 = ptrFirstTriangles + 4*tr +3;
			polyData1->GetPoint(*ptr1, v1);
			polyData1->GetPoint(*ptr2, v2);
			polyData1->GetPoint(*ptr3, v3);
			vtkTriangle::ComputeNormal(v1, v2, v3, normal);

			for(int ce=0; ce<idList->GetNumberOfIds(); ce++)
			{
				ptr1 = ptrSecondTriangles + 4*idList->GetId(ce) +1;
				ptr2 = ptrSecondTriangles + 4*idList->GetId(ce) +2;
				ptr3 = ptrSecondTriangles + 4*idList->GetId(ce) +3;
				polyData2->GetPoint(*ptr1, otherv1);
				polyData2->GetPoint(*ptr2, otherv2);
				polyData2->GetPoint(*ptr3, otherv3);
				vtkTriangle::ComputeNormal(otherv1, otherv2, otherv3, otherNormal);
				if( abs(vtkMath::Dot(normal, otherNormal)) > absNormalDotThresh )
				{
					nearTriangles1[tr] = true;
					nearTriangles2[idList->GetId(ce)] = true;
					break;
				}
			}
			if(nearTriangles1[tr])
			{
				//find triangle area
				area = vtkTriangle::TriangleArea(v1, v2, v3);
				overFirstArea += area;
			}
		}


		ptrTriangleCenters1Points+=3;
	}

	//find every triangle in transformedPoly that is near to first
	for(int tr=0; tr<polyData2->GetNumberOfPolys(); tr++)
	{
		if(!nearTriangles2[tr])
		{
			idList = firstKdTree.FindPointsWithinRadius(ptrTriangleCenters2Points, maxDistance);
			if(idList->GetNumberOfIds()>0)
			{
				ptr1 = ptrSecondTriangles + 4*tr +1;
				ptr2 = ptrSecondTriangles + 4*tr +2;
				ptr3 = ptrSecondTriangles + 4*tr +3;
				polyData2->GetPoint(*ptr1, v1);
				polyData2->GetPoint(*ptr2, v2);
				polyData2->GetPoint(*ptr3, v3);
				vtkTriangle::ComputeNormal(v1, v2, v3, normal);

				for(int ce=0; ce<idList->GetNumberOfIds(); ce++)
				{
					ptr1 = ptrFirstTriangles + 4*idList->GetId(ce) +1;
					ptr2 = ptrFirstTriangles + 4*idList->GetId(ce) +2;
					ptr3 = ptrFirstTriangles + 4*idList->GetId(ce) +3;
					polyData1->GetPoint(*ptr1, otherv1);
					polyData1->GetPoint(*ptr2, otherv2);
					polyData1->GetPoint(*ptr3, otherv3);
					vtkTriangle::ComputeNormal(otherv1, otherv2, otherv3, otherNormal);
					if( abs(vtkMath::Dot(normal, otherNormal)) > absNormalDotThresh )
					{
						nearTriangles2[tr] = true;
						break;
					}
				}
				if(nearTriangles2[tr])
				{
					//find triangle area
					area = vtkTriangle::TriangleArea(v1, v2, v3);
					overSecondArea += area;
				}
			}
		}
		else
		{
			//find triangle area
			ptr1 = ptrSecondTriangles + 4*tr +1;
			ptr2 = ptrSecondTriangles + 4*tr +2;
			ptr3 = ptrSecondTriangles + 4*tr +3;
			polyData2->GetPoint(*ptr1, v1);
			polyData2->GetPoint(*ptr2, v2);
			polyData2->GetPoint(*ptr3, v3);
			area = vtkTriangle::TriangleArea(v1, v2, v3);
			overSecondArea += area;
		}


		ptrTriangleCenters2Points+=3;
	}

	double overlapPerc = Max( overFirstArea/totFirstArea, overSecondArea/totSecondArea);


	if(getOverlappedPoints)
	{
		memset(overlappedPoints1, 0, sizeof(bool)*polyData1->GetNumberOfPoints());
		memset(overlappedPoints2, 0, sizeof(bool)*polyData2->GetNumberOfPoints());

		ptrFirstTriangles = polyData1->GetPolys()->GetPointer();
		ptrSecondTriangles = polyData2->GetPolys()->GetPointer();

		for(int tr=0; tr<polyData1->GetNumberOfPolys(); tr++)
		{
			if(nearTriangles1[tr])
			{
				ptr1 = ptrFirstTriangles + 4*tr +1;
				ptr2 = ptrFirstTriangles + 4*tr +2;
				ptr3 = ptrFirstTriangles + 4*tr +3;
				overlappedPoints1[*ptr1] = true;
				overlappedPoints1[*ptr2] = true;
				overlappedPoints1[*ptr3] = true;
			}
		}

		for(int tr=0; tr<polyData2->GetNumberOfPolys(); tr++)
		{
			if(nearTriangles2[tr])
			{
				ptr1 = ptrSecondTriangles + 4*tr +1;
				ptr2 = ptrSecondTriangles + 4*tr +2;
				ptr3 = ptrSecondTriangles + 4*tr +3;
				overlappedPoints2[*ptr1] = true;
				overlappedPoints2[*ptr2] = true;
				overlappedPoints2[*ptr3] = true;
			}
		}
	}


	triangleCenters1->Delete();
	triangleCenters2->Delete();

	delete[] nearTriangles1;
	delete[] nearTriangles2;

	return overlapPerc;
}


vtkIdType* Pairwise3DRegistrationEvaluation::GetPolyDataTrianglesPointer(vtkPolyData* polyData, const unsigned int iTriangle)
{
	vtkCellArray* Array = polyData->GetPolys();

	if(Array->GetNumberOfCells() == 0 )
		return NULL;
	else
		return Array->GetPointer() + iTriangle*4;

	//return polyData->GetPolys()->GetPointer() + iTriangle*4;
}

vtkIdType Pairwise3DRegistrationEvaluation::GetPolyDataNumberOfPoints(vtkPolyData* polyData)
{
	if(polyData == NULL)
	{
		return -1;
	}
	return polyData->GetPoints()->GetNumberOfPoints();
}

size_t Pairwise3DRegistrationEvaluation::Sample_random(std::vector<bool> &vData, const float samplingFactor, const unsigned int seed)
{
	if(seed != std::numeric_limits<unsigned int>::max())
	{
		srand(seed);
	}

	size_t sampledSize =  size_t( samplingFactor * (float)(vData.size()) );

	if(sampledSize < 0)
	{
		throw runtime_error("ERROR (Sample_random): (sampledSize < 0)");
	}

	if(sampledSize > vData.size())
	{
		throw runtime_error("ERROR (Sample_random): (sampledSize > vData.size())");
	}

	for(size_t sa=0; sa<sampledSize; sa++)
	{
		vData[sa] = true;
	}
	for(size_t sa=sampledSize; sa<vData.size(); sa++)
	{
		vData[sa] = false;
	}
	
	random_shuffle(vData.begin(), vData.end(), NextRandom);

	return sampledSize;
}