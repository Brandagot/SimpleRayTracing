/*

Copyright (c) 2020, Dr Franck P. Vidal (f.vidal@bangor.ac.uk),
http://www.fpvidal.net/
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the Bangor University nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/


/**
********************************************************************************
*
*   @file       main.cxx
*
*   @brief      A simple ray-tracer without parallelism.
*
*   @version    1.0
*
*   @date       14/10/2020
*
*   @author     Dr Franck P. Vidal
*
*   License
*   BSD 3-Clause License.
*
*   For details on use and redistribution please refer
*   to http://opensource.org/licenses/BSD-3-Clause
*
*   Copyright
*   (c) by Dr Franck P. Vidal (f.vidal@bangor.ac.uk),
*   http://www.fpvidal.net/, Oct 2020, 2020, version 1.0, BSD 3-Clause License
*
********************************************************************************
*/


//******************************************************************************
//  Include
//******************************************************************************
#include <iostream>  // for cerr
#include <exception> // to catch exceptions
#include <algorithm> // for max
#include <cmath>     // for pow
#include <limits>    // for inf
#include <stdexcept> // for exceptions
#include <sstream>   // to format error messages
#include <string>
#include <chrono>	 // For timestamps
#include <pthread.h>
#include <vector>
#include <algorithm>

#include <assimp/Importer.hpp>  // C++ importer interface
#include <assimp/scene.h>       // Output data structure
#include <assimp/postprocess.h> // Post processing flags

#ifndef __Vec3_h
#include "Vec3.h"
#endif

#ifndef __Ray_h
#include "Ray.h"
#endif

#ifndef __Ray_h
#include "Ray.h"
#endif

#ifndef __TriangleMesh_h
#include "TriangleMesh.h"
#endif

#ifndef __Material_h
#include "Material.h"
#endif

#ifndef __Image_h
#include "Image.h"
#endif


//******************************************************************************
//  Namespace
//******************************************************************************
using namespace std;

//******************************************************************************
//  Classes
//******************************************************************************

/**
 * Class for storing data required by threads.
 * 
 * Using a class instead of a struct so we can work with references, as they need to be initialized by the constructor.
 */
class PThreadData{
	// Thread information
	// Using a class constructor to initialize the parameters to we can use references

	private:
		PThreadData(); // Private so we can't call the default instructor.
	public:
		// Initialising the values of this data class - Required when using references
		PThreadData(Image& anOutputImage,
					const vector<TriangleMesh>& aTriangleMeshSet,
					const Vec3& aDetectorPosition,
					const Vec3& aRayOrigin,
					const Vec3& anUpVector,
					const Vec3& aRightVector,
					const Light& aLight):
					m_output_image(anOutputImage),
					m_triangle_mesh_set(aTriangleMeshSet),
					m_detector_position(aDetectorPosition),
					m_ray_origin(aRayOrigin),
					m_up_vector(anUpVector),
					m_right_vector(aRightVector),
					m_light(aLight),
					m_thread_int_id(-1),
					m_start_pixel(-1),
					m_end_pixel(-1){}

		pthread_t m_thread_id;
		unsigned int m_thread_int_id;

		// Reference to variables needed in the render loop.
		Image& m_output_image;
		const vector<TriangleMesh>& m_triangle_mesh_set;
		const Vec3& m_detector_position;
		const Vec3& m_ray_origin;
		const Vec3& m_up_vector;
		const Vec3& m_right_vector;
		const Light& m_light;

		// Start and end point for data assigned to thread - Thread will work on start pixel through to end pixel
		unsigned int m_start_pixel;
		unsigned int m_end_pixel;
};

// Callback function for PThreads - a thread renders from pixel Ni to Nj, loads is distributed in main
void* renderLoopCallBack(void* apData);

// Help message for command-line args
void showUsage(const std::string& aProgramName);

// Retrieves command-line args and parses them
void processCmd(
	int argc, char** argv,
	string& aFileName,
	unsigned int& aWidth,
	unsigned int& aHeight,
	unsigned char& r,
	unsigned char& g,
	unsigned char& b,
	unsigned int& number_of_threads
);

// Not yet sure how this works - but need to apply shading based on ray-tracing data 
Vec3 applyShading(
	const Light& aLight,
	const Material& aMaterial,
	const Vec3& aNormalVector,
	const Vec3& aPosition,
	const Vec3& aViewPostion
);

// Retrieves a .ply mesh from file and adds it to the set
void loadMeshes(
	const std::string& aFileName,
	vector<TriangleMesh>& aMeshSet
);

// Might be obsolete now that we're going to create a view plane instead
TriangleMesh createBackground(
	const Vec3& anUpperBBoxCorner,
	const Vec3& aLowerBBoxCorner
);

// Gets the boudning box
void getBBox(
	const vector<TriangleMesh>& aMeshSet,
	Vec3& anUpperBBoxCorner,
	Vec3& aLowerBBoxCorner
);

//******************************************************************************
//  Constant global variables
//******************************************************************************
const unsigned int g_default_image_width = 2048;
const unsigned int g_default_image_height = 2048;

// Color constants
const Vec3 g_black(0, 0, 0);
const Vec3 g_white(1, 1, 1);

const Vec3 g_red(1, 0, 0);
const Vec3 g_green(0, 1, 0);
const Vec3 g_blue(0, 0, 1);

const Vec3 g_background_color = g_black; // But the background behind the viewplane appreas more grey, so im not sure aout this

int main(int argc, char** argv){
	try{
		// defualt output file
		string output_file_name = "test.jpg";

		// default image size - changed later if passed processCMD
		unsigned int image_width = g_default_image_width;
		unsigned int image_height= g_default_image_height;

		// Default background color is set to grey - might try size the viewplane to match the entire view
		unsigned char r = 128;
		unsigned char g = 128;
		unsigned char b = 128;
		
		// Default number of threads
		unsigned int number_of_threads = 1;

		// updates above defaults if passed command args
		processCmd(argc, argv, 
			output_file_name,
			image_width, image_height,
			r, g, b, number_of_threads);

		// Loading polygon meshes
		cout << "Loading polygon meshes" << endl;
		
		vector<TriangleMesh> p_mesh_set;
		loadMeshes("./dragon.ply", p_mesh_set);



	} 
	// Catch exceptions and error messages
	// Catches thrown exceptions, strings and char*
	catch(const std::exception& e){
		std::cerr << "ERROR: " << e.what() << std::endl;
	}
	catch(const std::string& e)
	{
		std::cerr << "ERROR: " << e << std::endl;
	}
	catch(const char* e)
	{
		std::cerr << "ERROR: " << e << std::endl;
	}

	return 0;
}

//---------------------------------------------
void showUsage(const std::string& aProgramName)
//---------------------------------------------
{
	std::cerr << "Usage: " << aProgramName << " <option(s)>" << endl <<
		"Options:" << endl <<
		"\t-h,--help\t\t\tShow this help message" << endl <<
		"\t-t,--threads T\tSpecify the number of threads (default value: 4)" << endl << 
		"\t-s,--size IMG_WIDTH IMG_HEIGHT\tSpecify the image size in number of pixels (default values: 2048 2048)" << endl << 
		"\t-b,--background R G B\t\tSpecify the background colour in RGB, acceptable values are between 0 and 255 (inclusive) (default values: 128 128 128)" << endl << 
		"\t-j,--jpeg FILENAME\t\tName of the JPEG file (default value: test.jpg)" << endl << 
		std::endl;
}

//-------------------------------------------------------------------
void processCmd(int argc, char** argv,
				string& aFileName,
				unsigned int& aWidth, unsigned int& aHeight,
				unsigned char& r, unsigned char& g, unsigned char& b,
				unsigned int& number_of_threads)
//-------------------------------------------------------------------
{
	// Process the command line
	int i = 1;
	while (i < argc)
	{
		std::string arg = argv[i];
		
		if (arg == "-h" || arg == "--help")
		{
			showUsage(argv[0]);
			exit(EXIT_SUCCESS);
		}
		else if (arg == "-s" || arg == "--size")
		{
			++i;
			if (i < argc)
			{
				aWidth = stoi(argv[i]);
			}
			else
			{
				showUsage(argv[0]);
				exit(EXIT_FAILURE);
			}
			
			++i;
			if (i < argc)
			{
				aHeight = stoi(argv[i]);
			}
			else
			{
				showUsage(argv[0]);
				exit(EXIT_FAILURE);
			}
		}
		else if (arg == "-b" || arg == "--background")
		{
			++i;
			if (i < argc)
			{
				r = stoi(argv[i]);
			}
			else
			{
				showUsage(argv[0]);
				exit(EXIT_FAILURE);
			}

			++i;
			if (i < argc)
			{
				g = stoi(argv[i]);
			}
			else
			{
				showUsage(argv[0]);
				exit(EXIT_FAILURE);
			}
			
			++i;
			if (i < argc)
			{
				b = stoi(argv[i]);
			}
			else
			{
				showUsage(argv[0]);
				exit(EXIT_FAILURE);
			}
		}
		else if (arg == "-j" || arg == "--jpeg")
		{                
			++i;
			if (i < argc)
			{
				aFileName = ("./out/" +  (string)argv[i]);
			}
			else
			{
				showUsage(argv[0]);
				exit(EXIT_FAILURE);
			}
		}
		else if (arg == "-t" || arg == "--threads")
		{
			++i;
			if (i < argc)
			{
				number_of_threads = stoi(argv[i]);
			}
			else
			{
				showUsage(argv[0]);
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			showUsage(argv[0]);
			exit(EXIT_FAILURE);
		}
		++i;
	}
}

void loadMeshes(const std::string& aFileName, vector<TriangleMesh>& aMeshSet){
	// Create an instance of the importer class
	Assimp::Importer importer;

	// read the file with some example postprocessing
	// Usually - if speed is no the most important aspect for you
	// - you'll probably to request more postprocessing than we do in this example
	// I dont need to change anything about the postprocessing at the moment
	const aiScene* scene = importer.ReadFile(aFileName, 
		aiProcess_CalcTangentSpace 		|
		aiProcess_Triangulate			| // Triangulates all the faces in a mesh
		aiProcess_JoinIdenticalVertices	|
		aiProcess_SortByPType);

	// check if importer failed
	if(!scene){
		// Stringstream makes it easier to build the error message this way
		std::stringstream error_message;
		error_message << importer.GetErrorString() << ", in File" << __FILE__ << 
				", in Function" << __FUNCTION__ <<
				", at Line" << __LINE__;
		throw std::runtime_error(error_message.str());
	}

	// Now we can access the file's contents
	if(scene->HasMeshes()){
		
		aMeshSet.clear(); // Clearing vector of any other mesh sets
		
		for(int mesh_id = 0; mesh_id < scene->mNumMeshes; ++mesh_id){
			aiMesh* p_mesh = scene->mMeshes[mesh_id];
			TriangleMesh mesh;

			// Checks it is a triangle mesh and load vertices
			if(p_mesh->mPrimitiveTypes == aiPrimitiveType_TRIANGLE){
				// Load materials here if needed
				// Don't need materials in an x-ray, at least not graphics materials

				// Stored using floats instead of Vec3 because its more efficient?
				// Loading a face at a time, rather than a vertex at a time

				std::vector<float> p_vertices;
				for(unsigned int vertex_id = 0; vertex_id < p_mesh->mNumVertices; ++vertex_id){
					// Appends the x,y,z of a vertex to the vertices vector
					p_vertices.push_back(p_mesh->mVertices[vertex_id].x);
					p_vertices.push_back(p_mesh->mVertices[vertex_id].y);
					p_vertices.push_back(p_mesh->mVertices[vertex_id].z);					
				}

				// Load indices
				std::vector<unsigned int> p_index_set;
				for(unsigned int index_id = 0; index_id < p_mesh->mNumFaces; ++index_id){
					// Face must have 3 indices, since we're working with triangles
					if(p_mesh->mFaces[index_id].mNumIndices == 3){
						p_index_set.push_back(p_mesh->mFaces[index_id].mIndices[0]);
						p_index_set.push_back(p_mesh->mFaces[index_id].mIndices[1]);
						p_index_set.push_back(p_mesh->mFaces[index_id].mIndices[2]);
					}
				}
				mesh.setGeometry(p_vertices, p_index_set);
			}
			aMeshSet.push_back(mesh);
		}
	}
}

// Changes i've made to this
// Removed material setting code - dont need it in this version
// Read and understood how meshes are loaded and how geometry is built - this will help when the time comes in Vulkan
