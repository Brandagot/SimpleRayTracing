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
#include <glm/glm.hpp>

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

// Dont need typedef in c++
// Decided to cache the bbox corners for each thread to access later
struct RayTracerInfo
{
	Vec3 detector_position;
	Vec3 origin;
	Vec3 up;
	Vec3 right;
	Light light;
	Vec3 upper_bbox_corner;
	Vec3 lower_bbox_corner;
	Vec3 range;
};

class PThreadData{
	// Thread information
	// Using a class constructor to initialize the parameters to we can use references

	private:
		PThreadData(); // Private so we can't call the default instructor.
	public:
		// Initialising the values of this data class - Required when using references
		PThreadData(Image& anOutputImage,
					const vector<TriangleMesh>& aTriangleMeshSet,
					const RayTracerInfo& rayTracerInfo):
					m_output_image(anOutputImage),
					m_triangle_mesh_set(aTriangleMeshSet),
					m_ray_tracer_info(rayTracerInfo),
					m_thread_int_id(-1),
					m_start_pixel(-1),
					m_end_pixel(-1){}

		pthread_t m_thread_id;
		unsigned int m_thread_int_id;

		// Reference to variables needed in the render loop.
		Image& m_output_image;
		const vector<TriangleMesh>& m_triangle_mesh_set;
		const RayTracerInfo& m_ray_tracer_info;

		// Start and end point for data assigned to thread - Thread will work on start pixel through to end pixel
		unsigned int m_start_pixel;
		unsigned int m_end_pixel;
};

// How to define a matrix
// Rotation matrix, 90 degrees
glm::mat4x4 rotateMat(
	std::cos(0.785398),	0,std::sin(0.785398),0,
	0,					1,		0,		   0,
	-std::sin(0.785398),	0,std::cos(0.785398),0,
	0,					0,		0,  	   1
);

// Base transform matrix
glm::mat4x4 transformMatrix(
	1,0,0,0,
	0,1,0,0,
	0,0,1,0,
	0,0,0,1
);

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
	unsigned int& number_of_threads
);

// Retrieves a .ply mesh from file and adds it to the set
void loadMeshes(
	const std::string& aFileName,
	vector<TriangleMesh>& aMeshSet,
	glm::mat4x4 tranformMatrix
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

RayTracerInfo initialiseRayTracing(
	vector<TriangleMesh>& aMeshSet,
	const Vec3& anUpperBBoxCorner,
	const Vec3& aLowerBBoxCorner,
	const unsigned int image_height,
	const unsigned int image_width,
	Image& output_image,
	float lut
);

void pthreadWorkLoadAllocation(
	vector<PThreadData>& p_thread_data,
	unsigned int number_of_threads,
	Image& output_image
);

template <typename T> int signum(T val);

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

std::vector<float> L_buffer; // dont need to lock this, threads only access data based on pixels they have been assigned 



int main(int argc, char** argv){
	try{
		// defualt output file
		string output_file_name = "test.jpg";

		// default image size - changed later if passed processCMD
		unsigned int image_width = g_default_image_width;
		unsigned int image_height= g_default_image_height;

		// Default background color is set to grey - might try size the viewplane to match the entire view
		float lut = 0.0f;
		
		// Default number of threads
		unsigned int number_of_threads = 1;

		glm::mat4x4 transformMatrix(
			1,0,0,0,
			0,1,0,0,
			0,0,1,0,
			0,0,0,1
		);
		// transformMatrix = transformMatrix * rotateMat;

		// updates above defaults if passed command args
		processCmd(argc, argv, 
			output_file_name,
			image_width, image_height,
			number_of_threads);

		// Loading polygon meshes
		cout << "Loading polygon meshes" << endl << endl;
		
		vector<TriangleMesh> p_mesh_set;
		loadMeshes("./dragon_monkey.ply", p_mesh_set, transformMatrix);

		cout <<  "Retreiving scenes bbox" << endl;
 
		// Which corner in a 3D object?
		Vec3 lower_bbox_corner;
		Vec3 upper_bbox_corner;

		getBBox(p_mesh_set, upper_bbox_corner, lower_bbox_corner);

		// getting bounding box of the object
		// dont really get it. Gets the direction between the start anbd the 
		// Initialise ray-tracer properties
		Image output_image;
		RayTracerInfo rayTracerInfo = initialiseRayTracing(p_mesh_set, upper_bbox_corner, lower_bbox_corner, image_height, image_width, output_image, lut);

		// Allocate work for Pthreads
		// Need to change PThread data to accept RayTracerInfo
		vector<PThreadData> p_thread_data(number_of_threads, PThreadData(output_image, p_mesh_set, rayTracerInfo));
		pthreadWorkLoadAllocation(p_thread_data, number_of_threads, output_image);
		
		cout << "Creating and Executing threads" << endl << endl;

    	L_buffer = std::vector<float>(output_image.getWidth() * output_image.getHeight(), 80.000f);

		// Create threads, executes callback on each thread
		for(int i = 0; i < number_of_threads; i++){
			// Need to create the render loop callback - the render loop is going to be reworked
			pthread_create(&p_thread_data[i].m_thread_id, NULL, renderLoopCallBack, &p_thread_data[i]);
		}

		// Join threads, waiting for them all to finish
		for(int i = 0; i < number_of_threads; i++){
			pthread_join(p_thread_data[i].m_thread_id, NULL);
		}

		// image output loop
		// Placing error correction here
		for(int pixel = 0; pixel < output_image.getWidth() * output_image.getHeight(); pixel++){
			// x,y of pixel in the image
			unsigned int row = pixel / output_image.getWidth();
			unsigned int col = pixel % output_image.getWidth();

			float photonOut = L_buffer[row * output_image.getWidth() + col];

			// Need to get the average of all non-zero numbers
			// This might be a rough method for now - probably super ineffecient. Only does a cross, not square around the pixel -- this works in a + pattern for now, need to change this to be square. Its not a true average
			if(photonOut == -1){
				uint maxDepth = 4;
				std::vector<float> averagingValues;

				float value1 = 0;
				for(int i = 1; i <= maxDepth; i++){
					if((row * output_image.getWidth() + (col + i)) > L_buffer.size()){
						break;
					}

					if(L_buffer[row * output_image.getWidth() + (col + i)] != -1){
						value1 = L_buffer[row * output_image.getWidth() + (col + i)];
						break;
					}
				}
				if(value1 != 0) averagingValues.push_back(value1);

				float value2 = 0;
				for(int i = 1; i <= maxDepth; i++){
					if(((row - i) * output_image.getWidth() + col) > L_buffer.size()){
						break;
					}

					if(L_buffer[(row - i) * output_image.getWidth() + col] != -1){
						value2 = L_buffer[(row - i) * output_image.getWidth() + col];
						break;
					}
				}			
				if(value2 != 0) averagingValues.push_back(value2);
				
				float value3 = 0;
				for(int i = 1; i <= maxDepth; i++){
					if((row * output_image.getWidth() + (col - i)) > L_buffer.size()){
						break;
					}

					if(L_buffer[row * output_image.getWidth() + (col - i)] != -1){
						value3 = L_buffer[row * output_image.getWidth() + (col - i)];
						break;
					}
				}
				if(value3 != 0) averagingValues.push_back(value3);
				
				float value4 = 0;
				for(int i = 1; i <= maxDepth; i++){
					if(((row + i) * output_image.getWidth() + col) > L_buffer.size()){
						break;
					}

					if(L_buffer[(row + i) * output_image.getWidth() + col] != -1){
						value4 = L_buffer[(row + i) * output_image.getWidth() + col];
						break;
					}
				}
				if(value4 != 0) averagingValues.push_back(value4);

				float averageSum = 0;
				for(float value : averagingValues){
					averageSum += value;
				}

				photonOut = averageSum / averagingValues.size();

			}

			output_image.setPixel(col, row, photonOut);
		}

		output_image.saveTextFile(output_file_name); // TODO: Read JPEG exporter and check if gamma corrected?

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
		else if (arg == "-f" || arg == "--filename")
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

// If i want to transform the stuff, i know to apply an transform metrix to the stuff here right?
void loadMeshes(const std::string& aFileName, vector<TriangleMesh>& aMeshSet, glm::mat4x4 transformMatrix){
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
					// Apply transform to vertices (if there is one)
					p_vertices.push_back((transformMatrix[0][0] * p_mesh->mVertices[vertex_id].x) + (transformMatrix[0][1] * p_mesh->mVertices[vertex_id].y) + (transformMatrix[0][2] * p_mesh->mVertices[vertex_id].z));			
					p_vertices.push_back((transformMatrix[1][0] * p_mesh->mVertices[vertex_id].x) + (transformMatrix[1][1] * p_mesh->mVertices[vertex_id].y) + (transformMatrix[1][2] * p_mesh->mVertices[vertex_id].z));
					p_vertices.push_back((transformMatrix[2][0] * p_mesh->mVertices[vertex_id].x) + (transformMatrix[2][1] * p_mesh->mVertices[vertex_id].y) + (transformMatrix[2][2] * p_mesh->mVertices[vertex_id].z));
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

// Gets the boudning box
void getBBox(
	const vector<TriangleMesh>& aMeshSet,
	Vec3& anUpperBBoxCorner,
	Vec3& aLowerBBoxCorner
){
	float inf = std::numeric_limits<float>::infinity();

	aLowerBBoxCorner = Vec3(inf, inf, inf);
	anUpperBBoxCorner = Vec3(-inf, -inf, -inf);

	for(std::vector<TriangleMesh>::const_iterator ite = aMeshSet.begin(); ite != aMeshSet.end(); ++ite){
		Vec3 mesh_lower_bbox_corner = ite->getLowerBBoxCorner();
		Vec3 mesh_upper_bbox_corner = ite->getUpperBBoxCorner();

		// Not entirely sure here
		aLowerBBoxCorner[0] = std::min(aLowerBBoxCorner[0], mesh_lower_bbox_corner[0]);
		aLowerBBoxCorner[1] = std::min(aLowerBBoxCorner[1], mesh_lower_bbox_corner[1]);
		aLowerBBoxCorner[2] = std::min(aLowerBBoxCorner[2], mesh_lower_bbox_corner[2]);

		anUpperBBoxCorner[0] = std::max(anUpperBBoxCorner[0], mesh_upper_bbox_corner[0]);
		anUpperBBoxCorner[1] = std::max(anUpperBBoxCorner[1], mesh_upper_bbox_corner[1]);
		anUpperBBoxCorner[2] = std::max(anUpperBBoxCorner[2], mesh_upper_bbox_corner[2]);
	}
}

// Might be too many things going on in here, could pass it a struct instead?
RayTracerInfo initialiseRayTracing(	
	vector<TriangleMesh>& aMeshSet,
	const Vec3& anUpperBBoxCorner,
	const Vec3& aLowerBBoxCorner,
	const unsigned int image_height,
	const unsigned int image_width,
	Image& output_image,
	float lut){
	
	Vec3 range = anUpperBBoxCorner - aLowerBBoxCorner;
	Vec3 bbox_centre = aLowerBBoxCorner + range / 2.0;


	// Helps me visual the difference in each corner. Which is top, right and forward, etc
	#ifndef NDEBUG
		cout << bbox_centre.getX() << ", " << bbox_centre.getY() << ", " << bbox_centre.getZ() << endl;
		cout << aLowerBBoxCorner.getX() << ", " << aLowerBBoxCorner.getY() << ", " << aLowerBBoxCorner.getZ() << endl;
		cout << anUpperBBoxCorner.getX() << ", " << anUpperBBoxCorner.getY() << ", " << anUpperBBoxCorner.getZ() << endl;
	#endif

	float diagonal = range.getLength(); // Diagonal distance between the upper and lower corner

	Vec3 up(0.0, 0.0, -1.0);

	Vec3 origin(bbox_centre - Vec3(diagonal * 1, 0, 0)); // originall * 1
	Vec3 detector_position(bbox_centre + Vec3(diagonal * 0.6, 0, 0)); // was 0.6

	Vec3 direction((detector_position - origin));
	direction.normalize();

	output_image = Image(image_width, image_height, lut);

	Vec3 light_position = origin;
	Vec3 light_direction = bbox_centre - light_position;
	light_direction.normalise();
	Light light(g_white, light_direction, light_position);

	direction.normalise();
	Vec3 right(direction.crossProduct(up));	

	// Creates mesh that will go behind the scene
	// Dont need this really
	// aMeshSet.push_back(createBackground(anUpperBBoxCorner, aLowerBBoxCorner));

	RayTracerInfo info {
		detector_position,
		origin,
		up,
		right,
		light,
		anUpperBBoxCorner,
		aLowerBBoxCorner,
		range
	};

	return info;
}

// Might be obsolete now that we're going to create a view plane instead
TriangleMesh createBackground(
	const Vec3& anUpperBBoxCorner,
	const Vec3& aLowerBBoxCorner
){
	Vec3 range = anUpperBBoxCorner - aLowerBBoxCorner;

	// four vertex positions? for each corner of the detector / background?
	std::vector<float> vertices = {
		anUpperBBoxCorner[0] + range[0] * 0.1f, aLowerBBoxCorner[1] - range[1] * 0.5f, aLowerBBoxCorner[2] - range[2] * 0.5f,
		anUpperBBoxCorner[0] + range[0] * 0.1f, anUpperBBoxCorner[1] + range[1] * 0.5f, aLowerBBoxCorner[2] - range[2] * 0.5f,
		anUpperBBoxCorner[0] + range[0] * 0.1f, anUpperBBoxCorner[1] + range[1] * 0.5f, anUpperBBoxCorner[2] + range[2] * 0.5f,
		anUpperBBoxCorner[0] + range[0] * 0.1f, aLowerBBoxCorner[1] - range[1] * 0.5f, anUpperBBoxCorner[2] + range[2] * 0.5f,
	};

	std::vector<unsigned int> indices = {
		0, 1, 2,
		0, 2, 3
	};

	TriangleMesh background_mesh(vertices, indices);
	return (background_mesh);
}

//Only ints and smaller objects should be passed by value, because it's cheaper to copy them than to take the dereferencing hit within the function.
void pthreadWorkLoadAllocation(
	vector<PThreadData>& p_thread_data,
	unsigned int number_of_threads,
	Image& output_image
){
	cout << "Allocating PThread data" << endl << endl;
	int last_element = -1;
	unsigned int total_pixels = output_image.getWidth() * output_image.getHeight();
	unsigned int pixels_per_thread = total_pixels / number_of_threads;
	unsigned int remainder = total_pixels % number_of_threads;

	cout << "Number of cells per thread (1D array) " << pixels_per_thread << endl << endl;

	for(int i = 0; i < number_of_threads; i++){
		// Setting thread ID, and pixel range to perform work on
		p_thread_data[i].m_thread_int_id = i;
		p_thread_data[i].m_start_pixel = ++last_element;
		p_thread_data[i].m_end_pixel = last_element + pixels_per_thread - 1;

		// Adding a pixel from the remainders to the current thread
		if(remainder > 0){
			p_thread_data[i].m_end_pixel++;
			--remainder;
		}

		last_element = p_thread_data[i].m_end_pixel;
		
		cout << "Thread: " << p_thread_data[i].m_thread_int_id << endl;
	}

	cout << endl;
}

// typesafe signum function
// https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
template <typename T> int signum(T val){
	return(T(0) < val) - (val < T(0));
}

void* renderLoopCallBack(void* apData){
	PThreadData* p_thread_data = static_cast<PThreadData*>(apData);
	// cout << p_thread_data->m_ray_tracer_info.upper_bbox_corner << endl;
	// cout << p_thread_data->m_ray_tracer_info.range << endl;
    // maybe a local L buffer instead?

	float res1 = p_thread_data->m_ray_tracer_info.range[2] / p_thread_data->m_output_image.getWidth();
	float res2 = p_thread_data->m_ray_tracer_info.range[1] / p_thread_data->m_output_image.getHeight();
	float pixel_spacing[] = {2 * std::max(res1, res2), 2 * std::max(res1, res2)}; // ?

	// Process everyrow
	float inf  = std::numeric_limits<float>::infinity();

	// Need to replace this with an L-buffer 
    // This needs to go outside of this callback, to be shared with other threads

	// Process all allocated pixels
	for(int pixel = p_thread_data->m_start_pixel; pixel <= p_thread_data->m_end_pixel; ++pixel){
		
		// x,y of pixel in the image
		unsigned int row = pixel / p_thread_data->m_output_image.getWidth();
		unsigned int col = pixel % p_thread_data->m_output_image.getWidth();
		
		float v_offset = pixel_spacing[1] * (0.5 + row - p_thread_data->m_output_image.getHeight() / 2.0);
		float u_offset = pixel_spacing[0] * (0.5 + col - p_thread_data->m_output_image.getWidth() / 2.0);

		// Initialise the ray direction for this pixel
		Vec3 direction = p_thread_data->m_ray_tracer_info.detector_position + p_thread_data->m_ray_tracer_info.up * v_offset + p_thread_data->m_ray_tracer_info.right * u_offset - p_thread_data->m_ray_tracer_info.origin;
		direction.normalise();
		Ray ray(p_thread_data->m_ray_tracer_info.origin, direction);

		const TriangleMesh* p_intersected_object = 0;
		const Triangle* p_intersected_triangle = 0;

		// Process each mesh
		for(std::vector<TriangleMesh>::const_iterator mesh_ite = p_thread_data->m_triangle_mesh_set.begin();
				mesh_ite != p_thread_data->m_triangle_mesh_set.end();
				++mesh_ite){
			
			if(L_buffer[row * p_thread_data->m_output_image.getWidth() + col] == -1) break; // If the pixel has been flagged before, don't run it for the next mesh
			
			// The ray intersects with the mesh's bbox	
			if(mesh_ite->intersectBBox(ray)){
				
				// Process all triangles of the mesh
				float distance = 0.0f;
				int signSum = 0; // Used for validation
				for(unsigned int triangle_id = 0; triangle_id < mesh_ite->getNumberOfTriangles(); ++triangle_id){
					const Triangle& triangle = mesh_ite->getTriangle(triangle_id);
					
					// Retrieve intersect if any
					float t;
					bool intersect = ray.intersect(triangle, t);

					// checks if the intersect is with the dragon
					if(intersect && &*mesh_ite == &p_thread_data->m_triangle_mesh_set[0] && t > 0.0000001){
                        // float dotProduct = direction.dotProduct(triangle.getNormal());
						// distance += (dotProduct < 0) ? -t : t;			
						int sign = signum(direction.dotProduct(triangle.getNormal()));
						distance += (sign * t);
						signSum += sign; // Method of finding artefacts taken from francks paper: http://www.fpvidal.net/research/pdf/Vidal2016ComputMedImagingGraph.pdf 
					}
				}
				// Buffer save code here
				// mju of bone: 0.3971f - go back and find which of Franck's papers I used for this http://www.fpvidal.net/research/pdf/Vidal2016ComputMedImagingGraph.pdf
				float attenuationCoefficient = 0;
				if(&(*mesh_ite) == &(p_thread_data->m_triangle_mesh_set[0])){
					attenuationCoefficient = 0.1037f; // Soft tissue
				} else {
					attenuationCoefficient = 0.3971f; // Bone
				}
		
				if(signSum != 0){
					L_buffer[row * p_thread_data->m_output_image.getWidth() + col] = -1;
				} else{
					L_buffer[row * p_thread_data->m_output_image.getWidth() + col] = L_buffer[row * p_thread_data->m_output_image.getWidth() + col] * exp(-(attenuationCoefficient * (distance * 0.1))); // Converting distance to cm
				}
			}
		}
	}
}


// Changes i've made to this
// Removed material setting code - dont need it in this version
// Read and understood how meshes are loaded and how geometry is built - this will help when the time comes in Vulkan
// Used debud printing to visualise stuff. Hard to picture it in my brain
// They work out the bounding box by setting a point with the biggest or lower x, y, z value out of all the points in all the triangles
// Store ray-tracing info in a strcut for threads to access
// removed the additional shadowray being shot - no longer needed it
// resulting in a removal of the seconds loop - making things more efficient

// a 2048 x 2048 image with 6 threads took 2238.94 seconds to execute / 37m18.992s overall using the old method
// takes 23 on this new version, and uses more of the image
// Going to see if it is considerably shorter using the Lbuffer alg or not
// Now implemented to Lbuffer algo, not sure how efficient it is, seeing as each thread is creating an L-buffer - using a shit ton of memory
// Benefits are negligible on CPU - L buffer coomes in when GPU acceleration is available