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

//******************************************************************************
//  Function declarations
//******************************************************************************
void showUsage(const std::string& aProgramName);

void processCmd(int argc, char** argv,
				string& aFileName,
				unsigned int& aWidth, unsigned int& aHeight,
				unsigned char& r, unsigned char& g, unsigned char& b, string& inputFileName);

Vec3 applyShading(const Light& aLight,
				  const Material& aMaterial,
				  const Vec3& aNormalVector,
				  const Vec3& aPosition,
				  const Vec3& aViewPosition);

void loadMeshes(const std::string& aFileName,
				vector<TriangleMesh>& aMeshSet);

TriangleMesh createBackground(const Vec3& anUpperBBoxCorner,
							  const Vec3& aLowerBBoxCorner);

void getBBox(const vector<TriangleMesh>& aMeshSet,
			 Vec3& anUpperBBoxCorner,
			 Vec3& aLowerBBoxCorner);

RayTracerInfo initialiseRayTracing(
	vector<TriangleMesh>& aMeshSet,
	const Vec3& anUpperBBoxCorner,
	const Vec3& aLowerBBoxCorner,
	const unsigned int image_height,
	const unsigned int image_width,
	Image& output_image,
	float lut
);

void renderLoop(Image& anOutputImage,
				const vector<TriangleMesh>& aTriangleMeshSet,
				RayTracerInfo& rayTracerInfo);


//******************************************************************************
//  Constant global variables
//******************************************************************************
const unsigned int g_default_image_width  = 2048;
const unsigned int g_default_image_height = 2048;

const Vec3 g_black(0, 0, 0);
const Vec3 g_white(1, 1, 1);

const Vec3 g_red(1, 0, 0);
const Vec3 g_green(0, 1, 0);
const Vec3 g_blue(0, 0, 1);

const Vec3 g_background_colour = g_black;
RayTracerInfo rayTracerInfo;


//-----------------------------
int main(int argc, char** argv)
//-----------------------------
{
	try
	{
		// output file
		string output_file_name = "test.jpg";
		
		// Update the image size if needed
		unsigned int image_width = g_default_image_width;
		unsigned int image_height = g_default_image_height;

		string input_mesh_file = "./dragon.ply"; // Default file is dragon

		// Update the background colour if needed
		unsigned char r = 128;
		unsigned char g = 128;
		unsigned char b = 128;

		processCmd(argc, argv,
				   output_file_name,
				   image_width, image_height,
				   r, g, b, input_mesh_file);
						   
		// Start timer for loadMeshes
		cout << "Loading polygon meshes... " << endl;
		auto start = chrono::high_resolution_clock::now();
		
		// Load the polygon meshes
		vector<TriangleMesh> p_mesh_set;
		loadMeshes(input_mesh_file, p_mesh_set);

		// Stop and output time taken for loadMeshes
		auto end = chrono::high_resolution_clock::now();
		cout << "Loading meshes took: " << chrono::duration<double>(end - start).count() << " seconds" << endl << endl;


		// Start timer for getting the scene box
		cout << "Getting scene bbox... " << endl;
		start = chrono::high_resolution_clock::now();
		
		// Get the scene's bbox
		Vec3 lower_bbox_corner;
		Vec3 upper_bbox_corner;

		getBBox(p_mesh_set, upper_bbox_corner, lower_bbox_corner);

		// Stop and output time taken for getBBox
		float lut = 0.0f;

		Image output_image(image_width, image_height, lut);
		rayTracerInfo = initialiseRayTracing(p_mesh_set, upper_bbox_corner, lower_bbox_corner, image_height, image_width, output_image, lut);

		
		// Rendering loop
		renderLoop(output_image, p_mesh_set, rayTracerInfo);

		// Save the image
		output_image.saveTextFile(output_file_name);
	}
	// Catch exceptions and error messages
	catch (const std::exception& e)
	{
		std::cerr << "ERROR: " << e.what() << std::endl;
		return 1;
	}
	catch (const std::string& e)
	{
		std::cerr << "ERROR: " << e << std::endl;
		return 2;
	}
	catch (const char* e)
	{
		std::cerr << "ERROR: " << e << std::endl;
		return 3;
	}

	return 0;
}


//******************************************************************************
//  Function definitions
//******************************************************************************


//---------------------------------------------
void showUsage(const std::string& aProgramName)
//---------------------------------------------
{
	std::cerr << "Usage: " << aProgramName << " <option(s)>" << endl <<
		"Options:" << endl <<
		"\t-h,--help\t\t\tShow this help message" << endl <<
		"\t-s,--size IMG_WIDTH IMG_HEIGHT\tSpecify the image size in number of pixels (default values: 2048 2048)" << endl << 
		"\t-b,--background R G B\t\tSpecify the background colour in RGB, acceptable values are between 0 and 255 (inclusive) (default values: 128 128 128)" << endl << 
		"\t-j,--jpeg FILENAME\t\tName of the JPEG file (default value: test.jpg)" << endl << 
		std::endl;
}


//-------------------------------------------------------------------
void processCmd(int argc, char** argv,
				string& aFileName,
				unsigned int& aWidth, unsigned int& aHeight,
				unsigned char& r, unsigned char& g, unsigned char& b, string& inputFileName)
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
		else if (arg == "-f" || arg == "--filename")
		{                
			++i;
			if (i < argc)
			{
				aFileName = ("./out/" + string(argv[i]));
			}
			else
			{
				showUsage(argv[0]);
				exit(EXIT_FAILURE);
			}
		}
		else if (arg == "-i" || arg == "--input")
		{
			++i;
			if (i < argc)
			{
				inputFileName = argv[i];
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


//------------------------------------------
Vec3 applyShading(const Light& aLight,
				  const Material& aMaterial,
				  const Vec3& aNormalVector,
				  const Vec3& aPosition,
				  const Vec3& aViewPosition)
//------------------------------------------
{
	Vec3 ambient, diffuse, specular;

	// ambient
	ambient = aLight.getColour() * aMaterial.getAmbient();

	// diffuse
	Vec3 lightDir = (aLight.getPosition() - aPosition);
	lightDir.normalize();
	float diff = std::max(std::abs(aNormalVector.dotProduct(lightDir)), 0.0f);
	diffuse = aLight.getColour() * (diff * aMaterial.getDiffuse());

	// specular
	Vec3 viewDir(aViewPosition - aPosition);
	viewDir.normalize();

	Vec3 reflectDir = reflect(-viewDir, aNormalVector);
	float spec = std::pow(std::max(dot(viewDir, reflectDir), 0.0f), aMaterial.getShininess());
	specular = aLight.getColour() * (spec * aMaterial.getSpecular());

	return ambient + diffuse + specular;
}


//---------------------------------------------
void loadMeshes(const std::string& aFileName,
								vector<TriangleMesh>& aMeshSet)
//-----------------------------==--------------
{
	// Create an instance of the Importer class
	Assimp::Importer importer;

	// And have it read the given file with some example postprocessing
	// Usually - if speed is not the most important aspect for you - you'll
	// probably to request more postprocessing than we do in this example.
	const aiScene* scene = importer.ReadFile( aFileName,
			aiProcess_CalcTangentSpace       |
			aiProcess_Triangulate            |
			aiProcess_JoinIdenticalVertices  |
			aiProcess_SortByPType);

	// If the import failed, report it
	if( !scene)
	{
		std::stringstream error_message;
		error_message << importer.GetErrorString() << ", in File " << __FILE__ <<
				", in Function " << __FUNCTION__ <<
				", at Line " << __LINE__;

		throw std::runtime_error(error_message.str());
	}

	// Now we can access the file's contents.
	if (scene->HasMeshes())
	{
		aMeshSet.clear();

		for (int mesh_id = 0; mesh_id < scene->mNumMeshes; ++mesh_id)
		{
			aiMesh* p_mesh = scene->mMeshes[mesh_id];
			TriangleMesh mesh;

			// This is a triangle mesh
			if (p_mesh->mPrimitiveTypes == aiPrimitiveType_TRIANGLE)
			{
				aiMaterial* p_mat = scene->mMaterials[p_mesh->mMaterialIndex];

				aiColor3D ambient, diffuse, specular;
				float shininess;

				p_mat->Get(AI_MATKEY_COLOR_AMBIENT, ambient);
				p_mat->Get(AI_MATKEY_COLOR_DIFFUSE, diffuse);
				p_mat->Get(AI_MATKEY_COLOR_SPECULAR, specular);
				p_mat->Get(AI_MATKEY_SHININESS, shininess);

				Material material;
				material.setAmbient(Vec3(ambient.r, ambient.g, ambient.b));
				material.setDiffuse(Vec3(diffuse.r, diffuse.g, diffuse.b));
				material.setSpecular(Vec3(specular.r, specular.g, specular.b));
				material.setShininess(shininess);

				mesh.setMaterial(material);

				// Load the vertices
				std::vector<float> p_vertices;
				for (unsigned int vertex_id = 0; vertex_id < p_mesh->mNumVertices; ++vertex_id)
				{
					p_vertices.push_back(p_mesh->mVertices[vertex_id].x);
					p_vertices.push_back(p_mesh->mVertices[vertex_id].y);
					p_vertices.push_back(p_mesh->mVertices[vertex_id].z);
				}

				// Load indices
				std::vector<unsigned int> p_index_set;
				for (unsigned int index_id = 0; index_id < p_mesh->mNumFaces; ++index_id)
				{
					if (p_mesh->mFaces[index_id].mNumIndices == 3)
					{
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


//----------------------------------------------------------
TriangleMesh createBackground(const Vec3& anUpperBBoxCorner,
								const Vec3& aLowerBBoxCorner)
//----------------------------------------------------------
{
	Vec3 range = anUpperBBoxCorner - aLowerBBoxCorner;

	std::vector<float> vertices = {
		anUpperBBoxCorner[0] + range[0] * 0.1f, aLowerBBoxCorner[1] - range[1] * 0.5f, aLowerBBoxCorner[2] - range[2] * 0.5f,
		anUpperBBoxCorner[0] + range[0] * 0.1f, anUpperBBoxCorner[1] + range[1] * 0.5f, aLowerBBoxCorner[2] - range[2] * 0.5f,
		anUpperBBoxCorner[0] + range[0] * 0.1f, anUpperBBoxCorner[1] + range[1] * 0.5f, anUpperBBoxCorner[2] + range[2] * 0.5f,
		anUpperBBoxCorner[0] + range[0] * 0.1f, aLowerBBoxCorner[1] - range[1] * 0.5f, anUpperBBoxCorner[2] + range[2] * 0.5f,
	};
	
	std::vector<unsigned int> indices = {
		0, 1, 2,
		0, 2, 3,
	};

	TriangleMesh background_mesh(vertices, indices);
	return (background_mesh);
}


//------------------------------------------------
void getBBox(const vector<TriangleMesh>& aMeshSet,
						 Vec3& anUpperBBoxCorner,
						 Vec3& aLowerBBoxCorner)
//------------------------------------------------
{
	float inf = std::numeric_limits<float>::infinity();

	aLowerBBoxCorner = Vec3( inf,  inf,  inf);
	anUpperBBoxCorner = Vec3(-inf, -inf, -inf);

	for (std::vector<TriangleMesh>::const_iterator ite = aMeshSet.begin();
			ite != aMeshSet.end();
			++ite)
	{
		Vec3 mesh_lower_bbox_corner = ite->getLowerBBoxCorner();
		Vec3 mesh_upper_bbox_corner = ite->getUpperBBoxCorner();

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
	Vec3 detector_position(bbox_centre + Vec3(diagonal * 0.6, 0, 0));

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


//-------------------------------------------------------------
void renderLoop(Image& anOutputImage,
				  const vector<TriangleMesh>& aTriangleMeshSet,
				  RayTracerInfo& rayTracerInfo)
//-------------------------------------------------------------
{
	// Initialise some parameters
	Vec3 upper_bbox_corner;
	Vec3 lower_bbox_corner;
	getBBox(aTriangleMeshSet, upper_bbox_corner, lower_bbox_corner);

	// Initialise the ray-tracer properties
	Vec3 range = upper_bbox_corner - lower_bbox_corner;

	float res1 = range[2] / anOutputImage.getWidth();
	float res2 = range[1] / anOutputImage.getHeight();
	float pixel_spacing[] = {2 * std::max(res1, res2), 2 * std::max(res1, res2)};

	// Process every row
	float inf = std::numeric_limits<float>::infinity();

	std::vector<float> z_buffer(anOutputImage.getWidth() * anOutputImage.getHeight(), inf);

	// Process all allocated pixels
	for(int pixel = 0; pixel < anOutputImage.getWidth() * anOutputImage.getHeight(); ++pixel){
		
		// x,y of pixel in the image
		unsigned int row = pixel / anOutputImage.getWidth();
		unsigned int col = pixel % anOutputImage.getWidth();
		
		float v_offset = pixel_spacing[1] * (0.5 + row - anOutputImage.getHeight() / 2.0);
		float u_offset = pixel_spacing[0] * (0.5 + col - anOutputImage.getWidth() / 2.0);

		// Initialise the ray direction for this pixel
		Vec3 direction = rayTracerInfo.detector_position + rayTracerInfo.up * v_offset + rayTracerInfo.right * u_offset - rayTracerInfo.origin;
		direction.normalise();
		Ray ray(rayTracerInfo.origin, direction);

		const TriangleMesh* p_intersected_object = 0;
		const Triangle* p_intersected_triangle = 0;

		std::vector<float> intersect_points;
		intersect_points.reserve(10); // can change

		// Process each mesh
		for(std::vector<TriangleMesh>::const_iterator mesh_ite = aTriangleMeshSet.begin();
				mesh_ite != aTriangleMeshSet.end();
				++mesh_ite){
			
			// The ray intersects with the mesh's bbox	
			if(mesh_ite->intersectBBox(ray)){
				
				// Process all triangles of the mesh
				for(unsigned int triangle_id = 0; triangle_id < mesh_ite->getNumberOfTriangles(); ++triangle_id){
					const Triangle& triangle = mesh_ite->getTriangle(triangle_id);
					
					// Retrieve intersect if any
					float t;
					bool intersect = ray.intersect(triangle, t);

					// checks if the intersect is with the dragon
					// &*?? lol
					if(intersect && &*mesh_ite == &aTriangleMeshSet[0] && t > 0.0000001){
						intersect_points.push_back(t);
						p_intersected_object = &(*mesh_ite);
						p_intersected_triangle = &triangle;
					}
				}

				intersect_points.shrink_to_fit(); // this doesnt effect time too much
			}
		}
		// cout << intersect_points << endl;

		// Sets the distance values for all the pixels
		float distance = 0;
		if(p_intersected_object && p_intersected_triangle){

			if(intersect_points.size() % 2 == 0){
				sort(intersect_points.begin(), intersect_points.end());
				
				for(int i = 0; i < intersect_points.size(); i+=2){
					distance += intersect_points[i + 1] - intersect_points[i];
				}
			} else{
				cout << "Only one intersect on this ray" << endl;
				// Not sure how to deal with artifact
				// distance = z_buffer[row * anOutputImage.getWidth() + col - 1];
			}

			// some of these values be assigned, need to give them a default color
			// Do I even need a buffer? dont it all per pixel in one loop means probably not
			z_buffer[row * anOutputImage.getWidth() + col] = distance;
		}

		// Visualisation the pixel

		// Abitrary value to clamp distance
		// keeping the zeros because they represent the min values for the clamping

		distance = distance * 0.1;
		// cout << distance << endl;
		// distance = (distance - 0) / 50 * (1 - 0) + 0;
		// cout << distance << endl;

		// Divide distance by 10 to cm - the distance is considered to be in mm

		// 80 keV input
		// use 80.000 for kev output
		// 0 smallest output
		// I incident output

		// 0.3971 is the mju of bone - lin attenuation coefficient
		// add this to variables later
		float photonOut = 80.000f * exp(-(0.3971f * distance));

		anOutputImage.setPixel(col, row, photonOut);
	}
}

