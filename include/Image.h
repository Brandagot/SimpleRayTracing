#ifndef __Image_h
#define __Image_h


/**
 ********************************************************************************
 *
 *   @file       Image.h
 *
 *   @brief      Class to manage 2D greyscale images in float.
 *
 *   @version    1.0
 *
 *   @date       12/08/2021
 *
 *   @author     Dr Franck P. Vidal
 *
 ********************************************************************************
 */


//******************************************************************************
//  Include
//******************************************************************************
#include <string>
#include <vector>


//==============================================================================
/**
*   @class  Image
*   @brief  Image is a class to manage 2D RGB images.
*/
//==============================================================================
class Image
//------------------------------------------------------------------------------
{
//******************************************************************************
public:
    //--------------------------------------------------------------------------
    /// Default Constructor
    //--------------------------------------------------------------------------
    Image();


    //--------------------------------------------------------------------------
    /// Copy constructor
    /*
    *   @param anImage  the image to copy
    */
    //--------------------------------------------------------------------------
    Image(const Image& anImage);


    //--------------------------------------------------------------------------
    /// Create a blank image
    /*
     *   @param aWidth   the image width (in number of pixels)
     *   @param aHeight  the image height (in number of pixels)
     */
    //--------------------------------------------------------------------------
    Image(unsigned int aWidth,
    		unsigned int aHeight,
    		float aPixelValue = 0);


    //--------------------------------------------------------------------------
    /// Destructor
    //--------------------------------------------------------------------------
    ~Image();


    //--------------------------------------------------------------------------
    /// Copy operator
    /*
    *   @param anImage  the image to copy
    */
    //--------------------------------------------------------------------------
    Image& operator=(const Image& anImage);


    //--------------------------------------------------------------------------
    /// Release memory
    //--------------------------------------------------------------------------
    void destroy();


    //--------------------------------------------------------------------------
    /// Save the current image into a JPEG file
    /*
     *   @param aFileName    the name of the JPEG file
     */
    //--------------------------------------------------------------------------
    void saveJPEGFile(const char* aFileName,
        float vmin = 0.0,
        float vmax = 1.0);


    //--------------------------------------------------------------------------
    /// Save the current image into a JPEG file
    /*
     *   @param aFileName    the name of the JPEG file
     */
    //--------------------------------------------------------------------------
    void saveJPEGFile(const std::string& aFileName,
        float vmin = 0.0,
        float vmax = 1.0);


    //--------------------------------------------------------------------------
    /// Save the current image into a JPEG file
    /*
     *   @param aFileName    the name of the JPEG file
     */
    //--------------------------------------------------------------------------
    void saveTGAFile(const char* aFileName,
        float vmin = 0.0,
        float vmax = 1.0);


    //--------------------------------------------------------------------------
    /// Save the current image into a JPEG file
    /*
     *   @param aFileName    the name of the JPEG file
     */
    //--------------------------------------------------------------------------
    void saveTGAFile(const std::string& aFileName,
        float vmin = 0.0,
        float vmax = 1.0);


    //--------------------------------------------------------------------------
    /// Save the current image into a text file
    /*
     *   @param aFileName    the name of the JPEG file
     */
    //--------------------------------------------------------------------------
    void saveTextFile(const std::string& aFileName);


    //--------------------------------------------------------------------------
    /// Accessor on the image size (in number of pixels)
    /*
     *   @param aWidth   the image width
     *   @param aHeight   the image height
     */
    //--------------------------------------------------------------------------
    void getSize(unsigned int& aWidth, unsigned int& aHeight) const;


    //--------------------------------------------------------------------------
    /// Accessor on the image width (in number of pixels)
    /*
     *   @return   the image width
     */
    //--------------------------------------------------------------------------
    unsigned int getWidth() const;


    //--------------------------------------------------------------------------
    /// Accessor on the image height (in number of pixels)
    /*
     *   @return   the image height
     */
    //--------------------------------------------------------------------------
    unsigned int getHeight() const;


    //--------------------------------------------------------------------------
    /// Accessor on the raw pixel values
    /*
    *   @return the raw pixel values
    */
    //--------------------------------------------------------------------------
    float* getData() const;


    void setPixel(unsigned int i,
                  unsigned int j,
                  float aPixelValue);


    void getPixel(unsigned int i,
                  unsigned int j,
                  float& aPixelValue) const;
                  
                  
    std::vector<unsigned char> applyLUT(float vmin, float vmax);


//******************************************************************************
protected:
    //--------------------------------------------------------------------------
    /// Allocate memory
    /*
    *   @param aWidth   the image width (in number of pixels)
    *   @param aHeight  the image height (in number of pixels)
    */
    //--------------------------------------------------------------------------
    void setSize(unsigned int aWidth, unsigned int aHeight,
    		float aPixelValue = 0);


    /// The pixel data
    float*         m_p_pixel_data;


    /// The image width (in number of pixels)
    unsigned int   m_width;


    /// The image height (in number of pixels)
    unsigned int   m_height;
};


#include "Image.inl"


#endif
