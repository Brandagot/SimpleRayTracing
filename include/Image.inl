/**
 ********************************************************************************
 *
 *   @file       Image.inl
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
#include <stdexcept>
#include <sstream>


//------------------------
inline Image::Image():
//------------------------
        m_p_pixel_data(0),
        m_width(0),
        m_height(0)
//------------------------
{
}


//----------------------------------------
inline Image::Image(unsigned int aWidth,
                    unsigned int aHeight,
                    float aPixelValue):
//----------------------------------------
        m_p_pixel_data(0),
        m_width(0),
        m_height(0)
//----------------------------------------
{
    setSize(aWidth, aHeight, aPixelValue);
}


//--------------------
inline Image::~Image()
//--------------------
{
    // Release the memory
    destroy();
}


//--------------------------
inline void Image::destroy()
//--------------------------
{
    // Release the memory
    if (m_p_pixel_data)
    {
        delete [] m_p_pixel_data;
        m_p_pixel_data = 0;
    }

    // Reset parameters to their default values
    m_width  = 0;
    m_height = 0;
}


//-----------------------------------------------------------
inline void Image::saveJPEGFile(const std::string& aFileName,
                                float vmin,
                                float vmax)
//-----------------------------------------------------------
{
    saveJPEGFile(aFileName.data(), vmin, vmax);
}


//----------------------------------------------------------
inline void Image::saveTGAFile(const std::string& aFileName,
                               float vmin,
                               float vmax)
//----------------------------------------------------------
{
    saveTGAFile(aFileName.data(), vmin, vmax);
}


//-----------------------------------------------------------
inline void saveTextFile(const std::string& aFileName)
//-----------------------------------------------------------
{
    saveTextFile(aFileName.data());
}


//-----------------------------------------------------
inline void Image::getSize(unsigned int& aWidth,
                           unsigned int& aHeight) const
//-----------------------------------------------------
{
    aWidth  = m_width;
    aHeight = m_height;
}


//-----------------------------------------
inline unsigned int Image::getWidth() const
//-----------------------------------------
{
    return m_width;
}


//------------------------------------------
inline unsigned int Image::getHeight() const
//------------------------------------------
{
    return m_height;
}


//----------------------------------
inline float* Image::getData() const
//----------------------------------
{
    return (m_p_pixel_data);
}


//--------------------------------------------
inline void Image::setPixel(unsigned int i,
                            unsigned int j,
                            float aPixelValue)
//--------------------------------------------
{
    // The 2D index is valid
    if (i < m_width && j < m_height)
    {
        unsigned int index = j * m_width + i;
        m_p_pixel_data[index] = aPixelValue;
    }
    // The 2D index is not valid
    else
    {
        std::stringstream error_message;
        error_message << " in File " << __FILE__ <<
            ", in Function " << __FUNCTION__ <<
            ", at Line " << __LINE__;

        throw std::out_of_range(error_message.str());
    }
}


//---------------------------------------------------
inline void Image::getPixel(unsigned int i,
                            unsigned int j,
                            float& aPixelValue) const
//---------------------------------------------------
{
    // The 2D index is valid
    if (i < m_width && j < m_height)
    {
        unsigned int index = j * m_width + i;
        aPixelValue = m_p_pixel_data[index];
    }
    // The 2D index is not valid
    else
    {
        std::stringstream error_message;
        error_message << " in File " << __FILE__ <<
            ", in Function " << __FUNCTION__ <<
            ", at Line " << __LINE__;

        throw std::out_of_range(error_message.str());
    }
}


//----------------------------------------------------------------
inline std::vector<unsigned char> Image::applyLUT(float vmin, float vmax)
//----------------------------------------------------------------
{
    std::vector<unsigned char> p_pixel_data(3 * m_width * m_height);
    for (unsigned int i = 0; i < m_width * m_height; ++i)
    {
        if (m_p_pixel_data[i] < vmin)
        {
            p_pixel_data[i * 3] = 0;
            p_pixel_data[i * 3 + 1] = 0;
            p_pixel_data[i * 3 + 2] = 0;
        }
        else if (m_p_pixel_data[i] > vmax)
        {
            p_pixel_data[i * 3] = 255;
            p_pixel_data[i * 3 + 1] = 255;
            p_pixel_data[i * 3 + 2] = 255;
        }
        else
        {
            p_pixel_data[i] = round(255.0 * (m_p_pixel_data[i * 3] - vmin) / (vmax - vmin));
            p_pixel_data[i] = round(255.0 * (m_p_pixel_data[i * 3 + 1] - vmin) / (vmax - vmin));
            p_pixel_data[i] = round(255.0 * (m_p_pixel_data[i * 3] - vmin) / (vmax - vmin));
        }
    }
    
    return p_pixel_data;
}
