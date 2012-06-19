/*
* Copyright (c) ICG. All rights reserved.
*
* Institute for Computer Graphics and Vision
* Graz University of Technology / Austria
*
*
* This software is distributed WITHOUT ANY WARRANTY; without even
* the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
* PURPOSE.  See the above copyright notices for more information.
*
*
* Project     : Tools
* Module      : Volume Reader
* Class       : $RCSfile$
* Language    : C++/CUDA
* Description :
*
* Author     :
* EMail      :
*
*/


#ifndef VOLUMEHANDLER_H
#define VOLUMEHANDLER_H


#include "Definitions.h"

class VolumeHandler
{
  public:

    ////////////////////////////////////////////////////////////////////////////////
    //! Constructor of VolumeHandler
    ////////////////////////////////////////////////////////////////////////////////
    VolumeHandler();

    ////////////////////////////////////////////////////////////////////////////////
    //! Destructor of VolumeHandler
    ////////////////////////////////////////////////////////////////////////////////
    ~VolumeHandler();

    ////////////////////////////////////////////////////////////////////////////////
    //! Returns a host memory reference to floating point data (might be normalized)
    //! @return ITK pointer to image data
    ////////////////////////////////////////////////////////////////////////////////
    FloatImageType::Pointer getHostMem();

    ////////////////////////////////////////////////////////////////////////////////
    //! Returns the minimum and maximum values of the image, directly after reading
    //! (before the scaling to [0,1])
    ////////////////////////////////////////////////////////////////////////////////
    void getOrigMinMaxValues(float& orig_min_value, float& orig_max_value) const;

    ////////////////////////////////////////////////////////////////////////////////
    //! Reads an ITK volume file
    //! @param file_name     image to be read
    //! @param normalize     bool indicating if the volume should be transformed to
    //!                      gray values between [0,1]
    //! @return              reading was successful
    ////////////////////////////////////////////////////////////////////////////////
    bool readVolume(const std::string &file_name, const bool normalize);

    ////////////////////////////////////////////////////////////////////////////////
    //! Attaches a given host memory volume to the VolumeHandler to prepare for
    //! writing. Assumes volume with gray value range between [0,1]. Needs original 
    //! gray value extremas
    //!
    //! @param mem     image to be written later-on, expected to lie between [0,1]
    //! @return        setting was successful
    ////////////////////////////////////////////////////////////////////////////////
    bool setVolume(FloatImageType::Pointer, const float min_orig_value, const float max_orig_value);

    ////////////////////////////////////////////////////////////////////////////////
    //! Writes a volume file
    //! @param file_name     image to be written
    //! @param unnormalize   bool indicating if the volume should be transformed back
    //!                      from gray values between [0,1] to its original gray values
    //! @return              image was successfully written
    ////////////////////////////////////////////////////////////////////////////////
    bool writeVolume(const std::string &file_name, const bool unnormalize = true) const;

    ////////////////////////////////////////////////////////////////////////////////
    //! Writes a binary volume to a file
    //! @param file_name     image to be written, float labels 0,1 are used for back-
    //!                      and foreground -> writes as signed short volumes with 0,1
    //! @return              image was successfully written
    ////////////////////////////////////////////////////////////////////////////////
    bool writeVolumeBinary(const std::string &file_name) const;

  private:

    // copy constructor and assignment operator not implemented on purpose!
    VolumeHandler(const VolumeHandler& rhs);
    VolumeHandler& operator=(const VolumeHandler& rhs);

    int splitString(const std::string& input, const std::string& delimiter,
       std::vector<std::string>& results) const;

    bool isVolumeFileExtensionSupported(const std::string& file_name) const;

    FloatImageType::Pointer m_host_mem;
    float m_min_orig_gray_value;
    float m_max_orig_gray_value;
};

#endif //VOLUMEHANDLER_H
