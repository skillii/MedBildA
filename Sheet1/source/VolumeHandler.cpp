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


#include "VolumeHandler.h"



#define ITK_IMAGE_ALLOCATE_MACRO( NewImageType, newImage, \
  PrototypeImageType, prototypeImage ) \
  NewImageType::Pointer newImage = NewImageType::New(); \
  newImage->SetLargestPossibleRegion( \
  prototypeImage->GetLargestPossibleRegion() ); \
  newImage->SetRequestedRegion( \
  prototypeImage->GetRequestedRegion() ); \
  newImage->SetBufferedRegion( \
  prototypeImage->GetBufferedRegion() ); \
  newImage->SetSpacing( prototypeImage->GetSpacing() ); \
  newImage->SetOrigin( prototypeImage->GetOrigin() ); \
  newImage->SetDirection( prototypeImage->GetDirection() ); \
  newImage->Allocate();

#define ITK_EXCEPTION_CHECKED( message, arg, return_value ) \
  try \
{ \
  if (message != std::string("")) \
  std::cout << message << std::endl; \
  arg; \
  } \
  catch( itk::MemoryAllocationError& err ) \
{ \
  std::cout << "ITK Memory Allocation Error caught!" << std::endl; \
  std::cout << err << std::endl; \
  return return_value; \
  } \
  catch( itk::ExceptionObject& err ) \
{ \
  std::cout << "ITK Exception Object caught!" << std::endl; \
  std::cout << err << std::endl; \
  return return_value; \
  } \
  catch( std::exception& err ) \
{ \
  std::cout << "std Exception Object caught!" << std::endl; \
  std::cout << err.what() << std::endl; \
  return return_value; \
  } \
  catch( ... ) \
{ \
  std::cout << "unknown Exception Object caught!" << std::endl; \
  return return_value; \
  }


////////////////////////////////////////////////////////////////////////////////
VolumeHandler::VolumeHandler()
{
}

////////////////////////////////////////////////////////////////////////////////
VolumeHandler::~VolumeHandler()
{
}


////////////////////////////////////////////////////////////////////////////////
FloatImageType::Pointer VolumeHandler::getHostMem()
{
  return m_host_mem;
}

void VolumeHandler::getOrigMinMaxValues(float& orig_min_value, float& orig_max_value) const
{
  orig_min_value = m_min_orig_gray_value;
  orig_max_value = m_max_orig_gray_value;
}

////////////////////////////////////////////////////////////////////////////////
bool VolumeHandler::readVolume(const std::string& file_name, const bool normalize)
{
  // check if our file extension is one of the supported volume file extensions
  if (!isVolumeFileExtensionSupported(file_name))
    return false;

  VolumeReaderType::Pointer reader = VolumeReaderType::New();
  reader->SetFileName( file_name.c_str() );
  std::string message = std::string("Reading volume ") + file_name;

  // this macro will catch exceptions and return false if one was caught
  ITK_EXCEPTION_CHECKED( message, reader->Update(), false );
  m_host_mem = reader->GetOutput();

  MinMaxCalculator::Pointer calculator= MinMaxCalculator::New();

  calculator->SetImage(m_host_mem);
  calculator->Compute();

  const float max_orig_val = calculator->GetMaximum();
  const float min_orig_val = calculator->GetMinimum();

  if (normalize)
  {
    itk::ImageRegionIterator<FloatImageType> it(m_host_mem, m_host_mem->GetRequestedRegion());
    
	  // normalize the volume
    for(it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
      FloatImageType::ValueType value = it.Get();

      FloatImageType::ValueType normalized_value = (value-min_orig_val) / (max_orig_val-min_orig_val);
      it.Set(normalized_value);
    }
  }

  m_min_orig_gray_value = min_orig_val;
  m_max_orig_gray_value = max_orig_val;

  return true;
}

////////////////////////////////////////////////////////////////////////////////
bool VolumeHandler::setVolume(FloatImageType::Pointer mem, const float min_orig_value, const float max_orig_value)
{
  m_min_orig_gray_value = min_orig_value;
  m_max_orig_gray_value = max_orig_value;

  m_host_mem = mem;
  return true;
}

////////////////////////////////////////////////////////////////////////////////
bool VolumeHandler::writeVolume(const std::string &file_name, const bool unnormalize) const
{
  // check if our file extension is one of the supported volume file extensions
  if (!isVolumeFileExtensionSupported(file_name))
    return false;

  FloatImageType::Pointer image = m_host_mem;

  if (!unnormalize)
  {
    FloatVolumeWriterType::Pointer writer = FloatVolumeWriterType::New();
    writer->SetFileName( file_name.c_str() );
    writer->SetInput( image );
    ITK_EXCEPTION_CHECKED("", writer->Update(), false);
  }
  else
  {
    ITK_IMAGE_ALLOCATE_MACRO(ShortImageType, image_to_write, FloatImageType, image);

    // revert the normalization
    itk::ImageRegionConstIterator<FloatImageType> it(image, image->GetRequestedRegion());
    itk::ImageRegionIterator<ShortImageType> it2(image_to_write,
                                                 image_to_write->GetRequestedRegion());
    for(it.GoToBegin(),it2.GoToBegin(); !it.IsAtEnd(); ++it,++it2)
    {
      FloatImageType::ValueType value = m_min_orig_gray_value +
          it.Get() * (m_max_orig_gray_value - m_min_orig_gray_value);

      it2.Set( static_cast<ShortImageType::ValueType>(value) );
    }

    ShortVolumeWriterType::Pointer writer = ShortVolumeWriterType::New();
    writer->SetFileName( file_name.c_str() );
    writer->SetInput( image_to_write );
    ITK_EXCEPTION_CHECKED("", writer->Update(), false);
  }

  return true;
}

bool VolumeHandler::writeVolumeBinary(const std::string &file_name) const
{
  // check if our file extension is one of the supported volume file extensions
  if (!isVolumeFileExtensionSupported(file_name))
    return false;

  FloatImageType::Pointer image = m_host_mem;

  ITK_IMAGE_ALLOCATE_MACRO(ShortImageType, image_to_write, FloatImageType, image);

  itk::ImageRegionConstIterator<FloatImageType> it(image, image->GetRequestedRegion());
  itk::ImageRegionIterator<ShortImageType> it2(image_to_write,
                                               image_to_write->GetRequestedRegion());
  for(it.GoToBegin(),it2.GoToBegin(); !it.IsAtEnd(); ++it,++it2)
  {
    it2.Set( static_cast<ShortImageType::ValueType>( it.Get() ) );
  }

  ShortVolumeWriterType::Pointer writer = ShortVolumeWriterType::New();
  writer->SetFileName( file_name.c_str() );
  writer->SetInput( image_to_write );
  ITK_EXCEPTION_CHECKED("", writer->Update(), false);

  return true;
}


int VolumeHandler::splitString(const std::string& input,
                               const std::string& delimiter, std::vector<std::string>& results) const
{
  int iPos = 0;
  int newPos = -1;
  int sizeS2 = delimiter.size();
  int isize = input.size();

  std::vector<int> positions;

  newPos = input.find (delimiter, 0);

  if( newPos < 0 ) { return 0; }

  bool skip_first_character = false;

  if (newPos == 0)
  {
    newPos=input.find(delimiter, 1);
    skip_first_character = true;
  }


  int numFound = 0;

  while( newPos > iPos )
  {
    numFound++;
    positions.push_back(newPos);
    iPos = newPos;
    newPos = input.find (delimiter, iPos+sizeS2+1);
  }

  for( unsigned int i=0; i <= positions.size(); i++ )
  {
    std::string s;
    if( i == 0 )
    {
      if (skip_first_character)
        s = input.substr( i+1, positions[i]-1 );
      else
        s = input.substr( i, positions[i] );
    }
    else
    {

      int offset = positions[i-1] + sizeS2;
      if( offset < isize )
      {
        if( i == positions.size() )
        {
          s = input.substr(offset);
        }
        else if( i > 0 )
        {
          s = input.substr( positions[i-1] + sizeS2,
                            positions[i] - positions[i-1] - sizeS2 );
        }
      }
    }


    if( s.size() > 0 )
    {
      results.push_back(s);
    }
  }

  return numFound;
}

bool VolumeHandler::isVolumeFileExtensionSupported(const std::string& file_name) const
{
  std::vector<std::string> validFilenameExtensions;
  validFilenameExtensions.push_back(std::string("hdr")); // Analyze volume data
  validFilenameExtensions.push_back(std::string("mha")); // ITK volume data
  validFilenameExtensions.push_back(std::string("mhd")); // ITK volume data
  validFilenameExtensions.push_back(std::string("nifti")); // Nifti volume data
  validFilenameExtensions.push_back(std::string("vtk")); // Nifti volume data
  validFilenameExtensions.push_back(std::string("nii")); // Nifti volume data

  // extract the extension
  //std::cout << "filename: " << file_name.c_str() << std::endl;

  std::vector<std::string> paths;
  const int numPathSplits = splitString(file_name,std::string("/"),paths);
  // last one is the filename
  const std::string filenameRelative = (numPathSplits > 1) ? paths[paths.size()-1] : file_name;

  std::vector<std::string> fileParts;
  const int numExtensionSplits = splitString(filenameRelative,std::string("."),fileParts);

  if (numExtensionSplits > 0)
  {
    const std::string extension = fileParts[fileParts.size()-1];
    //std::cout << "extension: " << extension << std::endl;

    for (std::vector<std::string>::const_iterator it = validFilenameExtensions.begin();
         it != validFilenameExtensions.end(); ++it)
    {
      if (*it == extension)
      {
        //std::cout << "found valid extension!" << std::endl;
        return true;
      }
    }

    std::cout << "VolumeHandler::isVolumeFilenameValid: " << extension << " is not a supported" <<
                 " ITK volume data extension!" << std::endl;
    return false;
  }
  else
  {
    std::cout << "VolumeHandler::isVolumeFilenameValid: could not find an extension in " <<
                 filenameRelative.c_str() << std::endl;
    return false;
  }
}

