

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

// itk includes

#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkMinimumMaximumImageCalculator.h"

typedef itk::Image< float, 3 > FloatImageType;
typedef itk::Image< float, 4 > PrimalImageType;

// used in VolumeHandler
typedef itk::Image< signed short, 3 > ShortImageType;
typedef itk::Image< unsigned char, 3 > UCharImageType;
typedef itk::ImageFileReader< FloatImageType > VolumeReaderType;
typedef itk::ImageFileWriter< FloatImageType > FloatVolumeWriterType;
typedef itk::ImageFileWriter< ShortImageType > ShortVolumeWriterType;
typedef itk::ImageFileWriter< UCharImageType > UCharVolumeWriterType;
typedef itk::MinimumMaximumImageCalculator< FloatImageType > MinMaxCalculator;

#endif
