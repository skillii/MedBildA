

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

// itk includes

#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkImageDuplicator.h"

#include "itkMultiResolutionPyramidImageFilter.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkGradientImageFilter.h"



typedef itk::Image< float, 3 > FloatImageType;
typedef itk::Image< float, 5 > HessianImageType;
typedef itk::Image< float, 4 > PrimalImageType;
typedef itk::Image< float, 4 > EigenVImageType;


// used in VolumeHandler
typedef itk::Image< signed short, 3 > ShortImageType;
typedef itk::Image< unsigned char, 3 > UCharImageType;
typedef itk::ImageFileReader< FloatImageType > VolumeReaderType;
typedef itk::ImageFileWriter< FloatImageType > FloatVolumeWriterType;
typedef itk::ImageFileWriter< ShortImageType > ShortVolumeWriterType;
typedef itk::ImageFileWriter< UCharImageType > UCharVolumeWriterType;
typedef itk::MinimumMaximumImageCalculator< FloatImageType > MinMaxCalculator;



typedef itk::MultiResolutionPyramidImageFilter<FloatImageType, FloatImageType> PyramidFilter;
typedef itk::HessianRecursiveGaussianImageFilter<FloatImageType> HessianFilter;
typedef itk::GradientImageFilter< FloatImageType, float> GradientImageFilter;
typedef itk::ImageDuplicator< FloatImageType > DuplicatorType;


//Attention!! Only one or the other!
//#define USE_PRIMAL_DUAL_ROF 1
#define USE_TVL1 1

#endif
