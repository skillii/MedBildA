/*
 * CenterlineExtraction.h
 *
 *  Created on: Jun 20, 2012
 *      Author: tom
 */

#ifndef CENTERLINEEXTRACTION_H_
#define CENTERLINEEXTRACTION_H_


#include "VolumeHandler.h"
#include <queue>

class CenterlineExtraction
{

   VolumeHandler volumeHandler;
   FloatImageType::Pointer eigenvector[3];
   FloatImageType::Pointer maxMedialnessOverScales;
   FloatImageType::Pointer centerlineImage;

   float min_vals_eigenvector[3];
   float max_vals_eigenvector[3];
   float min_val_maxMedialness;
   float max_val_maxMedialness;

   static const float threshold_low;
   static const float threshold_high;
   static const int skippixels;


   std::queue<FloatImageType::IndexType> reconnectionQueue;


   void performReconnection(FloatImageType::IndexType index, itk::Vector<float, 3> direction);
   void SetAllPixels(FloatImageType::Pointer img, float value);
   //float getMaxNeighbour(FloatImageType::IndexType index, FloatImageType::IndexType &neighbour);
   float getMaxNeighbour(FloatImageType::IndexType index, itk::Vector<float, 3> &direction);

public:

   CenterlineExtraction()
   {

   }

   int readSavedFiles();

   void performNonMaximaSurpression();
   void performReconnection();

   FloatImageType::Pointer getMedialnessImage()
   {
       return maxMedialnessOverScales;
   }

   FloatImageType::Pointer getCenterlineImage()
   {
       return centerlineImage;
   }

};


#endif /* CENTERLINEEXTRACTION_H_ */
