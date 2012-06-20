/*
 * CenterlineExtraction.h
 *
 *  Created on: Jun 20, 2012
 *      Author: tom
 */

#ifndef CENTERLINEEXTRACTION_H_
#define CENTERLINEEXTRACTION_H_


#include "VolumeHandler.h"
//#include <queu>

class CenterlineExtraction
{

   VolumeHandler volumeHandler;
   FloatImageType::Pointer eigenvector[3];
   FloatImageType::Pointer maxMedialnessOverScales;

   float min_vals_eigenvector[3];
   float max_vals_eigenvector[3];
   float min_val_maxMedialness;
   float max_val_maxMedialness;

   static const float threshold_low;
   static const float threshold_high;


   //std::queu

public:

   CenterlineExtraction()
   {

   }

   int readSavedFiles();

   void performNonMaximaSurpression();

   FloatImageType::Pointer getMedialnessImage()
   {
     return maxMedialnessOverScales;
   }

};


#endif /* CENTERLINEEXTRACTION_H_ */
